// (C) 2020 Riad S. Wahby <rsw@cs.stanford.edu>
//
// This file is part of eccons.
//
// Licensed under the Apache License, Version 2.0 (see
// LICENSE or https://www.apache.org/licenses/LICENSE-2.0).
// This file may not be copied, modified, or distributed
// except according to those terms.

// eccons: elliptic curve construction toolkit

// on input q, find (p, D) s.t. F_p has an elliptic curve of order q with discriminant D
//
// This is an implementation of Alg 2.2, steps 1-3a plus the primality testing of 3b,
// from Br\"{o}ker and Stevenhagen, "Constructing Elliptic Curves of Prime Order."
//     https://arxiv.org/abs/0712.2022
//
// I have extended this algorithm to support finding curves of order 4*q, q prime by
// specializing the algorithm of Br\"{o}ker and Stevenhagen, "Efficient CM-Constructions
// of Elliptic Curves over Finite Fields." Math. Comp. vol 76 no 260, Oct. 2007.
//     https://www.jstor.org/stable/40234483
//
// The above specialization leverages the fact that the algorithm for finding curves
// of arbitrary order N decays to almost exactly the algorithm for finding curves of
// prime order q when N=4q. Specifically, to find curves for order 4*q, run the
// prime-order algorithm for order q, but for the primality test in step 3b, rather
// than test p = q + 1 +/- x (for x a solution to x^2 + (-D)y^2 = 4*q), one should
// test p' = 4*q + 1 +/- 2*x. Then if p' is prime, then F_p' has a curve of order N
// with j-invariant a root of H_D.
//
// When constructing elliptic curves with specified embedding degree, this program
// uses the Cocks-Pinch method. See Boneh, Rubin, and Silverberg, "Finding composite
// order ordinary elliptic curves using the Cocks-Pinch method." J. Number Theory
// vol 131 issue 5, May 2011.
//     https://www.sciencedirect.com/science/article/pii/S0022314X10001344
// See also Freeman, "Methods for constructing pairing-friendly elliptic curves,"
// 10th Workshop on Elliptic Curve Cryptography, Sept. 2006.
//     http://cacr.uwaterloo.ca/conferences/2006/ecc2006/freeman.pdf
//
// When constructing Montgomery curves from the Weierstrass form, this program uses
// the work of Okeya, Kurumatani, Sakurai, "Elliptic curves with the Montgomery-form
// and their cryptographic applications," PKC, January 2000.
//     https://link.springer.com/chapter/10.1007/978-3-540-46588-1_17
// The Edwards form follows from the birational equivalence given by Bernstein,
// Birkner, Joye, Lange, and Peters, "Twisted Edwards Curves," AFRICACRYPT, June 2008.
//     https://eprint.iacr.org/2008/013.pdf
// The authors also give the completeness condition for twisted Edwards curves.
//
// Edwards curves constructed as described above are not complete. To construct a
// complete curve, this program uses work from Morain, "Edwards curves and CM curves."
//     https://arxiv.org/abs/0904.2243
// The algorithm used here to trace a descending path into the isogeny volcano is
// slightly adapted from the one given by Fouquet and Morain, "Isogeny Volcanoes
// and the SEA Algorithm," ANTS, July 2002.
//     https://link.springer.com/chapter/10.1007/3-540-45455-1_23
//
// Finally, the alternative form of the Edwards curve that this program computes
// (when possible) is the one suggested by Hisil, Wong, Carter, and Dawson,
// "Twisted Edwards curves revisited," ASIACRYPT, Dec 2008.
//     https://link.springer.com/chapter/10.1007/978-3-540-89255-7_20
// It has a fast point addition algorithm when using extended coordinates, which
// the above paper also introduces.

#include <array>
#include <cstdlib>
#include <cstring>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <random>
#include <stack>
#include <tuple>
#include <vector>

unsigned get_two_adicity(const mpz_class &p);

bool maybe_print_entry(const mpz_class &D, mpz_class &p, long two_adicity, bool found) {
    unsigned count = two_adicity ? get_two_adicity(p) : 0;
    if (count < two_adicity) {
        return found;
    }

    if (!found) {       // first one, so print out column header
        gmp_printf("DPvalues = [\n#     =D=        =p=\n");
    }

    gmp_printf("    [%Zd, %Zd],", D.get_mpz_t(), p.get_mpz_t());
    if (two_adicity) {
        printf("    # 2-adicity: %d\n", count);
    } else {
        printf("\n");
    }

    return true;
}

// modular square roots via Tonelli-Shanks
int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p) {
    mpz_class w_c, y_c, y_save_c;
    auto *w = w_c.get_mpz_t();
    auto *y = y_c.get_mpz_t();
    auto *y_save = y_save_c.get_mpz_t();
    unsigned int i = 0, s = 0;

    // if n is a multiple of p, then sqrt is 0
    if (mpz_divisible_p(n, p)) {
        mpz_set_ui(q, 0);
        return 1;
    }

    // if p=3 mod 4, q = n ^ ((p+1) / 4) mod p
    if(mpz_tstbit(p, 1) == 1) {
        mpz_set(q, p);
        mpz_add_ui(q, q, 1);
        mpz_fdiv_q_2exp(q, q, 2);
        mpz_powm(q, n, q, p);
        goto return_principle_root;
    }

    // factor out 2^s from p-1
    mpz_set(q, p);
    mpz_sub_ui(q, q, 1);
    for (s = 0; mpz_tstbit(q, s) == 0; s++);
    mpz_fdiv_q_2exp(q, q, s);

    // Search for a non-residue mod p by picking the first w s.t. (w/p) = -1
    for (mpz_set_ui(w, 2); mpz_legendre(w, p) != -1; mpz_add_ui(w, w, 1));

    // w = w^Q mod p
    mpz_powm(w, w, q, p);
    // y = n^Q mod p
    mpz_powm(y, n, q, p);
    // q = n^((Q+1)/2) mod p
    mpz_add_ui(q, q, 1);
    mpz_fdiv_q_2exp(q, q, 1);
    mpz_powm(q, n, q, p);

    // Tonelli-Shanks main loop
    while (true) {
        // find the order of y mod p
        // if y = 1 mod p, then we've found the sqrt
        mpz_set(y_save, y);
        for (i = 0; (i < s) && (mpz_cmp_ui(y, 1) != 0); i++, mpz_powm_ui(y, y, 2, p));
        if (i == 0) {
            goto return_principle_root;
        } else if (i == s) {
            return 0;
        }

        // otherwise, let t = w^(2^(s-i-1))
        mpz_powm_ui(w, w, 1 << (s-i-1), p);

        // update values for next iteration
        s = i;
        // q = q * w
        mpz_mul(q, q, w);
        mpz_mod(q, q, p);
        // w = w^2
        mpz_mul(w, w, w);
        mpz_mod(w, w, p);
        // y = y w
        mpz_mul(y, y_save, w);
        mpz_mod(y, y, p);
    }

return_principle_root:
    // let q be the principal root
    mpz_fdiv_q_2exp(y, p, 1);
    if (mpz_cmp(q, y) == 1) {
        mpz_sub(q, p, q);
    }
    // check that q is indeed the sqrt
    mpz_mul(y, q, q);
    mpz_mod(y, y, p);
    mpz_mod(w, n, p);
    if (mpz_cmp(y, w) != 0) {
        return 0;
    }
    return 1;
}

// we return x s.t. x^2 + |D|y^2 = 4q
//
// This follows section 1.5.3 of Cohen, "A Graduate Course in Computational Number Theory"
int cornacchia(mpz_class &x, const mpz_class &q, const mpz_class &D, const mpz_class &sqD) {
    // variables
    mpz_class tmp;
    mpz_class a;
    mpz_class b;
    mpz_class l;
    auto *mtmp = tmp.get_mpz_t();
    auto *ml = l.get_mpz_t();
    const auto *mD = D.get_mpz_t();

    // set up temporary variables
    a = 2 * q;                          // r_0 = N = 2 * q
    b = sqD;                            // r_1 = sqrt(D) - NOTE: main() ensures that sqD is sqrt(D) mod 4q
    l = 4 * q;
    mpz_sqrt(ml, ml);                   // floor(sqrt(4*q))

    // Euclidean algorithm
    while (b > l) {
        tmp = a % b;
        a = b;
        b = tmp;
    }

    // primitive solutions have (4q - r_k^2) divisible by D and square
    tmp = (4 * q) - (b * b);        // 4*q - r_k^2
    mpz_abs(ml, mD);
    if (mpz_divisible_p(mtmp, ml)) {
        mpz_divexact(mtmp, mtmp, ml);
        if (mpz_perfect_square_p(mtmp)) {
            x = b;
            return 1;
        }
    }

    return 0;
}

void convert_next_value(char *current, unsigned &len, char *&next, char &sgn) {
    // find first non-numeric character
    for (; (len > 0) && (current[0] >= '0') && (current[0] <= '9'); len--, current++);
    sgn = current[0];
    // the below is safe even when len == 0 because len was calculated by strlen
    current[0] = '\0';
    next = current;
}

bool convert_formula (char *current, mpz_class &q) {
    unsigned len = std::strlen(current);
    char *next = nullptr;
    long tmp_l = 0;
    bool op_is_add = true;
    char sgn = '+';
    mpz_class tmp_z;
    mpz_class one = 1;

    q = 0;
    while (len > 0) {
        // grab either 2^x or x
        tmp_z = 0;
        if ((len > 1) && (std::strncmp("2^", current, 2) == 0)) {
            // found 2^, now get the exponent
            current += 2;
            len -= 2;
            if (len == 0) {
                return false;
            }
            convert_next_value(current, len, next, sgn);
            tmp_l = std::strtol(current, nullptr, 10);
            next[0] = sgn;      // (for error reporting)
            current = next;

            // add or subtract 2^exponent
            mpz_mul_2exp(tmp_z.get_mpz_t(), one.get_mpz_t(), tmp_l);
        } else {
            convert_next_value(current, len, next, sgn);
            tmp_z.set_str(current, 0);
            next[0] = sgn;      // (for error reporting)
            current = next;
        }

        // do the operation
        if (op_is_add) {
            q += tmp_z;
        } else {
            q -= tmp_z;
        }

        // grab the next sign, if there is one
        if (len > 0) {
            if (sgn == '+') {
                op_is_add = true;
            } else if (sgn == '-') {
                op_is_add = false;
            } else {
                return false;
            }
            current += 1;
            len -= 1;
        }
    }
    return true;
}

unsigned get_two_adicity(const mpz_class &p) {
    unsigned count = 0;
    mpz_class pm1 = p - 1;
    auto *pm1t = pm1.get_mpz_t();
    while (mpz_divisible_ui_p(pm1t, 2)) {
        count++;
        mpz_divexact_ui(pm1t, pm1t, 2);
    }

    return count;
}

void print_usage (const char *name) {
    std::cout << "Usage:" << std::endl
              << "  " << name << " <q>          [options]" << std::endl
              << R"EOF(
You can specify the target either as an integer or as a formula, e.g., '2^255-19'.

  option      description                                       default 
  --          --                                                --
  -r <max_r>  sets how large a D to search for                  6

  -D <maxD>   Do not consider discriminants |D| > maxD.         ((r+1)*log q)^2

  -k <deg>    Search for a curve with embedding degree <deg>    unconstrained
              via Cocks-Pinch method. Cannot be used with -4.

  -R <seed>   Use <seed> to search for candidate trace values   random
              when -k is set.

  -4          Search for a curve of order 4*q rather than q.    false
              Cannot be used with -k.

  -E          Construct complete twisted Edwards curves only.   false
              Implies -s; should be used with -4 or -k.

  -2 <min>    Print 2-adicity of each identified p; reject      unconstrained
              any p with less than <min> 2-adicity.

  -s          after identifying appropriate D and p,            false
              output a Sage script to construct curves.

  -j <nproc>  Use at most <nproc> parallel processes.           # of cpus

Example:
  )EOF"
              << name << " 2^251-9 -r 10 -s -4" << std::endl;
}

void check_nargs(int nargs, int expect, const char *name, const char* argname) {
    if (nargs < expect) {
        print_usage(name);
        std::cout << std::endl << argname << " requires an argument!" << std::endl;
        std::exit(-1);
    }
}

int main (int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        std::exit(-1);
    }

    // STEP 0: parse argumnts
    bool emit_script = false;
    bool times_four = false;
    bool compl_edw = false;
    int arg_curr = 1;
    long dmax_user = 0;
    long k_target = 0;
    long r = 6;
    long two_adicity = 0;
    long rseed_user = 0;
    long max_nproc = 0;
    mpz_class q;
    while (arg_curr < argc) {
        if (std::strncmp("-R", argv[arg_curr], 2) == 0) {
            check_nargs(argc - arg_curr, 2, argv[0], "-R");
            rseed_user = std::strtol(argv[arg_curr+1], nullptr, 10);
            arg_curr += 2;
        } else if (std::strncmp("-D", argv[arg_curr], 2) == 0) {
            check_nargs(argc - arg_curr, 2, argv[0], "-D");
            dmax_user = std::strtol(argv[arg_curr+1], nullptr, 10);
            arg_curr += 2;
        } else if (std::strncmp("-k", argv[arg_curr], 2) == 0) {
            check_nargs(argc - arg_curr, 2, argv[0], "-k");
            k_target = std::strtol(argv[arg_curr+1], nullptr, 10);
            arg_curr += 2;
        } else if (std::strncmp("-4", argv[arg_curr], 2) == 0) {
            times_four = true;
            arg_curr += 1;
        } else if (std::strncmp("-E", argv[arg_curr], 2) == 0) {
            compl_edw = true;
            emit_script = true;
            arg_curr += 1;
        } else if (std::strncmp("-2", argv[arg_curr], 2) == 0) {
            check_nargs(argc - arg_curr, 2, argv[0], "-2");
            two_adicity = std::strtol(argv[arg_curr+1], nullptr, 10);
            arg_curr += 2;
        } else if (std::strncmp("-s", argv[arg_curr], 2) == 0) {
            emit_script = true;
            arg_curr += 1;
        } else if (std::strncmp("-r", argv[arg_curr], 2) == 0) {
            check_nargs(argc - arg_curr, 2, argv[0], "-r");
            r = std::strtol(argv[arg_curr+1], nullptr, 10);
            arg_curr += 2;
        } else if (std::strncmp("-j", argv[arg_curr], 2) == 0) {
            check_nargs(argc - arg_curr, 2, argv[0], "-j");
            max_nproc = std::strtol(argv[arg_curr+1], nullptr, 10);
            arg_curr += 2;
        } else if (argv[arg_curr][0] == '-') {
            print_usage(argv[0]);
            std::cout << std::endl << "Unrecognized option " << argv[arg_curr] << std::endl;
            std::exit(-1);
        } else {
            if (!convert_formula(argv[arg_curr], q)) {
                print_usage(argv[0]);
                std::cout << std::endl << "Failed to convert '" << argv[arg_curr] << "' to an integer." << std::endl;
                std::exit(-1);
            }
            arg_curr += 1;
        }
    }

    // sanity checks and preparatory stuff
    q = abs(q);
    k_target = std::abs(k_target);
    r = std::abs(r);
    if (rseed_user) {
        rseed_user = std::abs(rseed_user);
    } else if (k_target != 0) {
        std::random_device rd;
        rseed_user = rd();
    }
    dmax_user = std::abs(dmax_user);
    if (q <= 3 || !mpz_probab_prime_p(q.get_mpz_t(), 128)) {
        print_usage(argv[0]);
        std::cout << std::endl << "This program requires q > 3 and prime." << std::endl;
        std::exit(-1);
    } else if (k_target != 0) {
        constexpr const char *const ERR1 = "This program does not support -4 and -k at the same time, but often the\nresults from -k will have cofactor divisible by four (hint: try -R).";
        constexpr const char *const ERR2 = "This program requires q to be congruent to 1 mod k.";
        const char *err = nullptr;
        if (times_four) {
            err = ERR1;
        } else if (q % k_target != 1) {
            err = ERR2;
        }
        if (err != nullptr) {
            print_usage(argv[0]);
            std::cout << std::endl << err << std::endl;
            std::exit(-1);
        }
    }
    mpz_class N = q * (times_four ? 4 : 1);
    gmp_printf("N = %Zd\n", N.get_mpz_t());

    // now begins the algorithm
    using DTableEnt = std::array<mpz_class, 2>;
    using DTable = std::vector<DTableEnt>;
    DTable dtab{};

    unsigned log_q = mpz_sizeinbase(q.get_mpz_t(), 2);
    mpz_class p{3};
    mpz_class tmp_z;
    mpz_class tmp_z2;
    mpz_class tmp_z3;
    mpz_class kq_div;
    if (k_target != 0) {
        kq_div = (q - 1) / k_target;
        printf("# random seed was %ld\n", rseed_user);
    } else {
        printf("\n");
    }
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(rseed_user);

    // STEP 1: Create a table of (D, sqrt(D)) candidates
    while (p < (r+1) * log_q) {
        if (mpz_legendre(q.get_mpz_t(), p.get_mpz_t()) == 1) {
            if ((p-1) % 4 == 0) {
                tmp_z2 = p;
            } else {
                tmp_z2 = -p;
            }

            if (mpz_sqrtm(tmp_z.get_mpz_t(), tmp_z2.get_mpz_t(), q.get_mpz_t()) != 1) {
                gmp_printf("# ERROR detected computing sqrt(%Zd) mod %Zd\n", tmp_z2.get_mpz_t(), q.get_mpz_t());
                std::exit(-1);
            } else {
                dtab.emplace_back(DTableEnt{tmp_z2, tmp_z});
            }
        }

        // go to the next prime
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    }

    // STEP 2: Compute products of the table elements to find D and sqrt(D)
    // max D we will consider is ((r+1)log_q)^2
    mpz_class x;
    mpz_class Dmax;
    if (dmax_user == 0) {
        Dmax = (r+1) * log_q;
        Dmax *= Dmax;
    } else {
        Dmax = dmax_user;
    }

    // walk the D candidates, finding candidates
    using DStackEnt = std::tuple<size_t, DTableEnt>;
    using DStack = std::stack<DStackEnt>;
    DStack dstk{};

    dstk.push(DStackEnt{0, {1, 1}});
    auto last_ent = dtab.size();
    bool found = false;

    while (!dstk.empty()) {
        DStackEnt curr = std::move(dstk.top());
        dstk.pop();
        auto next = std::get<0>(curr);
        auto& Dent = std::get<1>(curr);

        while (true) {
            // if we've run out of dtab, give up
            if (next >= last_ent) {
                break;
            }

            // if we've gotten too big, give up
            tmp_z = Dent[0] * dtab[next][0];
            tmp_z = abs(tmp_z);
            if (tmp_z > Dmax) {
                break;
            }

            // if there are two paths, push one onto the stack
            if (next + 1 < last_ent) {
                dstk.push(DStackEnt{next+1, Dent});
            }

            // extend current path
            Dent[0] *= dtab[next][0];
            Dent[1] = (Dent[1] * dtab[next][1]) % q;
            next++;

            // for order q or 4q, unconstrained k, check that candidate D is -3 mod 8
            if ((k_target == 0) && (Dent[0] % 8 == -3)) {
                // make sure that we have the odd sqrt of D s.t. it's also a sqrt mod 4q
                if (mpz_even_p(Dent[1].get_mpz_t())) {
                    Dent[1] = q - Dent[1];
                }

                // STEP 3: run Cornacchia's algorithm; if q + 1 +/- x is prime, we've found p
                //         If we're looking for a curve of order 4q, check 4q + 1 +/- 2x instead.
                if (cornacchia(x, q, Dent[0], Dent[1])) {
                    for (unsigned i = 0; i < 2; i++) {
                        p = N + 1 + (times_four ? 2 : 1) * (i ? x : -x);
                        if (mpz_probab_prime_p(p.get_mpz_t(), 128) > 0) {
                            found = maybe_print_entry(Dent[0], p, two_adicity, found);
                        }
                    }
                }
            // for constrained k, just need D < 0 and =0 or =1 mod 4
            } else if ((k_target != 0) && (Dent[0] < 0) && ((Dent[0] % 4 == 0) || Dent[0] % 4 == -3)) {
                unsigned maxrep = log_q * log_q;
                for (unsigned i = 0; i < maxrep; i++) {
                    mpz_invert(tmp_z3.get_mpz_t(), Dent[1].get_mpz_t(), q.get_mpz_t());
                    do {
                        // X of order k mod q
                        tmp_z2 = rand.get_z_range(q);
                        mpz_powm_ui(tmp_z.get_mpz_t(), tmp_z2.get_mpz_t(), k_target, q.get_mpz_t());
                    } while (tmp_z == 1);
                    for (unsigned j = 0; j < 2; j++) {
                        // Y = (X-1) * sqrt(-D) mod q
                        tmp_z2 = (tmp_z - 1);
                        if (j) {
                            tmp_z2 *= -tmp_z3;
                        } else {
                            tmp_z2 *= tmp_z3;
                        }
                        tmp_z2 %= q;

                        // p = ((X+1)^2 + DY^2)
                        p = tmp_z + 1;
                        p *= p;
                        tmp_z2 *= tmp_z2;
                        tmp_z2 *= -Dent[0];
                        p += tmp_z2;

                        // make sure p is divisible by 4
                        if (!mpz_divisible_ui_p(p.get_mpz_t(), 4)) {
                            continue;
                        }
                        mpz_divexact_ui(p.get_mpz_t(), p.get_mpz_t(), 4);

                        if (mpz_probab_prime_p(p.get_mpz_t(), 128) > 0) {
                            found = maybe_print_entry(Dent[0], p, two_adicity, found);
                            if (found) {
                                // break two levels
                                i = maxrep;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // STEP 4: Print out info for the user.
    if (!found) {
        std::cout << "# Sorry! I didn't find any candidates. Try running with higher -r." << std::endl;
    } else {
        std::cout << ']' << std::endl;
        if (emit_script) {
            std::cout << std::endl
                      << R"EOF(
import heapq
import multiprocessing as mp
import sage.schemes.elliptic_curves.isogeny_small_degree as isd

MIN_DVAL = -100000000 ### change this value to control max abs(D)
NPROC = )EOF";
            if (max_nproc) {
                std::cout << max_nproc;
            } else {
                std::cout << "mp.cpu_count()";
            }
            std::cout << " ### change this value to increase or decrease number of parallel threads\nCOMPL_EDW = "
                      << (compl_edw ? "True" : "False")
                      << R"EOF( ### when True, generate only complete twisted Edwards curves

# Construct a curve by computing a j-invariant given the CM discriminant.
@parallel(ncpus=NPROC)
def compute_curve(dpj):
    (D, p) = dpj[:2]
    F = GF(p)
    if len(dpj) == 3:
        j = dpj[2]
    else:
        j = hilbert_class_polynomial(D).any_root(F)
    (E, npoints) = maybe_twist(EllipticCurve(F, j=j), D, p)
    if E and COMPL_EDW:
        E = to_complete_edwards(E, D, p)
    return (E, npoints)

# decide whether to use the curve or the twist by testing for correct point order
MT_NUM_REPS = 16
def maybe_twist(E, D, p):
    npoints = get_npoints(D, p)
    if npoints is None:
        trace = E.trace_of_frobenius()
        e_np = p + 1 - trace
        et_np = p + 1 + trace
        use_curve = False
        use_twist = False
        if e_np % N == 0:
            use_curve = True
            npoints = e_np
        elif et_np % N == 0:
            use_twist = True
            npoints = et_np
    else:
        Et = E.quadratic_twist()
        zero = E(0, 1, 0)
        zeroT = Et(0, 1, 0)
        use_curve = all( zero == npoints * E.random_point() for _ in range(0, MT_NUM_REPS) )
        use_twist = all( zeroT == npoints * Et.random_point() for _ in range(0, MT_NUM_REPS) )
    if use_curve and use_twist:
        print("ERROR: supersingular curve?")
        return (None, None)
    elif use_curve:
        return (E, npoints)
    elif use_twist:
        return (Et, npoints)
    print("ERROR: neither curve nor twist have expected order")
    return (None, None)

# convert a twisted Edwards curve to a complete curve via 2-isogenies.
# This is from Morain, F. "Edwards curves and CM curves." arXiv 0904.2243, \S 4.2.
def to_complete_edwards(E, D, p):
    M = to_montgomery(E)
    if not M:
        return None
    if is_complete(E, M):
        return E
    nsteps = num_levels(D, p)
    if nsteps is None:
        return None
    step = 0
    ctr = 1
    q = [(0, 0, E)]
    new_E = None
    while q:
        (prio, _, curve) = heapq.heappop(q)
        if prio <= -nsteps:
            if is_complete(curve):
                new_E = curve
                break
            if prio < -nsteps:
                continue
        for new_iso in set( iso.codomain() for iso in isd.isogenies_2(curve) ):
            heapq.heappush(q, (prio - 1, ctr, new_iso))
            ctr += 1
    if new_E is None:
        print("ERROR: no complete twisted Edwards curves found for D=%d" % D)
    return new_E

# This is from K. Okeya, H. Kurumatani, and K. Sakurai. "Elliptic Curves with
# the Montgomery-Form and Their Cryptographic Applications." PKC, January 2000.
def to_montgomery(E):
    a = E.a4()
    b = E.a6()
    F = E.base_field()
    p = F.order()
    R.<x> = F[]
    poly = x^3 + a*x + b
    salvals = [ (1/F(3*alpha^2+a).sqrt(),alpha) for (alpha, _) in poly.roots()[:1] if F(3*alpha^2+a).is_square() ]
    if not salvals:
        return None
    (s, alpha) = salvals[0]
    (A, B) = (F(3*alpha*s), F(s))
    (aP, d) = (F((A+2)/B), F((A-2)/B))
    return [A, B, aP, d]

# solve 4p = U^2 - DV^2
R.<u,v> = PolynomialRing(ZZ)
def solve_uv(D, p):
    solns = solve_diophantine(u**2-D*v**2-4*p)
    if not solns:
        print("ERROR: expected solutions to U^2-DV^2=4p, but found none\n")
        return None
    return (int(abs(solns[0][0])), int(abs(solns[0][1])))

# given CM discriminant D, can easily solve for trace
def get_npoints(D, p):
    solns = solve_uv(D, p)
    if not solns:
        return None
    (U, _) = solns
    npoints = p + 1 - U
    if npoints % N != 0:
        npoints = p + 1 + U
    if npoints % N != 0:
        return None # fall back to point counting
    return npoints

# How many levels of the isogeny volcano do we need to descend to find a complete curve?
# This is from Morain, F. "Edwards curves and CM curves." arXiv 0904.2243, \S 4.2
def num_levels(D, p):
    solns = solve_uv(D, p)
    if not solns:
        return None
    (_, V) = solns
    nsteps = 0
    while V % 2 == 0:
        nsteps += 1
        V /= 2
    if nsteps == 0:
        print("ERROR: expected V to be even, but it was odd\n")
        return None
    return nsteps

# An Edwards curve is complete if a is square and d is nonsquare.
# This is from Bernstein, D.J., Birkner, P., Joye, M., Lange, T., and Peters, C.
# "Twisted Edwards Curves." Proc. AFRICACRYPT 2008.
def is_complete(E, M=None):
    if M is None:
        M = to_montgomery(E)
    if M is None:
        return None
    F = E.base_field()
    if F(M[2]).is_square() and not F(M[3]).is_square():
        return True
    return False

for (curve_input, (E, npoints)) in compute_curve( x for x in DPvalues if x[0] >= MIN_DVAL ):
    if E is None:
        continue
    F = E.base_field()
    D = curve_input[0][0][0]
    p = curve_input[0][0][1]
    assert p == F.order()
    curve_info = [p, D, E.j_invariant(), npoints, 2 * p + 1 - npoints, E.a4(), E.a6()]
    monty_info = to_montgomery(E)
    if not COMPL_EDW or monty_info:
        print("p = %d\nD = %s\nj = %d\n  #E  = %d\n  #Et = %d\nWeierstrass form: y^2 = x^3 + a4 x + a6\n  a4  = %d\n  a6  = %d" % tuple(curve_info))
        if monty_info:
            print("Montgomery form: B y^2 = x^3 + A x^2 + x\n  A   = %d\n  B   = %d\nEdwards form: a x^2 + y^2 = 1 + d x^2 y^2\n  a   = %d\n  d   = %d" % tuple(monty_info))
            if F(-monty_info[2]).is_square():
                print("Alt Edwards form: -x^2 + y^2 = 1 + d' x^2 y^2\n  d'  = %d\n" % F(-monty_info[3]/monty_info[2]))
            else:
                print()
        else:
            print("Cannot convert to Montgomery/Edwards.\n")
)EOF";
        }
    }
}
