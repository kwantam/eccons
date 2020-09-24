#!/usr/bin/sage
# vim: syntax=python

import heapq
import multiprocessing as mp
import sage.schemes.elliptic_curves.isogeny_small_degree as isd
import sys

# change this value to control max abs(D)
MIN_DVAL = -100000000

# Construct a curve by computing a j-invariant given the CM discriminant.
@parallel
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

def main(DPvalues, COMPL_EDW):
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
            print("p = %d" % curve_info[0])
            print("D = %s" % curve_info[1])
            print("j = %d" % curve_info[2])
            print("  #E  = %d" % curve_info[3])
            print("  #Et = %d" % curve_info[4])
            print("Weierstrass form: y^2 = x^3 + a4 x + a6")
            print("  a4  = %d" % curve_info[5])
            print("  a6  = %d" % curve_info[6])
            if monty_info:
                print("Montgomery form: B y^2 = x^3 + A x^2 + x")
                print("  A   = %d" % monty_info[0])
                print("  B   = %d" % monty_info[1])
                print("Edwards form: a x^2 + y^2 = 1 + d x^2 y^2")
                print("  a   = %d" % monty_info[2])
                print("  d   = %d" % monty_info[3])
                if F(-monty_info[2]).is_square():
                    print("Alt Edwards form: -x^2 + y^2 = 1 + d' x^2 y^2")
                    print("  d'  = %d" % F(-monty_info[3]/monty_info[2]))
                print()
            else:
                print("Cannot convert to Montgomery/Edwards.\n")

if __name__ == "__main__":
    exec(sys.stdin.read())
    if not 'DPvalues' in vars():
        print("Did not find any candidates.")
    else:
        compute_curve.parallel.p_iter.ncpus = NPROC
        main(DPvalues, COMPL_EDW)
