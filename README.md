# eccons: elliptic curve construction

This program is a toolkit for constructing a few different kinds of elliptic curves.

For technical details, see the comment at the top of `eccons.cc`.

## building, etc.

To build eccons, just type `make`.

You should be able to build eccons with any reasonably modern C++ compiler
(specifically, one that supports C++17). I have tested on several recent
versions of clang and gcc.

You will need the GMP and GMPXX libraries and headers installed in a place your
compiler can find them. On Archlinux, install the `gmp` package; on recent
Debian-alikes (including Ubuntu), install `libgmp10` and `libgmp-dev` packages.

Finally, you will need a reasonably recent version of Sage to generate elliptic
curve parameters. I have tested with versions 8.x and 9.x.

On M1 Mac, you have to run these commands to install gmp before running make:

```
arch -arm64 brew install gmp
find /usr /opt -name "gmp.h"
export LDFLAGS="-I/opt/homebrew/include/ -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib" # Change paths to match the above command's results
```

## usage

The `eccons` executable produces output that can be written to the filesystem
or piped into the `eccons.sage` script.

`eccons --help` gives usage info.

As an example, if we want to generate a complete Edwards curve
having a subgroup of order 2^251 - 9, we might say:

    ./eccons -4 -E 2^251-9 | sage eccons.sage

This should print out the parameters for two such elliptic curves.

## license

    Copyright 2020 Riad S. Wahby

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
