# Matching distance between bifiltrations and 2-persistence modules.

## Accompanying paper
M. Kerber, A. Nigmetov, 
Efficient Approximation of the Matching Distance for 2-parameter persistence.
SoCG 2020.

Bug reports can be sent to "anigmetov EMAIL SIGN lbl DOT gov".

## Dependencies

* Your compiler must support C++11.
* Boost.

## Usage:

1. To use a standalone command-line utility matching_dist:

`./matching_dist -d dimension -e relative_error file1 file2`

If `relative_error` is not specified, the default 0.1 is used.
If `dimension` is not specified, the default value is 0.
Run `./matching_dist` without parameters to see other options.

The output is an approximation of the exact distance (which is assumed to be non-zero).
Precisely: if  *d_exact* is the true distance and *d_approx* is the output, then 

>     | d_exact - d_approx | / d_exact < relative_error.

Files file1 and file2 must contain 1-critical bi-filtrations in a plain text format which is similar to PHAT. The first line of a file must say *bifiltration_phat_like*. The second line contains the total number of simplices *N*. The next *N* lines contain simplices in the format *dim x y boundary*.
* *dim*: the dimension of the simplex
* *x, y*: coordinates of the critical value
* *boundary*: indices of simplices forming the boundary of the current simplex. Indices are separated by space.
* Simplices are indexed starting from 0.

For example, the bi-filtration of a segment with vertices appearing at (0,0) and the 1-segment appearing at (3,4) shall be written as:

>    bifiltration_phat_like
>
>    3
>
>    \# lines starting with \# are ignored
>
>    \# vertex A has dimension 0, hence no boundary, its index is 0
>
>    0 0 0
>
>    \# vertex B has index 1
>
>    0 0 0
>
>    \# 1-dimensional simplex {A, B}
>
>    1 3 4 0 1

2. To use from your code.

Here you can compute the matching distance either between bi-filtrations or between  persistence modules.
First, you need to include `#include "matching_distance.h"` Practically every class you need is parameterized by Real type, which should be either float or double.  The header provides two functions called `matching_distance.` 
See `example/module_example.cpp` for additional details.

## License

See `licence.txt` in the repository root folder.

## Building

CMakeLists.txt can be used to build the command-line utility in the standard
way.  On Linux/Mac/Windows with Cygwin:
>    `mkdir build`
>    `cd build`
>   `cmake ..`
>    `make`

On Windows with Visual Studio: use `cmake-gui` to create the solution in build directory and build it with VS.

The library itself is header-only and does not require separate compilation.
