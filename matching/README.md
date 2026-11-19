# Matching distance between bifiltrations and 2-persistence modules.

## Accompanying paper
M. Kerber, A. Nigmetov, 
Efficient Approximation of the Matching Distance for 2-parameter persistence.
SoCG 2020.

Bug reports can be sent to "anigmetov EMAIL SIGN lbl DOT gov".

## Dependencies

* Your compiler must support C++14.
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

> bifiltration_phat_like  
> 3  
> \# lines starting with \# are ignored  
> \# vertex A has dimension 0, hence no boundary, its index is 0  
> 0 0 0  
> \# vertex B has index 1  
> 0 0 0  
> \# 1-dimensional simplex {A, B}  
> 1 3 4 0 1  

2. To use from your code.

Here you can compute the matching distance either between bi-filtrations or between  persistence modules.
First, you need to include `#include <hera/matching_distance.h>` Practically every class you need is parameterized by Real type, which should be either float or double.  The header provides two functions called `matching_distance.` 
See `example/module_example.cpp` for additional details.

## Parallel evaluation

To enable TBB evaluation in a C++ application, define `MD_USE_TBB` and link TBB
(for CMake consumers, `TBB::tbb`). Parallel evaluation is supported by
`TraverseStrategy::breadth_first` and `TraverseStrategy::breadth_first_value`.
Set `CalculationParams::n_jobs` to 1 for serial execution, a positive value for
the requested concurrency, or 0 (the default) for automatic concurrency.
Without `MD_USE_TBB`, evaluation remains serial.

`n_hera_calls` counts executed bottleneck calls, including prefetched calls whose
results are discarded or unused after early termination.

## Tests

From the repository root:

```sh
cmake -S . -B build -DHERA_BUILD_TESTS=ON
cmake --build build --target matching_tests
ctest --test-dir build -R '^matching_(module|distance)$' --output-on-failure
```

With CMake 4, add `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` to the configure command
for the bundled Catch2 dependency.

If CMake finds TBB, it also provides the parallel regression target:

```sh
cmake --build build --target matching_tests_tbb
ctest --test-dir build -R '^matching_parallel$' --output-on-failure
```

Use `CMAKE_PREFIX_PATH` when TBB is installed outside the usual search paths.
The existing exhaustive bound checks are registered separately as
`matching_bounds`.

## License

See `licence.txt` in the repository root folder.
