## Description

This repository contains software to compute bottleneck and Wasserstein
distances between persistence diagrams, and matching distance
between 2-parameter persistence modules and (1-critical) bi-filtrations.

## License

The software is licensed under
BSD license, see license.txt file.
If you are going to use this software for research purposes,
you probably do not need to worry about that.

## Dependencies

Hera requires C++14 standard support, Boost and TBB.
Catch2, PHAT and pybind11 are included in Hera itself.

## Usage

Hera can be used as a header-only C++ library. It also contains standalone executables and Python bindings.  By default, the standard sequence of commands  

```bash
git clone git@github.com:anigmetov/hera.git
cd hera
mkdir build
cd build
cmake ..
make
```
will not compile any binaries. There are CMake options HERA_BUILD_EXAMPLES and HERA_BUILD_PYTHON_BINDINGS, which are OFF. For tests, there is option HERA_BUILD_TESTS.

- If you want to use python bindings:
```bash
git clone git@github.com:anigmetov/hera.git
cd hera
mkdir build
cd build
cmake .. -DHERA_BUILD_PYTHON_BINDINGS=ON
make -j8
```
The python library `pyhera` will be located in `build/bindings/python`.
You may also need to add `-DPYTHON_EXECUTABLE=path_to_python` as an option
for `cmake`, if you notice that it detects a wrong version of python.

- If you want to run executables from command-line giving them text files as an input:
```bash
git clone git@github.com:anigmetov/hera.git
cd hera
mkdir build
cd build
cmake .. -DHERA_BUILD_EXAMPLES=ON
make -j8
 ```
This will build command-line utilities `bottleneck/bottleneck_dist`,
`wasserstein/wasserstein_dist`, `matching/matching_dist`. See
README files in the corresponding directories.
It will also build `wasserstein/wasserstein_dist_dipha` to read diagrams in
binary DIPHA format and `wasserstein/wasserstein_dist_point_cloud` to
compute Wasserstein distance between sets of points of equal cardinality.

- If you want to use Hera in your code, copy the Hera repository somewhere inside your project. Assuming that it is `your_project/extern/hera`, you only need two lines in your CMakeLists:
```cmake
# to discover hera library
add_subdirectory(extern/hera)

# to add paths that Hera needs to target's include directories
target_link_libraries(your_target PRIVATE hera)
```
Inside the code of `your_target` you `#include<hera/{bottleneck,wasserstein,matching}.h>`, depending on which functions you need.
See examples in [bottleneck_dist.cpp](bottleneck/bottleneck_dist.cpp), [wasserstein_dist.cpp](wasserstein/wasserstein_dist.cpp) and [matching_dist.cpp](matching/matching_dist.cpp).


## References

If you use Hera in your project, we would appreciate if you
cite the corresponding paper.

Bottleneck or Wasserstein distance:

Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov,
"Geometry Helps to Compare Persistence Diagrams.",
Journal of Experimental Algorithmics, vol. 22, 2017, pp. 1--20.
(conference version: ALENEX 2016).
```
@article{jea_hera,
  title={Geometry helps to compare persistence diagrams},
  author={Kerber, Michael and Morozov, Dmitriy and Nigmetov, Arnur},
  journal={Journal of Experimental Algorithmics (JEA)},
  volume={22},
  pages={1--20},
  year={2017},
  publisher={ACM New York, NY, USA}
}
```

Matching distance:

Michael Kerber, Arnur Nigmetov, "Efficient Approximation of the Matching
Distance for 2-parameter persistence.", SoCG 2020

```
@inproceedings{hera_matching,
  title={Efficient Approximation of the Matching Distance for 2-Parameter Persistence},
  author={Kerber, Michael and Nigmetov, Arnur},
  booktitle={36th International Symposium on Computational Geometry: SoCG 2020},
  pages={LIPIcs--SoCG},
  year={2020},
  organization={Schloss Dagstuhl-Leibniz-Zentrum f{\"u}r Informatik GmbH}
}
```

## Changes

09.03.2020. Add matching distance, change directory names.
WARNING: geom_bottleneck -> bottleneck
         geom_matching/wasserstein -> wasserstein

See [CHANGELOG.md](CHANGELOG.md) for all further changes.
