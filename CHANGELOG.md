# 2.0.0
    - Breaking change: structure of directories changed!
      C++ users can/should use Hera in a more standard way:
      either `#include <hera/bottleneck.h>`,
      `#include <hera/wasserstein.h>`, or `#include <hera/matching_distance.h>,
      depending on what function they need. Function names and signatures
      remain unchanged.
      NB: unfortunately, now it is necessary to add two entries to include_directories:
      `[path_to_hera_root_directory]/include` and `[path_to_hera_root_directory]/extern`.
      There is probably a way to fix this with CMake, this is on the TO-DO list.
    - Code for Wasserstein and bottleneck distance computation was unified.
    - Python bindings are now provided.
    - `DiagramPoint` now has id.
    - Fix: replace tbb::parallel_do with tbb::parallel_for_each.
