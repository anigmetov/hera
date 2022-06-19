# 2.0.0
    - Breaking change: structure of directories changed!  C++ users can/should use Hera in a more standard way: either `#include <hera/bottleneck.h>`, `#include <hera/wasserstein.h>`, or `#include <hera/matching_distance.h>, depending on what function(s) they need. Function names and signatures remain unchanged.
    - Code for Wasserstein and bottleneck distance computation was unified.
    - Python bindings are now provided.
    - `DiagramPoint` now has id.
    - Fix: replace tbb::parallel_do with tbb::parallel_for_each.
