# 2.0.0
    - Breaking change: structure of directories changed!  C++ users can/should use Hera in a more standard way: either `#include <hera/bottleneck.h>`, `#include <hera/wasserstein.h>`, or `#include <hera/matching_distance.h>, depending on what function(s) they need. Function names and signatures remain unchanged.
    - Code for Wasserstein and bottleneck distance computation was unified.
    - Python bindings are now provided.
    - `DiagramPoint` now has id.
    - Fix: replace tbb::parallel_do with tbb::parallel_for_each.

# 2.0.1
    - Fix: Wasserstein distance computation when one diagram has no finite
    off-diagonal points is computed without invoking 2D auction algorithm
    (previously only worked if there were no off-diagonal points at all).
    - License change: BSD for bottleneck and Wasserstein parts, GPL for the whole
    library and the matching distance part.
