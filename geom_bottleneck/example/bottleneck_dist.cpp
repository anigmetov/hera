#include <iomanip>
#include "bottleneck.h"

// any container of pairs of doubles can be used,
// we use vector in this example.

typedef std::vector<std::pair<double, double>>  PairVector;

// estimate initial guess on sampled diagram?
constexpr bool useSamplingHeur = false;
// if diagrams contain fewer points, don't use heuristic
constexpr int heurThreshold = 30000;

int main(int argc, char* argv[])
{
    if (argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " file1 file2 [relative_error]. Without relative_error calculate the exact distance." << std::endl;
        return 1;
    }

    PairVector diagramA, diagramB;
    int decPrecision { 0 };
    if (!geom_bt::readDiagramPointSet(argv[1], diagramA, decPrecision)) {
        std::exit(1);
    }

    if (!geom_bt::readDiagramPointSet(argv[2], diagramB, decPrecision)) {
        std::exit(1);
    }

    double res;
    if (argc >= 4) {
        // the third parameter is epsilon,
        // return approximate distance (faster)
        double approxEpsilon =  atof(argv[3]);
        if (approxEpsilon > 0.0) {
            if (useSamplingHeur && diagramA.size() > heurThreshold && diagramB.size() > heurThreshold) {
                res = geom_bt::bottleneckDistApproxHeur(diagramA, diagramB, approxEpsilon);
            } else {
                res = geom_bt::bottleneckDistApprox(diagramA, diagramB, approxEpsilon);
            }
        } else if (approxEpsilon == 0.0) {
            res = geom_bt::bottleneckDistExact(diagramA, diagramB, decPrecision);
        } else {
            std::cerr << "The third parameter (relative error) must be positive!" << std::endl;
            std::exit(1);
        }
    } else {
        // only filenames have been supplied, return exact distance
        res = geom_bt::bottleneckDistExact(diagramA, diagramB, decPrecision);
    }
    std::cout << std::setprecision(15) << res << std::endl;

    // Alternative could be to construct DiagramPointSet
    // using the constructor with iterators.
    // May be useful if the same diagram is used multiple times
    // to avoid copying data from user's container each time.
    
    //geom_bt::DiagramPointSet dA(diagramA.begin(), diagramA.end());
    //geom_bt::DiagramPointSet dB(diagramB.begin(), diagramB.end());
    //double result1 = geom_bt::bottleneckDistExact(dA, dB);
    //std::cout << std::setprecision(15) << result1 << std::endl;

    return 0;
}
