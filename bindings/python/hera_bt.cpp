#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <utility>

namespace py = pybind11;

#include <bottleneck.h>

double bottleneck_distance_approx(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2, double delta)
{
    hera::bt::MatchingEdge<double> longest_edge;
    return hera::bottleneckDistApprox(points_1, points_2, delta, longest_edge, false);
}

decltype(auto) bottleneck_distance_approx_with_edge(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2, double delta)
{
    hera::bt::MatchingEdge<double> longest_edge;
    double dist = hera::bottleneckDistApprox(points_1, points_2, delta, longest_edge, true);
    return std::make_pair(dist, longest_edge);
}


double bottleneck_distance_exact(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2)
{
    hera::bt::MatchingEdge<double> longest_edge;
    int dec_precision = 14;
    return hera::bottleneckDistExact(points_1, points_2, dec_precision, longest_edge, false);
}

decltype(auto) bottleneck_distance_exact_with_edge(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2)
{
    hera::bt::MatchingEdge<double> longest_edge;
    int dec_precision = 14;
    double dist = hera::bottleneckDistExact(points_1, points_2, dec_precision, longest_edge, true);
    return std::make_pair(dist, longest_edge);
}


void init_bt(py::module& m)
{
    m.def("bottleneck_distance_approx", bottleneck_distance_approx);
    m.def("bottleneck_distance_approx_with_edge", bottleneck_distance_approx_with_edge);
    m.def("bottleneck_distance_exact", bottleneck_distance_exact);
    m.def("bottleneck_distance_exact_with_edge", bottleneck_distance_exact_with_edge);
}
