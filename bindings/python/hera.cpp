#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <sstream>

#include <hera/common.h>
#include <hera/wasserstein.h>
#include <hera/bottleneck.h>

namespace py = pybind11;

void init_ws_geom(py::module&);

double wasserstein_dist(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2, const std::vector<double>& prices)
{
    hera::AuctionParams<double> params;
    return hera::wasserstein_dist(points_1, points_2, params);
}

void init_ws(py::module& m)
{
    m.def("wasserstein_dist", wasserstein_dist);
}

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

PYBIND11_MODULE(_hera, m)
{
    using DiagramPoint = hera::DiagramPoint<double>;

    m.doc() = "Hera python bindings";

    py::class_<DiagramPoint>(m, "DiagramPoint")
            .def(py::init<double, double, int, int>())
            .def("__repr__", [](const DiagramPoint& p) { std::stringstream ss; ss << p; return ss.str(); })
            .def("get_birth", &DiagramPoint::getRealX)
            .def("get_death", &DiagramPoint::getRealY)
            .def_readwrite("id", &DiagramPoint::id)
            .def_readwrite("user_tag", &DiagramPoint::user_tag)
            ;

    init_bt(m);
    init_ws(m);
    init_ws_geom(m);
}
