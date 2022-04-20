#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <sstream>

#include <hera/common.h>
#include <hera/wasserstein.h>
#include <hera/bottleneck.h>

namespace py = pybind11;

using DiagramPoint = hera::DiagramPoint<double>;
using AuctionParams = hera::AuctionParams<double>;
using PairVector = std::vector<std::pair<double, double>>;
using DgmPtVector = std::vector<DiagramPoint>;

void init_ws_geom(py::module&);

double wasserstein_dist_p(const PairVector& points_1, const PairVector& points_2, AuctionParams& params, const std::vector<double>& prices)
{
    return hera::wasserstein_dist(points_1, points_2, params);
}

double wasserstein_cost_p(const PairVector& points_1, const PairVector& points_2, AuctionParams& params, const std::vector<double>& prices)
{
    return hera::wasserstein_cost(points_1, points_2, params);
}

double wasserstein_dist_d(const DgmPtVector& points_1, const DgmPtVector& points_2, AuctionParams& params, const std::vector<double>& prices)
{
    return hera::wasserstein_dist(points_1, points_2, params);
}

double wasserstein_cost_d(const DgmPtVector& points_1, const DgmPtVector& points_2, AuctionParams& params, const std::vector<double>& prices)
{
    return hera::wasserstein_cost(points_1, points_2, params);
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

PYBIND11_MODULE(_hera, m)
{
    m.doc() = "Hera python bindings";

    py::class_<DiagramPoint>(m, "DiagramPoint")
            .def(py::init<double, double, int, int>())
            .def("__repr__", [](const DiagramPoint& p) { std::stringstream ss; ss << p; return ss.str(); })
            .def("get_birth", &DiagramPoint::getRealX)
            .def("get_death", &DiagramPoint::getRealY)
            .def_readwrite("id", &DiagramPoint::id)
            .def_readwrite("user_tag", &DiagramPoint::user_tag)
            ;

    py::class_<AuctionParams>(m, "AuctionParams")
            .def(py::init<>())
            .def_readwrite("wasserstein_power", &AuctionParams::wasserstein_power)
            .def_readwrite("delta", &AuctionParams::delta)
            .def_readwrite("internal_p", &AuctionParams::internal_p)
            .def_readwrite("epsilon_common_ratio", &AuctionParams::epsilon_common_ratio)
            .def_readwrite("inital_epsilon", &AuctionParams::initial_epsilon)
            .def_readwrite("max_num_phases", &AuctionParams::max_num_phases)
            .def_readwrite("max_bids_per_round", &AuctionParams::max_bids_per_round)
            .def_readwrite("dim", &AuctionParams::dim)
            .def_readwrite("final_relative_error", &AuctionParams::final_relative_error)
            .def_readwrite("tolerate_max_iter_exceeded", &AuctionParams::tolerate_max_iter_exceeded)
            .def_readwrite("return_matching", &AuctionParams::return_matching)
            .def_readwrite("match_inf_points", &AuctionParams::match_inf_points)
            .def_readwrite("matching_a_to_b", &AuctionParams::matching_a_to_b_)
            .def_readwrite("matching_b_to_a", &AuctionParams::matching_b_to_a_)
            .def("clear_matching", &AuctionParams::clear_matching)
            .def("add_to_matching", &AuctionParams::add_to_matching)
            ;

    // bottleneck
    m.def("bottleneck_distance_approx", bottleneck_distance_approx);
    m.def("bottleneck_distance_approx_with_edge", bottleneck_distance_approx_with_edge);
    m.def("bottleneck_distance_exact", bottleneck_distance_exact);
    m.def("bottleneck_distance_exact_with_edge", bottleneck_distance_exact_with_edge);

    // Wasserstein
    m.def("wasserstein_dist", wasserstein_dist_p);
    m.def("wasserstein_cost", wasserstein_cost_p);
    m.def("wasserstein_dist", wasserstein_dist_d);
    m.def("wasserstein_cost", wasserstein_cost_d);

    // Wasserstein - point clouds
    init_ws_geom(m);
}
