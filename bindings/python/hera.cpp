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
using AuctionResult = hera::AuctionResult<double>;
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

    py::class_<AuctionResult>(m, "WassersteinResult")
            .def(py::init<>())
            .def_readonly("num_rounds", &AuctionResult::num_rounds)
            .def_readonly("num_phases", &AuctionResult::num_phases)
            .def_readonly("distance", &AuctionResult::distance)
            .def_readonly("cost", &AuctionResult::cost)
            .def_readonly("start_epsilon", &AuctionResult::start_epsilon)
            .def_readonly("final_epsilon", &AuctionResult::final_epsilon)
            .def_readonly("final_relative_error", &AuctionResult::final_relative_error)
            .def_readonly("matching_a_to_b", &AuctionResult::matching_a_to_b_)
            .def_readonly("matching_b_to_a", &AuctionResult::matching_b_to_a_)
            .def_readonly("prices", &AuctionResult::prices)
            .def("compute_distance", &AuctionResult::compute_distance)
            .def("clear_matching", &AuctionResult::clear_matching)
            .def("add_to_matching", &AuctionResult::add_to_matching)
            .def("__str__", [](const AuctionResult& r) { std::stringstream ss; ss << r; return ss.str(); })
            .def("__repr__", [](const AuctionResult& r) { std::stringstream ss; ss << r; return ss.str(); })
            ;

    py::class_<AuctionParams>(m, "WassersteinParams")
            .def(py::init<>())
            .def_readwrite("wasserstein_power", &AuctionParams::wasserstein_power)
            .def_readwrite("delta", &AuctionParams::delta)
            .def_readwrite("internal_p", &AuctionParams::internal_p)
            .def_readwrite("epsilon_common_ratio", &AuctionParams::epsilon_common_ratio)
            .def_readwrite("initial_epsilon", &AuctionParams::initial_epsilon)
            .def_readwrite("max_num_phases", &AuctionParams::max_num_phases)
            .def_readwrite("max_bids_per_round", &AuctionParams::max_bids_per_round)
            .def_readwrite("dim", &AuctionParams::dim)
            .def_readwrite("tolerate_max_iter_exceeded", &AuctionParams::tolerate_max_iter_exceeded)
            .def_readwrite("return_matching", &AuctionParams::return_matching)
            .def_readwrite("match_inf_points", &AuctionParams::match_inf_points)
            .def("__repr__", [](const AuctionParams& p) { std::stringstream ss; ss << p; return ss.str(); })
            .def(py::pickle(
                // __getstate__
                [](const AuctionParams& p) { return py::make_tuple(p.wasserstein_power, p.delta, p.internal_p,
                    p.initial_epsilon, p.epsilon_common_ratio, p.max_num_phases, p.max_bids_per_round,
                    p.dim, p.tolerate_max_iter_exceeded, p.return_matching, p.match_inf_points); },
                // __setstate__
                [](py::tuple t) {
                    if (t.size() != 11)
                        throw std::runtime_error("Invalid tuple for AuctionParams");

                    AuctionParams p;

                    p.wasserstein_power          = t[0].cast<decltype(p.wasserstein_power)>();
                    p.delta                      = t[1].cast<decltype(p.delta)>();
                    p.internal_p                 = t[2].cast<decltype(p.internal_p)>();
                    p.initial_epsilon            = t[3].cast<decltype(p.initial_epsilon)>();
                    p.epsilon_common_ratio       = t[4].cast<decltype(p.epsilon_common_ratio)>();
                    p.max_num_phases             = t[5].cast<decltype(p.max_num_phases)>();
                    p.max_bids_per_round         = t[6].cast<decltype(p.max_bids_per_round)>();
                    p.dim                        = t[7].cast<decltype(p.dim)>();
                    p.tolerate_max_iter_exceeded = t[8].cast<decltype(p.tolerate_max_iter_exceeded)>();
                    p.return_matching            = t[9].cast<decltype(p.return_matching)>();
                    p.match_inf_points           = t[10].cast<decltype(p.match_inf_points)>();

                    return p;
                }))
            ;

    // bottleneck
    m.def("bottleneck_distance_approx", bottleneck_distance_approx);
    m.def("bottleneck_distance_approx_with_edge", bottleneck_distance_approx_with_edge);
    m.def("bottleneck_distance_exact", bottleneck_distance_exact);
    m.def("bottleneck_distance_exact_with_edge", bottleneck_distance_exact_with_edge);

    // Wasserstein
    m.def("wasserstein_dist_", wasserstein_dist_p);
    m.def("wasserstein_cost_", wasserstein_cost_p);
    m.def("wasserstein_dist_", wasserstein_dist_d);
    m.def("wasserstein_cost_", wasserstein_cost_d);

    // Wasserstein - point clouds
    init_ws_geom(m);
}
