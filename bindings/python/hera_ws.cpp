#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

namespace py = pybind11;

#include <hera/wasserstein.h>


double wasserstein_dist(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2, const std::vector<double>& prices)
{
    hera::AuctionParams<double> params;

    return hera::wasserstein_dist(points_1, points_2, params);
}


void init_ws(py::module& m)
{
    m.def("wasserstein_dist", wasserstein_dist);
}
