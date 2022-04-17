#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

namespace py = pybind11;

#include <wasserstein_pure_geom.hpp>

using DynamicPointVector = hera::ws::dnn::DynamicPointVector<double>;

using Traits = hera::ws::dnn::DynamicPointTraits<double>;


DynamicPointVector convert_2d_points_to_dnn(const std::vector<std::pair<double, double>>& points)
{
    constexpr int dim = 2;
    Traits traits(dim);

    DynamicPointVector result = traits.container(points.size());

    for(size_t i = 0; i < points.size(); ++i) {
        result[i][0] = points[i].first;
        result[i][1] = points[i].second;
    }

    return result;
}


double wasserstein_cost(const std::vector<std::pair<double, double>>& points_1, const std::vector<std::pair<double, double>>& points_2, const std::vector<double>& prices)
{
    using Traits = hera::ws::dnn::DynamicPointTraits<double>;
    constexpr int dim = 2;
    hera::ws::dnn::DynamicPointTraits<double> traits(dim);

    if (points_1.size() != points_2.size()) {
        std::cerr << "points_1.size = " << points_1.size() << " != points_2.size = " << points_2.size() << std::endl;
        throw std::runtime_error("Point clouds must have same cardinality");
    }

    if (prices.size() != points_1.size() and prices.size() != 0) {
        std::cerr << "prices.size = " << prices.size() << " != points_1.size = " << points_1.size() << std::endl;
        throw std::runtime_error("Price - points size mismatch");
    }

    auto dpoints_1 = convert_2d_points_to_dnn(points_1);
    auto dpoints_2 = convert_2d_points_to_dnn(points_2);

    hera::AuctionParams<double> params;
    params.dim = dim;

    return hera::ws::wasserstein_cost(dpoints_1, dpoints_2, params);
}


void init_ws_geom(py::module& m)
{
    m.def("wasserstein_cost", wasserstein_cost);
}
