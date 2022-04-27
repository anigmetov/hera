#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <tuple>
#include <utility>
#include <vector>


namespace py = pybind11;

#include <hera/wasserstein_pure_geom.hpp>

using DynamicPointVector = hera::ws::dnn::DynamicPointVector<double>;

using Pair2dPoint = std::pair<double, double>;
using Tuple3dPoint = std::tuple<double, double, double>;

using Vector2dPoints = std::vector<Pair2dPoint>;
using Vector3dPoints = std::vector<Tuple3dPoint>;

using Traits = hera::ws::dnn::DynamicPointTraits<double>;


DynamicPointVector convert_2d_points_to_dnn(const Vector2dPoints& points)
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


DynamicPointVector convert_3d_points_to_dnn(const Vector3dPoints& points)
{
    constexpr int dim = 3;
    Traits traits(dim);

    DynamicPointVector result = traits.container(points.size());

    for(size_t i = 0; i < points.size(); ++i) {
        result[i][0] = std::get<0>(points[i]);
        result[i][1] = std::get<1>(points[i]);
        result[i][2] = std::get<2>(points[i]);
    }

    return result;
}

double wasserstein_cost_geom(const Vector2dPoints& points_1, const Vector2dPoints& points_2, const std::vector<double>& prices)
{
    using Traits = hera::ws::dnn::DynamicPointTraits<double>;
    constexpr int dim = 2;
    hera::ws::dnn::DynamicPointTraits<double> traits(dim);

    if (points_1.size() != points_2.size()) {
        std::cerr << "points_1.size = " << points_1.size() << " != points_2.size = " << points_2.size() << std::endl;
        throw std::runtime_error("Point clouds must have same cardinality");
    }

    auto dpoints_1 = convert_2d_points_to_dnn(points_1);
    auto dpoints_2 = convert_2d_points_to_dnn(points_2);

    hera::AuctionParams<double> params;
    params.dim = dim;

    return hera::ws::wasserstein_cost(dpoints_1, dpoints_2, params, prices);
}


void init_ws_geom(py::module& m)
{
    m.def("wasserstein_cost_geom", wasserstein_cost_geom);
}
