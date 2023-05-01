#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace py = pybind11;

#include <hera/wasserstein_pure_geom.hpp>

using Real = double;

using DynamicPointVector = hera::ws::dnn::DynamicPointVector<Real>;

using Pair2dPoint = std::pair<Real, Real>;
using Tuple3dPoint = std::tuple<Real, Real, Real>;

using Vector2dPoints = std::vector<Pair2dPoint>;
using Vector3dPoints = std::vector<Tuple3dPoint>;

using Params = hera::AuctionParams<Real>;
using Result = hera::AuctionResult<Real>;

using Traits = hera::ws::dnn::DynamicPointTraits<Real>;

using NumpyArray = py::array_t<Real, py::array::c_style | py::array::forcecast>;

DynamicPointVector convert_numpy_to_dnn(NumpyArray pts)
{
    if (pts.ndim() != 2)
        throw std::runtime_error("Dimension mismatch, expected 2D array of shape n_points x dim, one point per row");

    const int dim = static_cast<int>(pts.shape(1));
    const size_t n_points = pts.shape(0);

    Traits traits(dim);
    DynamicPointVector result = traits.container(n_points);

    py::buffer_info pts_buf = pts.request();

    Real* pdata {static_cast<Real*>(pts_buf.ptr)};

    for(size_t i = 0 ; i < n_points ; ++i) {
        for(int d = 0 ; d < dim ; ++d) {
            result[i][d] = *pdata++;
        }
    }

    return result;
}

DynamicPointVector convert_2d_points_to_dnn(const Vector2dPoints& points)
{
    constexpr int dim = 2;
    Traits traits(dim);

    DynamicPointVector result = traits.container(points.size());

    for(size_t i = 0 ; i < points.size() ; ++i) {
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

    for(size_t i = 0 ; i < points.size() ; ++i) {
        result[i][0] = std::get<0>(points[i]);
        result[i][1] = std::get<1>(points[i]);
        result[i][2] = std::get<2>(points[i]);
    }

    return result;
}

//Result wasserstein_cost_geom_detailed(const Vector2dPoints& points_1, const Vector2dPoints& points_2, Params& params, const std::vector<Real>& prices)
//{
//    using Traits = hera::ws::dnn::DynamicPointTraits<Real>;
//    constexpr int dim = 2;
//    hera::ws::dnn::DynamicPointTraits<Real> traits(dim);
//
//    if (points_1.size() != points_2.size()) {
//        std::cerr << "points_1.size = " << points_1.size() << " != points_2.size = " << points_2.size() << std::endl;
//        throw std::runtime_error("Point clouds must have same cardinality");
//    }
//
//    auto dpoints_1 = convert_2d_points_to_dnn(points_1);
//    auto dpoints_2 = convert_2d_points_to_dnn(points_2);
//
//    return hera::ws::wasserstein_cost_detailed(dpoints_1, dpoints_2, params, prices);
//}

Result wasserstein_cost_geom_detailed(const NumpyArray& points_1, const NumpyArray& points_2, Params& params, const std::vector<Real>& prices)
{
    auto dpoints_1 = convert_numpy_to_dnn(points_1);
    auto dpoints_2 = convert_numpy_to_dnn(points_2);

    return hera::ws::wasserstein_cost_detailed(dpoints_1, dpoints_2, params, prices);
}

Real wasserstein_cost_geom(const NumpyArray & points_1, const NumpyArray& points_2, Params& params, const std::vector<Real>& prices)
{
    return wasserstein_cost_geom_detailed(points_1, points_2, params, prices).cost;
}

Real wasserstein_cost_geom_no_params(const NumpyArray& points_1, const NumpyArray& points_2, const std::vector<Real>& prices)
{
    hera::AuctionParams<Real> params;
    params.dim = 2;

    return wasserstein_cost_geom(points_1, points_2, params, prices);
}

void init_ws_geom(py::module& m)
{
    m.def("wasserstein_cost_geom_detailed_", wasserstein_cost_geom_detailed);
    m.def("wasserstein_cost_geom_", wasserstein_cost_geom_no_params);
    m.def("wasserstein_cost_geom_", wasserstein_cost_geom);
}
