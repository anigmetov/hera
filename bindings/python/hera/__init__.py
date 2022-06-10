#!/usr/bin/env python3

import typing
from icecream import ic

# from ._hera import bottleneck_distance_exact_with_edge, bottleneck_distance_exact, bottleneck_distance_approx_with_edge, \
#     bottleneck_distance_approx, WassersteinParams, WassersteinResult, DiagramPoint


# from ._hera import wasserstein_cost_geom_detailed as hera_wasserstein_cost_geom_detailed
# from ._hera import wasserstein_cost_geom as hera_wasserstein_cost_geom

from ._hera import *

def bottleneck_dist(dgm_a, dgm_b, delta: float = 0.01, return_bottleneck_edge: bool = False):
    if delta == 0.0:
        if return_bottleneck_edge:
            return bottleneck_distance_exact_with_edge(dgm_a, dgm_b)
        else:
            return bottleneck_distance_exact(dgm_a, dgm_b)
    elif delta > 0.0:
        if return_bottleneck_edge:
            return bottleneck_distance_approx_with_edge(dgm_a, dgm_b, delta)
        else:
            return bottleneck_distance_approx(dgm_a, dgm_b, delta)
    else:
        raise RuntimeError(f"Relative error delta must be non-negative, got {delta}")


# def wasserstein_cost(pts_1, pts_2, prices):
#     return hera_wasserstein_cost_geom(pts_1, pts_2, prices)
#
#
# def wasserstein_cost_g(pts_1, pts_2, prices):
#     return hera_wasserstein_cost_geom(pts_1, pts_2, params, prices)
#

def wasserstein_cost_geom(pts_a: typing.List[typing.Tuple[float, float]], pts_b: typing.List[typing.Tuple[float, float]],
                     params: typing.Optional[WassersteinParams] = None, prices=None, return_detailed=False):

    if type(pts_a) is not list:
        pts_a = [ (x, y) for x, y in pts_a ]
    if type(pts_b) is not list:
        pts_b = [ (x, y) for x, y in pts_b ]
    if params is None:
        params = WassersteinParams()
    if prices is None:
        prices = []
    if type(prices) is not list:
        prices = list(prices)
    if params.return_matching or return_detailed:
        return wasserstein_cost_geom_detailed_(pts_a, pts_b, params, prices)
    else:
        return wasserstein_cost_geom_(pts_a, pts_b, params, prices)


def wasserstein_dist(dgm_a: typing.List[typing.Tuple[float, float]], dgm_b: typing.List[typing.Tuple[float, float]],
                     params: typing.Optional[WassersteinParams] = None, prices=None):
    if params is None:
        params = WassersteinParams()
    if prices is None:
        prices = []
    if not params.return_matching:
        return wasserstein_dist_(dgm_a, dgm_b, params, prices)
    else:
        dgm_a = [DiagramPoint(b, d, i, i) for i, (b, d) in enumerate(dgm_a)]
        dgm_b = [DiagramPoint(b, d, i, i) for i, (b, d) in enumerate(dgm_b)]
        dist = wasserstein_dist(dgm_a, dgm_b, params, prices)
        raise RuntimeError("add matching")
        return dist
