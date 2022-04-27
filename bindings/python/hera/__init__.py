#!/usr/bin/env python3

import typing

from ._hera import *
import numpy as np

def bottleneck_dist(dgm_a, dgm_b, delta: float=0.01, return_bottleneck_edge: bool=False):
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


def wasserstein_cost_g(pts_1, pts_2, prices):
    if type(pts_1) == np.ndarray:
        pts_1 = [ [x, y] for [x, y] in pts_1 ]

    if type(pts_2) == np.ndarray:
        pts_2 = [ [x, y] for [x, y] in pts_2 ]

    params = AuctionParams()
    dist = wasserstein_cost_geom(pts_1, pts_2, params, prices)
    return dist, params

def wasserstein_dist_1(dgm_a: typing.List[typing.Tuple[float, float]], dgm_b: typing.List[typing.Tuple[float, float]],
                     params: typing.Optional[AuctionParams]=None, prices=None):
    if params is None:
        params = AuctionParams()
    if prices is None:
        prices = []
    if not params.return_matching:
        return wasserstein_dist(dgm_a, dgm_b, params, prices)
    else:
        dgm_a = [ DiagramPoint(b, d, i, i) for i, (b, d) in enumerate(dgm_a) ]
        dgm_b = [ DiagramPoint(b, d, i, i) for i, (b, d) in enumerate(dgm_b) ]
        dist = wasserstein_dist(dgm_a, dgm_b, params, prices)
        return (dist, params.matching_a_to_b, params.matching_b_to_a)
