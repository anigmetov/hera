#!/usr/bin/env python3

import typing

from ._hera import *

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


def wasserstein_dist_1(dgm_a: typing.List[typing.Tuple[float, float]], dgm_b: typing.List[typing.Tuple[float, float]],
                     params: typing.Optional[AuctionParams]=None, prices=None):
    if params is None:
        params = AuctionParams()
    if prices is None:
        prices = []

    dist = wasserstein_dist(dgm_a, dgm_b, params, prices)
    if params.return_matching:
        return (dist, params.matching_a_to_b, params.matching_b_to_a)
    else:
        return dist
