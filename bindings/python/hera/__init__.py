from ._hera import *
import numpy as np

def bottleneck_dist(dgm_a, dgm_b, delta: float=0.01, return_bottleneck_edge: bool=False):
    if delta < 0.0:
        raise RuntimeError(f"Relative error delta must be non-negative, got {delta}")
    if delta == 0.0:
        if return_bottleneck_edge:
            return bottleneck_distance_exact_with_edge(dgm_a, dgm_b)
        else:
            return bottleneck_distance_exact(dgm_a, dgm_b)
    else:
        if return_bottleneck_edge:
            return bottleneck_distance_approx_with_edge(dgm_a, dgm_b, delta)
        else:
            return bottleneck_distance_approx(dgm_a, dgm_b, delta)


def wasserstein_cost(pts_1, pts_2, prices):
    if type(pts_1) == np.ndarray:
        pts_1 = [ [x, y] for [x, y] in pts_1 ]

    if type(pts_2) == np.ndarray:
        pts_2 = [ [x, y] for [x, y] in pts_2 ]

    params = AuctionParams()
    dist = wasserstein_cost_geom(pts_1, pts_2, params, prices)
    return dist, params

