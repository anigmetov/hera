#!/usr/bin/env python3

import hera

def try_bt():
    dgm1 = [ [1, 2], [3, 3.2]]
    dgm2 = [ [ 1, 2.3], [3, 3.4], [5, 5.01]]
    dist, edge = hera.bottleneck_dist(dgm1, dgm2, return_bottleneck_edge=True)
    print(dist, edge)
    dist = hera.bottleneck_dist(dgm1, dgm2, delta=0)
    print(dist)


def try_ws():
    dgm1 = [ [1, 2], [3, 3.2]]
    dgm2 = [ [ 1, 2.3], [3, 3.4], [5, 5.01]]
    dist = hera.wasserstein_dist(dgm1, dgm2)
    print(dist)
    params = hera.WassersteinParams()
    params.delta = 0.1
    params.wasserstein_power = 2.0
    dist = hera.wasserstein_dist(dgm1, dgm2, params)
    print(dist)


def try_ws_geom():
    dgm1 = [ [1, 2], [3, 3.2], [4.8, 4.9]]
    dgm2 = [ [ 1, 2.3], [3, 3.4], [5, 5.01]]
    dist = hera.wasserstein_cost_geom(dgm1, dgm2)
    print(dist)
    params = hera.WassersteinParams()
    params.delta = 0.1
    params.wasserstein_power = 2.0
    dist = hera.wasserstein_cost_geom(dgm1, dgm2, params)
    print(dist)
    res = hera.wasserstein_cost_geom_detailed_(dgm1, dgm2, params, [])
    print(res)


if __name__ == "__main__":
    try_bt()
    try_ws()
    try_ws_geom()
