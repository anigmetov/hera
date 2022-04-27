#!/usr/bin/env python3

import hera

def try_bt():
    dgm1 = [ [1, 2], [3, 3.2]]
    dgm2 = [ [ 1, 2.3], [3, 3.4], [5, 5.01]]
    dist, edge = hera.bottleneck_dist(dgm1, dgm2, return_bottleneck_edge=True)
    print(dist, edge)

def try_ws():
    dgm1 = [ [1, 2], [3, 3.2]]
    dgm2 = [ [ 1, 2.3], [3, 3.4], [5, 5.01]]

    pass

if __name__ == "__main__":
    try_bt()