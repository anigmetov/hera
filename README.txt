This repository contains software to compute bottleneck and Wasserstein
distances between persistence diagrams.
See Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov,
"Geometry Helps to Compare Persistence Diagrams.", ALENEX 2016.
http://dx.doi.org/10.1137/1.9781611974317.9

Unfortunately, the license for the two parts is different. The bottleneck 
part (geom_bottleneck/) is licensed under LGPL (because it uses ANN library), 
the Wasserstein part (in geom_matching/) is licensed under less restrictive license. 
If you are going to use this software for research purposes,
you probably do not need to worry about that.

See README files in subdirectories for usage and building.