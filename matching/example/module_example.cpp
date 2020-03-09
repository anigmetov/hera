#include <iostream>
#include "matching_distance.h"

using namespace md;

int main(int argc, char** argv)
{
    // create generators.
    // A generator is a point in plane,
    // generators are stored in a vector of points:
    PointVec<double> gens_a;

    // module A will have one generator that appears at point (0, 0)
    gens_a.emplace_back(0, 0);

    // relations are stored in a vector of relations
    using RelationVec = ModulePresentation<double>::RelVec;
    RelationVec rels_a;

    // A relation is a struct with position and column
    using Relation = ModulePresentation<double>::Relation;

    // at this point the relation rel_a_1 will appear:
    Point<double> rel_a_1_position { 1, 1 };

    // vector IndexVec contains non-zero indices of the corresponding relation
    // (we work over Z/2). Since we have one generator only, the relation
    // contains only one entry, 0
    IndexVec rel_a_1_components { 0 };

    // construct a relation from position and column:
    Relation rel_a_1 { rel_a_1_position, rel_a_1_components };

    // and add it to a vector of relations
    rels_a.push_back(rel_a_1);

    // after populating vectors of generators and relations
    // construct a module:
    ModulePresentation<double> module_a { gens_a, rels_a };


    // same for module_b. It will also have just one
    // generator and one relation, but at different positions.

    PointVec<double> gens_b;
    gens_b.emplace_back(1, 1);

    RelationVec rels_b;

    Point<double> rel_b_1_position { 2, 2 };
    IndexVec rel_b_1_components { 0 };

    rels_b.emplace_back(rel_b_1_position, rel_b_1_components);

    ModulePresentation<double> module_b { gens_b, rels_b };

    // create CalculationParams
    CalculationParams<double> params;
    // set relative error to 10 % :
    params.delta = 0.1;
    // go at most 8 levels deep in quadtree:
    params.max_depth = 8;

    double dist = matching_distance(module_a, module_b, params);
    std::cout << "dist = " << dist << std::endl;

    return 0;
}
