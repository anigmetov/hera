#include <iostream>
#include <hera/matching_distance.h>

using namespace md;

int main(int /*argc*/, char** /*argv*/)
{
    using RelationVec = ModulePresentation<double>::RelVec;
    using Relation = ModulePresentation<double>::Relation;

    // Module A (three rectangles)
    PointVec<double> gens_a;
    RelationVec rels_a;


    // First rectangle
    gens_a.emplace_back(-3, -1);

    Point<double> rel_a_1_position { -3, 3 };
    IndexVec rel_a_1_components { 0 };
    Relation rel_a_1 { rel_a_1_position, rel_a_1_components };
    rels_a.push_back(rel_a_1);

    Point<double> rel_a_2_position { 1, -1 };
    IndexVec rel_a_2_components { 0 };
    Relation rel_a_2 { rel_a_2_position, rel_a_2_components };
    rels_a.push_back(rel_a_2);

    // Second rectangle
    gens_a.emplace_back(-1, -1);

    Point<double> rel_a_3_position { -1, 1 };
    IndexVec rel_a_3_components { 1 };
    Relation rel_a_3 { rel_a_3_position, rel_a_3_components };
    rels_a.push_back(rel_a_3);

    Point<double> rel_a_4_position { 1, -1 };
    IndexVec rel_a_4_components { 1 };
    Relation rel_a_4 { rel_a_4_position, rel_a_4_components };
    rels_a.push_back(rel_a_4);

    // Third rectangle
    gens_a.emplace_back(-1, -3);

    Point<double> rel_a_5_position { -1, 1 };
    IndexVec rel_a_5_components { 2 };
    Relation rel_a_5 { rel_a_5_position, rel_a_5_components };
    rels_a.push_back(rel_a_5);

    Point<double> rel_a_6_position { 3, -3 };
    IndexVec rel_a_6_components { 2 };
    Relation rel_a_6 { rel_a_6_position, rel_a_6_components };
    rels_a.push_back(rel_a_6);


    ModulePresentation<double> module_a { gens_a, rels_a };


    // Module B (one rectangle)
    PointVec<double> gens_b;
    RelationVec rels_b;

    // Rectangle
    gens_b.emplace_back(-2, -2);

    Point<double> rel_b_1_position { -2, 2 };
    IndexVec rel_b_1_components { 0 };
    Relation rel_b_1 { rel_b_1_position, rel_b_1_components };
    rels_b.push_back(rel_b_1);

    Point<double> rel_b_2_position { 2, -2 };
    IndexVec rel_b_2_components { 0 };
    Relation rel_b_2 { rel_b_2_position, rel_b_2_components };
    rels_b.push_back(rel_b_2);


    ModulePresentation<double> module_b { gens_b, rels_b };


    // Computations
    CalculationParams<double> params;
    params.delta = 0.1;
    params.max_depth = 5;
    params.initialization_depth = 0;

    double dist = matching_distance(module_a, module_b, params);
    std::cout << "dist = " << dist << std::endl;

    return 0;
}
