#include "catch/catch.hpp"

#include <sstream>
#include <iostream>
#include <string>

#define MD_TEST_CODE

#include "common_util.h"
#include "persistence_module.h"
#include "matching_distance.h"

using Real = double;
using Point = md::Point<Real>;
using Bifiltration = md::Bifiltration<Real>;
using BifiltrationProxy = md::BifiltrationProxy<Real>;
using CalculationParams = md::CalculationParams<Real>;
using CellWithValue = md::CellWithValue<Real>;
using DualPoint = md::DualPoint<Real>;
using DualBox = md::DualBox<Real>;
using BoundStrategy = md::BoundStrategy;
using TraverseStrategy = md::TraverseStrategy;
using AxisType = md::AxisType;
using AngleType = md::AngleType;
using ValuePoint = md::ValuePoint;
using Column = md::Column;
using PointVec = md::PointVec<Real>;
using Module = md::ModulePresentation<Real>;
using Relation = Module::Relation;
using RelationVec = Module::RelVec;
using IndexVec = md::IndexVec;

using md::k_corner_vps;

TEST_CASE("Module projection", "[module][projection]")
{
    PointVec gens;
    gens.emplace_back(1, 1); // A
    gens.emplace_back(2, 3); // B
    gens.emplace_back(3, 2); // C

    RelationVec rels;

    Point rel_x_position { 3.98, 2.47 };
    IndexVec rel_x_components { 0, 2 }; // X: A + C = 0

    Point rel_y_position { 2.5, 4 };
    IndexVec rel_y_components { 0, 1 }; // Y: A + B = 0

    Point rel_z_position { 5, 5 };
    IndexVec rel_z_components { 1 };   // Z: B = 0


    rels.emplace_back(rel_x_position, rel_x_components);
    rels.emplace_back(rel_y_position, rel_y_components);
    rels.emplace_back(rel_z_position, rel_z_components);

    Module module { gens, rels };

    {
        DualPoint slice_1(AxisType::x_type, AngleType::flat, 6.0 / 7.0, 3.0);
        std::vector<Real> gen_ps_1, rel_ps_1;
        phat::boundary_matrix<> matr_1;

        module.get_slice_projection_matrix(slice_1, matr_1, gen_ps_1, rel_ps_1);

        phat::column c_1_0, c_1_1, c_1_2;

        matr_1.get_col(0, c_1_0);
        matr_1.get_col(1, c_1_1);
        matr_1.get_col(2, c_1_2);

        phat::column c_1_0_correct { 0, 1};
        phat::column c_1_1_correct { 0, 2};
        phat::column c_1_2_correct { 2 };

        REQUIRE(c_1_0 == c_1_0_correct);
        REQUIRE(c_1_1 == c_1_1_correct);
        REQUIRE(c_1_2 == c_1_2_correct);
    }

    {

        DualPoint slice_2(AxisType::y_type, AngleType::flat, 2.0 / 9.0, 5.0);
        std::vector<Real> gen_ps_2, rel_ps_2;
        phat::boundary_matrix<> matr_2;

        module.get_slice_projection_matrix(slice_2, matr_2, gen_ps_2, rel_ps_2);

        phat::column c_2_0, c_2_1, c_2_2;

        matr_2.get_col(0, c_2_0);
        matr_2.get_col(1, c_2_1);
        matr_2.get_col(2, c_2_2);

        phat::column c_2_0_correct { 0, 1};
        phat::column c_2_1_correct { 0, 2};
        phat::column c_2_2_correct { 1 };

        //std::cerr << "gen_ps_2: " << md::container_to_string(gen_ps_2) << std::endl;
        //std::cerr << "rel_ps_2: " << md::container_to_string(rel_ps_2) << std::endl;

        REQUIRE(c_2_0 == c_2_0_correct);
        REQUIRE(c_2_1 == c_2_1_correct);
        REQUIRE(c_2_2 == c_2_2_correct);
    }


}
