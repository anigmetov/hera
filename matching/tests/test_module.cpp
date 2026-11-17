#include <catch2/catch.hpp>

#include <algorithm>
#include <limits>

#define MD_TEST_CODE
#include <hera/matching_distance.h>

using Real = double;
using Point = md::Point<Real>;
using Module = md::ModulePresentation<Real>;

TEST_CASE("Module slice barcode", "[module][barcode]")
{
    const Real inf = std::numeric_limits<Real>::infinity();
    const md::PointVec<Real> generators {Point(5, 5), Point(1, 1), Point(3, 3)};
    Module::RelVec relations {Module::Relation(Point(6, 6), md::IndexVec{0})};
    md::Diagram<Real> expected {{1, inf}, {3, inf}, {5, 6}};

    SECTION("Non-involutive generator order with fewer relations than generators") {}
    SECTION("Redundant and empty relations do not create module bars") {
        relations.emplace_back(Point(7, 7), md::IndexVec{0});
        relations.emplace_back(Point(0, 0), md::IndexVec{});
    }
    SECTION("Generators without relations remain essential") {
        relations.clear();
        expected = {{1, inf}, {3, inf}, {5, inf}};
    }
    SECTION("Zero-length bars are omitted") {
        relations = {Module::Relation(Point(5, 5), md::IndexVec{0})};
        expected = {{1, inf}, {3, inf}};
    }

    const Module module(generators, relations);
    const md::DualPoint<Real> diagonal(md::AxisType::y_type, md::AngleType::flat, 1, 0);
    auto diagram = module.weighted_slice_diagram(diagonal);
    std::sort(diagram.begin(), diagram.end());
    REQUIRE(diagram == expected);
}

TEST_CASE("Modules without generators have no slice bars", "[module][barcode]")
{
    const md::DualPoint<Real> diagonal(md::AxisType::y_type, md::AngleType::flat, 1, 0);
    const Module empty;
    REQUIRE(empty.weighted_slice_diagram(diagonal).empty());

    const Module relations_only({}, {Module::Relation(Point(1, 1), {})});
    REQUIRE(relations_only.weighted_slice_diagram(diagonal).empty());
}

TEST_CASE("Different essential ranks have infinite matching distance", "[module][matching_distance]")
{
    const Module rank_one({Point(0, 0)}, {});
    const Module rank_two({Point(0, 0), Point(1, 1)}, {});
    md::CalculationParams<Real> params;
    SECTION("Diagonal initialization") { params.initialization_depth = 0; }
    SECTION("Refined-grid initialization") { params.initialization_depth = 2; }
    params.max_depth = 0;
    params.hera_epsilon = 0;
    params.n_jobs = 1;
    REQUIRE(md::matching_distance(rank_one, rank_two, params) == std::numeric_limits<Real>::infinity());
    REQUIRE(params.actual_error == 0);
}
