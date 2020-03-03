#include "catch/catch.hpp"

#include <sstream>
#include <iostream>
#include <string>

#include "common_util.h"
#include "simplex.h"
#include "matching_distance.h"

//using namespace md;
using Real = double;
using Point = md::Point<Real>;
using Bifiltration = md::Bifiltration<Real>;
using BifiltrationProxy = md::BifiltrationProxy<Real>;
using CalculationParams = md::CalculationParams<Real>;
using CellWithValue = md::CellWithValue<Real>;
using DualPoint = md::DualPoint<Real>;
using DualBox = md::DualBox<Real>;
using Simplex = md::Simplex<Real>;
using AbstractSimplex = md::AbstractSimplex;
using BoundStrategy = md::BoundStrategy;
using TraverseStrategy = md::TraverseStrategy;
using AxisType = md::AxisType;
using AngleType = md::AngleType;
using ValuePoint = md::ValuePoint;
using Column = md::Column;


TEST_CASE("AbstractSimplex", "[abstract_simplex]")
{
    AbstractSimplex as;
    REQUIRE(as.dim() == -1);

    as.push_back(1);
    REQUIRE(as.dim() == 0);
    REQUIRE(as.facets().size() == 0);

    as.push_back(0);
    REQUIRE(as.dim() == 1);
    REQUIRE(as.facets().size() == 2);
    REQUIRE(as.facets()[0].dim() == 0);
    REQUIRE(as.facets()[1].dim() == 0);

}

TEST_CASE("Vertical line", "[vertical_line]")
{
    // line x = 1
    DualPoint l_vertical(AxisType::x_type, AngleType::steep, 0, 1);

    REQUIRE(l_vertical.is_vertical());
    REQUIRE(l_vertical.is_steep());

    Point p_1(0.5, 0.5);
    Point p_2(1.5, 0.5);
    Point p_3(1.5, 1.5);
    Point p_4(0.5, 1.5);
    Point p_5(1, 10);

    REQUIRE(l_vertical.x_from_y(10) == 1);
    REQUIRE(l_vertical.x_from_y(-10) == 1);
    REQUIRE(l_vertical.x_from_y(0) == 1);

    REQUIRE(not l_vertical.contains(p_1));
    REQUIRE(not l_vertical.contains(p_2));
    REQUIRE(not l_vertical.contains(p_3));
    REQUIRE(not l_vertical.contains(p_4));
    REQUIRE(l_vertical.contains(p_5));

    REQUIRE(l_vertical.goes_below(p_1));
    REQUIRE(not l_vertical.goes_below(p_2));
    REQUIRE(not l_vertical.goes_below(p_3));
    REQUIRE(l_vertical.goes_below(p_4));

    REQUIRE(not l_vertical.goes_above(p_1));
    REQUIRE(l_vertical.goes_above(p_2));
    REQUIRE(l_vertical.goes_above(p_3));
    REQUIRE(not l_vertical.goes_above(p_4));

}

TEST_CASE("Horizontal line", "[horizontal_line]")
{
    // line y = 1
    DualPoint l_horizontal(AxisType::y_type, AngleType::flat, 0, 1);

    REQUIRE(l_horizontal.is_horizontal());
    REQUIRE(l_horizontal.is_flat());
    REQUIRE(l_horizontal.y_slope() == 0);
    REQUIRE(l_horizontal.y_intercept() == 1);

    Point p_1(0.5, 0.5);
    Point p_2(1.5, 0.5);
    Point p_3(1.5, 1.5);
    Point p_4(0.5, 1.5);
    Point p_5(2, 1);

    REQUIRE((not l_horizontal.contains(p_1) and
            not l_horizontal.contains(p_2) and
            not l_horizontal.contains(p_3) and
            not l_horizontal.contains(p_4) and
            l_horizontal.contains(p_5)));

    REQUIRE(not l_horizontal.goes_below(p_1));
    REQUIRE(not l_horizontal.goes_below(p_2));
    REQUIRE(l_horizontal.goes_below(p_3));
    REQUIRE(l_horizontal.goes_below(p_4));
    REQUIRE(l_horizontal.goes_below(p_5));

    REQUIRE(l_horizontal.goes_above(p_1));
    REQUIRE(l_horizontal.goes_above(p_2));
    REQUIRE(not l_horizontal.goes_above(p_3));
    REQUIRE(not l_horizontal.goes_above(p_4));
    REQUIRE(l_horizontal.goes_above(p_5));
}

TEST_CASE("Flat Line with positive slope", "[flat_line]")
{
    // line y = x / 2 + 1
    DualPoint l_flat(AxisType::y_type, AngleType::flat, 0.5, 1);

    REQUIRE(not l_flat.is_horizontal());
    REQUIRE(l_flat.is_flat());
    REQUIRE(l_flat.y_slope() == 0.5);
    REQUIRE(l_flat.y_intercept() == 1);

    REQUIRE(l_flat.y_from_x(0) == 1);
    REQUIRE(l_flat.y_from_x(1) == 1.5);
    REQUIRE(l_flat.y_from_x(2) == 2);

    Point p_1(3, 2);
    Point p_2(-2, 0.01);
    Point p_3(0, 1.25);
    Point p_4(6, 4.5);
    Point p_5(2, 2);

    REQUIRE((not l_flat.contains(p_1) and
            not l_flat.contains(p_2) and
            not l_flat.contains(p_3) and
            not l_flat.contains(p_4) and
            l_flat.contains(p_5)));

    REQUIRE(not l_flat.goes_below(p_1));
    REQUIRE(l_flat.goes_below(p_2));
    REQUIRE(l_flat.goes_below(p_3));
    REQUIRE(l_flat.goes_below(p_4));
    REQUIRE(l_flat.goes_below(p_5));

    REQUIRE(l_flat.goes_above(p_1));
    REQUIRE(not l_flat.goes_above(p_2));
    REQUIRE(not l_flat.goes_above(p_3));
    REQUIRE(not l_flat.goes_above(p_4));
    REQUIRE(l_flat.goes_above(p_5));

}
