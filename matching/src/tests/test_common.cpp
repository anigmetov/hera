#include "catch/catch.hpp"

#include <sstream>
#include <iostream>
#include <string>

#include "common_util.h"
#include "simplex.h"
#include "matching_distance.h"

using namespace md;

TEST_CASE("Rational", "[common_utils][rational]")
{
    // gcd
    REQUIRE(gcd(10, 5) == 5);
    REQUIRE(gcd(5, 10) == 5);
    REQUIRE(gcd(5, 7) == 1);
    REQUIRE(gcd(7, 5) == 1);
    REQUIRE(gcd(13, 0) == 13);
    REQUIRE(gcd(0, 13) == 13);
    REQUIRE(gcd(16, 24) == 8);
    REQUIRE(gcd(24, 16) == 8);
    REQUIRE(gcd(16, 32) == 16);
    REQUIRE(gcd(32, 16) == 16);


    // reduce
    REQUIRE(reduce({2, 1}) == std::make_pair(2, 1));
    REQUIRE(reduce({1, 2}) == std::make_pair(1, 2));
    REQUIRE(reduce({2, 2}) == std::make_pair(1, 1));
    REQUIRE(reduce({0, 2}) == std::make_pair(0, 1));
    REQUIRE(reduce({0, 20}) == std::make_pair(0, 1));
    REQUIRE(reduce({35, 49}) == std::make_pair(5, 7));
    REQUIRE(reduce({35, 25}) == std::make_pair(7, 5));

    // midpoint
    REQUIRE(midpoint(Rational {0, 1}, Rational {1, 2}) == std::make_pair(1, 4));
    REQUIRE(midpoint(Rational {1, 4}, Rational {1, 2}) == std::make_pair(3, 8));
    REQUIRE(midpoint(Rational {1, 2}, Rational {1, 2}) == std::make_pair(1, 2));
    REQUIRE(midpoint(Rational {1, 2}, Rational {1, 1}) == std::make_pair(3, 4));
    REQUIRE(midpoint(Rational {3, 7}, Rational {5, 14}) == std::make_pair(11, 28));


    // arithmetic

    REQUIRE(Rational(1, 2) + Rational(3, 5) == Rational(11, 10));
    REQUIRE(Rational(2, 5) - Rational(3, 10) == Rational(1, 10));
    REQUIRE(Rational(2, 3) * Rational(4, 7) == Rational(8, 21));
    REQUIRE(Rational(2, 3) * Rational(3, 2) == Rational(1));
    REQUIRE(Rational(2, 3) / Rational(3, 2) == Rational(4, 9));
    REQUIRE(Rational(1, 2) * Rational(3, 5) == Rational(3, 10));

    // comparison
    REQUIRE(Rational(100000, 2000000) < Rational(100001, 2000000));
    REQUIRE(!(Rational(100001, 2000000) < Rational(100000, 2000000)));
    REQUIRE(!(Rational(100000, 2000000) < Rational(100000, 2000000)));
    REQUIRE(Rational(-100000, 2000000) < Rational(100001, 2000000));
    REQUIRE(Rational(-100001, 2000000) < Rational(100000, 2000000));
};

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
