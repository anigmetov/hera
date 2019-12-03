#include "catch/catch.hpp"

#include <sstream>
#include <iostream>

#include "common_util.h"
#include "box.h"
#include "bifiltration.h"

using namespace md;

//TEST_CASE("Small check", "[bifiltration][dim2]")
//{
//    Bifiltration bif("/home/narn/code/matching_distance/code/src/tests/test_bifiltration_full_triangle_rene.txt", BifiltrationFormat::rene);
//    auto simplices = bif.simplices();
//    bif.sanity_check();
//
//    REQUIRE( simplices.size() == 7 );
//
//    REQUIRE( simplices[0].dim() == 0 );
//    REQUIRE( simplices[1].dim() == 0 );
//    REQUIRE( simplices[2].dim() == 0 );
//    REQUIRE( simplices[3].dim() == 1 );
//    REQUIRE( simplices[4].dim() == 1 );
//    REQUIRE( simplices[5].dim() == 1 );
//    REQUIRE( simplices[6].dim() == 2);
//
//    REQUIRE( simplices[0].position() == Point(0, 0));
//    REQUIRE( simplices[1].position() == Point(0, 0));
//    REQUIRE( simplices[2].position() == Point(0, 0));
//    REQUIRE( simplices[3].position() == Point(3, 1));
//    REQUIRE( simplices[6].position() == Point(30, 40));
//
//    Line line_1(Line::pi / 2.0, 0.0);
//    auto dgm = bif.slice_diagram(line_1);
//}
