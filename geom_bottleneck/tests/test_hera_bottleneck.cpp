#include "catch/catch.hpp"

#include <sstream>
#include <iostream>

#include "bottleneck.h"

using PairVector = std::vector<std::pair<double, double>>;

std::vector<std::string> split_on_delim(const std::string& s, char delim)
{
    std::stringstream ss(s);
    std::string token;
    std::vector<std::string> tokens;
    while(std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}


// single row in a file with test cases
struct TestFromFileCase {

    std::string file_1;
    std::string file_2;
    double q;
    double internal_p;
    double answer;

    TestFromFileCase(std::string s)
    {
        auto tokens = split_on_delim(s, ' ');
        assert(tokens.size() == 5);

        file_1 = tokens.at(0);
        file_2 = tokens.at(1);
        q = std::stod(tokens.at(2));
        internal_p = std::stod(tokens.at(3));
        answer = std::stod(tokens.at(4));

        if ( q < 1.0 or std::isinf(q) or
            (internal_p != hera::get_infinity<double>() and internal_p < 1.0)) {
            throw std::runtime_error("Bad line in test_list.txt");
        }
    }
};

std::ostream& operator<<(std::ostream& out, const TestFromFileCase& s)
{
    out << "[" << s.file_1 << ", " << s.file_2 << ", q = " << s.q << ", norm = ";
    if (s.internal_p != hera::get_infinity()) {
        out << s.internal_p;
    } else {
        out << "infinity";
    }
    out << ", answer = " << s.answer << "]";
    return out;
}


TEST_CASE("simple cases", "bottleneckDistApprox")
{
    PairVector diagram_A, diagram_B;
    double delta = 0.01;
    //double internal_p = hera::get_infinity<double>();

    SECTION("trivial: two empty diagrams") {
        REQUIRE(  0.0 == hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta));
    }

    SECTION("trivial: one empty diagram, one single-point diagram") {

        diagram_A.emplace_back(1.0, 2.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        REQUIRE(  fabs(d1 - 0.5) <= 0.00000000001 );

        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        REQUIRE(  fabs(d2 - 0.5) <= 0.00000000001 );
    }

    SECTION("trivial: two single-point diagrams-1") {

        diagram_A.emplace_back(10.0, 20.0);  // (5, 5)
        diagram_B.emplace_back(13.0, 19.0);  // (3, 3)

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double correct_answer = 3.0;
        REQUIRE(  fabs(d1 - correct_answer) <= delta * correct_answer);
        REQUIRE(  fabs(d2 - correct_answer) <= delta * correct_answer);
    }

    SECTION("trivial: two single-point diagrams-2") {

        diagram_A.emplace_back(10.0, 20.0);  // (5, 5)
        diagram_B.emplace_back(130.0, 138.0);  // (4, 4)

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double correct_answer = 5.0;
        REQUIRE(  fabs(d1 - correct_answer) <= delta * correct_answer );
        REQUIRE(  fabs(d2 - correct_answer) <= delta * correct_answer );

    }

}


TEST_CASE("infinity points", "bottleneckDistApprox")
{
    PairVector diagram_A, diagram_B;
    double delta = 0.01;

    // do not use Hera's infinity! it is -1
    double inf = std::numeric_limits<double>::infinity();

    SECTION("two points at infinity, no finite points") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        double d = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double corr_answer = 1.0;
        REQUIRE(  fabs(d - corr_answer) <= delta * corr_answer);
    }

    SECTION("two points at infinity") {

        // edge cost 3.0
        diagram_A.emplace_back(10.0, 20.0);  // (5, 5)
        diagram_B.emplace_back(13.0, 19.0);  // (3, 3)

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        double d = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double corr_answer = 3.0;
        REQUIRE(  fabs(d - corr_answer) <= delta * corr_answer);
    }

    SECTION("three points at infinity, no finite points") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);
        diagram_B.emplace_back(2.0, inf);

        double d = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double corr_answer = inf;
        REQUIRE( d  == corr_answer );
    }

    SECTION("three points at infinity") {

        // edge cost 3.0
        diagram_A.emplace_back(10.0, 20.0);  // (5, 5)
        diagram_B.emplace_back(13.0, 19.0);  // (3, 3)

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        double d = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double corr_answer = inf;
        REQUIRE( d  == corr_answer );
    }


    SECTION("all four corners at infinity, no finite points, finite answer") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        // edge cost 1.0
        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        // edge cost 1.0
        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        // edge cost 1.0
        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        double d = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double corr_answer = 1.0;

        REQUIRE( d  == corr_answer );
    }

    SECTION("all four corners at infinity, no finite points, infinite answer-1") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        // edge cost 1.0
        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        // edge cost 1.0
        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        // edge cost 1.0
        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, no finite points, infinite answer-2") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        // edge cost 1.0
        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        // edge cost 1.0
        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        // edge cost 1.0
        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, no finite points, infinite answer-3") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        // edge cost 1.0
        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        // edge cost 1.0
        diagram_A.emplace_back(inf, 1.0);
        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        // edge cost 1.0
        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, no finite points, infinite answer-4") {

        // edge cost 1.0
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        // edge cost 1.0
        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        // edge cost 1.0
        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        // edge cost 1.0
        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);
        diagram_B.emplace_back(-inf, 2.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, with finite points, infinite answer-1") {

        diagram_A.emplace_back(1.0, inf);
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        // finite edge
        diagram_A.emplace_back(10.0, 20.0);
        diagram_B.emplace_back(13.0, 19.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, with finite points, infinite answer-2") {

        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        // finite edge
        diagram_A.emplace_back(10.0, 20.0);
        diagram_B.emplace_back(13.0, 19.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, with finite points, infinite answer-3") {

        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        diagram_A.emplace_back(inf, 1.0);
        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);

        // finite edge
        diagram_A.emplace_back(10.0, 20.0);
        diagram_B.emplace_back(13.0, 19.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }

    SECTION("all four corners at infinity, no finite points, infinite answer-4") {

        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        diagram_A.emplace_back(1.0, -inf);
        diagram_B.emplace_back(2.0, -inf);

        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        diagram_A.emplace_back(-inf, 1.0);
        diagram_B.emplace_back(-inf, 2.0);
        diagram_B.emplace_back(-inf, 2.0);

        // finite edge
        diagram_A.emplace_back(10.0, 20.0);
        diagram_B.emplace_back(13.0, 19.0);

        double d1 = hera::bottleneckDistApprox<>(diagram_A, diagram_B, delta);
        double d2 = hera::bottleneckDistApprox<>(diagram_B, diagram_A, delta);
        double corr_answer = inf;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );
    }


    SECTION("simple small example with finite answer") {
        diagram_A.emplace_back(1.0, inf);
        diagram_B.emplace_back(2.0, inf);

        diagram_A.emplace_back(1.9, inf);
        diagram_B.emplace_back(1.1, inf);

        // 1.1 - 1.0 +  2.0 - 1.9 = 0.2

        diagram_A.emplace_back(inf, 1.0);
        diagram_B.emplace_back(inf, 2.0);

        diagram_A.emplace_back(inf, 1.9);
        diagram_B.emplace_back(inf, 1.1);


        // finite edge
        diagram_A.emplace_back(10.0, 20.0);
        diagram_B.emplace_back(13.0, 19.0);

        double d1 = hera::bottleneckDistExact<>(diagram_A, diagram_B);
        double d2 = hera::bottleneckDistExact<>(diagram_B, diagram_A);
        double corr_answer = 3.0;

        REQUIRE( d1 == corr_answer );
        REQUIRE( d2 == corr_answer );


     }

}

