#include <map>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>

#include "opts/opts.h"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

#include "common_util.h"
#include "bifiltration.h"

using Index = md::Index;
using Point = md::Point;
using Column = md::Column;

int g_max_coord = 100;

using ASimplex = md::AbstractSimplex;

using ASimplexToBirthMap = std::map<ASimplex, Point>;

namespace spd = spdlog;

// random generator is global
std::random_device rd;
std::mt19937_64 gen(rd());

//std::mt19937_64 gen(42);

Point get_random_position(int max_coord)
{
    assert(max_coord > 0);
    std::uniform_int_distribution<int> distr(0, max_coord);
    return Point(distr(gen), distr(gen));

}

Point get_random_position_less_than(Point ub)
{
    std::uniform_int_distribution<int> distr_x(0, ub.x);
    std::uniform_int_distribution<int> distr_y(0, ub.y);
    return Point(distr_x(gen), distr_y(gen));
}

Point get_random_position_greater_than(Point lb)
{
    std::uniform_int_distribution<int> distr_x(lb.x, g_max_coord);
    std::uniform_int_distribution<int> distr_y(lb.y, g_max_coord);
    return Point(distr_x(gen), distr_y(gen));
}

// non-proper faces (empty and simplex itself) are also faces
bool is_face(const ASimplex& face_candidate, const ASimplex& coface_candidate)
{
    return std::includes(coface_candidate.begin(), coface_candidate.end(), face_candidate.begin(),
            face_candidate.end());
}

bool is_top_simplex(const ASimplex& candidate_simplex, const std::vector<ASimplex>& current_top_simplices)
{
    return std::none_of(current_top_simplices.begin(), current_top_simplices.end(),
            [&candidate_simplex](const ASimplex& ts) { return is_face(candidate_simplex, ts); });
}

void add_if_top(const ASimplex& candidate_simplex, std::vector<ASimplex>& current_top_simplices)
{
    // check that candidate simplex is not face of someone in current_top_simplices
    if (!is_top_simplex(candidate_simplex, current_top_simplices))
        return;

    spd::debug("candidate_simplex is top, will be added to top_simplices");
    // remove s from currrent_top_simplices, if s is face of candidate_simplex
    current_top_simplices.erase(std::remove_if(current_top_simplices.begin(), current_top_simplices.end(),
            [candidate_simplex](const ASimplex& s) { return is_face(s, candidate_simplex); }),
            current_top_simplices.end());

    current_top_simplices.push_back(candidate_simplex);
}

ASimplex get_random_simplex(int n_vertices, int dim)
{
    std::vector<int> all_vertices(n_vertices, 0);
    // fill in with 0..n-1
    std::iota(all_vertices.begin(), all_vertices.end(), 0);
    std::shuffle(all_vertices.begin(), all_vertices.end(), gen);
    return ASimplex(all_vertices.begin(), all_vertices.begin() + dim + 1, true);
}

void generate_positions(const ASimplex& s, ASimplexToBirthMap& simplex_to_birth, Point upper_bound)
{
    auto pos = get_random_position_less_than(upper_bound);
    auto curr_pos_iter = simplex_to_birth.find(s);
    if (curr_pos_iter != simplex_to_birth.end())
        pos = md::greatest_lower_bound(pos, curr_pos_iter->second);
    simplex_to_birth[s] = pos;
    for(const ASimplex& facet : s.facets()) {
        generate_positions(facet, simplex_to_birth, pos);
    }
}

md::Bifiltration get_random_bifiltration(int n_vertices, int max_dim, int n_top_simplices)
{
    ASimplexToBirthMap simplex_to_birth;

    // generate vertices
    for(int i = 0; i < n_vertices; ++i) {
        Point vertex_birth = get_random_position(g_max_coord / 10);
        ASimplex vertex;
        vertex.push_back(i);
        simplex_to_birth[vertex] = vertex_birth;
    }

    std::vector<ASimplex> top_simplices;
    // generate top simplices
    while((int)top_simplices.size() < n_top_simplices) {
        std::uniform_int_distribution<int> dimension_distr(1, max_dim);
        int dim = dimension_distr(gen);
        auto candidate_simplex = get_random_simplex(n_vertices, dim);
        spd::debug("candidate_simplex = {}", candidate_simplex);
        add_if_top(candidate_simplex, top_simplices);
    }

    Point upper_bound{static_cast<md::Real>(g_max_coord), static_cast<md::Real>(g_max_coord)};
    for(const auto& top_simplex : top_simplices) {
        generate_positions(top_simplex, simplex_to_birth, upper_bound);
    }

    std::vector<std::pair<ASimplex, Point>> simplex_birth_pairs{simplex_to_birth.begin(), simplex_to_birth.end()};
    std::vector<md::Column> boundaries{simplex_to_birth.size(), md::Column()};

// assign ids and save boundaries
    int id = 0;

    for(int dim = 0; dim <= max_dim; ++dim) {
        for(int i = 0; i < (int) simplex_birth_pairs.size(); ++i) {
            ASimplex& simplex = simplex_birth_pairs[i].first;
            if (simplex.dim() == dim) {
                simplex.id = id++;
                md::Column bdry;
                for(auto& facet : simplex.facets()) {
                    auto facet_iter = std::find_if(simplex_birth_pairs.begin(), simplex_birth_pairs.end(),
                            [facet](const std::pair<ASimplex, Point>& sbp) { return facet == sbp.first; });
                    assert(facet_iter != simplex_birth_pairs.end());
                    assert(facet_iter->first.id >= 0);
                    bdry.push_back(facet_iter->first.id);
                }
                std::sort(bdry.begin(), bdry.end());
                boundaries[i] = bdry;
            }
        }
    }

// create vector of Simplex-es
    std::vector<md::Simplex> simplices;
    for(int i = 0; i < (int) simplex_birth_pairs.size(); ++i) {
        int id = simplex_birth_pairs[i].first.id;
        int dim = simplex_birth_pairs[i].first.dim();
        Point birth = simplex_birth_pairs[i].second;
        Column bdry = boundaries[i];
        simplices.emplace_back(id, birth, dim, bdry);
    }

// sort by id
    std::sort(simplices.begin(), simplices.end(),
            [](const md::Simplex& s1, const md::Simplex& s2) { return s1.id() < s2.id(); });
    for(int i = 0; i < (int)simplices.size(); ++i) {
        assert(simplices[i].id() == i);
        assert(i == 0 || simplices[i].dim() >= simplices[i - 1].dim());
    }

    return md::Bifiltration(simplices.begin(), simplices.end());
}

int main(int argc, char** argv)
{
    spd::set_level(spd::level::info);
    int n_vertices;
    int max_dim;
    int n_top_simplices;
    std::string fname;

    using opts::Option;
    using opts::PosOption;
    opts::Options ops;

    bool help = false;

    ops >> Option('v', "n-vertices", n_vertices, "number of vertices")
        >> Option('d', "max-dim", max_dim, "maximal dim")
        >> Option('m', "max-coord", g_max_coord, "maximal coordinate")
        >> Option('t', "n-top-simplices", n_top_simplices, "number of top simplices")
        >> Option('h', "help", help, "show help message");

    if (!ops.parse(argc, argv) || help || !(ops >> PosOption(fname))) {
        std::cerr << "Usage: " << argv[0] << "\n" << ops << std::endl;
        return 1;
    }


    auto bif1 = get_random_bifiltration(n_vertices, max_dim, n_top_simplices);
    std::cout << "Generated bifiltration." << std::endl;
    bif1.save(fname, md::BifiltrationFormat::phat_like);
    std::cout << "Saved to file " << fname << std::endl;
    return 0;
}

