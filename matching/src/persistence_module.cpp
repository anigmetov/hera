#include <numeric>
#include <algorithm>
#include <unordered_set>

#include <phat/boundary_matrix.h>
#include <phat/compute_persistence_pairs.h>

#include "persistence_module.h"

namespace md {

    /**
     *
     * @param values vector of length n
     * @return [a_1,...,a_n] such that
     *           1) values[a_1] <= values[a_2] <= ... <= values[a_n]
     *           2) a_1,...,a_n is a permutation of 1,..,n
     */

    template<typename T>
    IndexVec get_sorted_indices(const std::vector<T>& values)
    {
        IndexVec result(values.size());
        std::iota(result.begin(), result.end(), 0);
        std::sort(result.begin(), result.end(),
                [&values](size_t a, size_t b) { return values[a] < values[b]; });
        return result;
    }

    // helper function to initialize const member positions_ in ModulePresentation
    PointVec
    concat_gen_and_rel_positions(const PointVec& generators, const ModulePresentation::RelVec& relations)
    {
        std::unordered_set<Point> ps(generators.begin(), generators.end());
        for(const auto& rel : relations) {
            ps.insert(rel.position_);
        }
        return PointVec(ps.begin(), ps.end());
    }


    void ModulePresentation::init_boundaries()
    {
        max_x_ = std::numeric_limits<Real>::max();
        max_y_ = std::numeric_limits<Real>::max();
        min_x_ = -std::numeric_limits<Real>::max();
        min_y_ = -std::numeric_limits<Real>::max();

        for(const auto& gen : positions_) {
            min_x_ = std::min(gen.x, min_x_);
            min_y_ = std::min(gen.y, min_y_);
            max_x_ = std::max(gen.x, max_x_);
            max_y_ = std::max(gen.y, max_y_);
        }

        bounding_box_ = Box(Point(min_x_, min_y_), Point(max_x_, max_y_));
    }


    ModulePresentation::ModulePresentation(const PointVec& _generators, const RelVec& _relations) :
        generators_(_generators),
        relations_(_relations)
    {
        init_boundaries();
    }

    void ModulePresentation::translate(md::Real a)
    {
        for(auto& g : generators_) {
            g.translate(a);
        }

        for(auto& r : relations_) {
            r.position_.translate(a);
        }

        positions_ = concat_gen_and_rel_positions(generators_, relations_);
        init_boundaries();
    }


    /**
     *
     * @param slice line on which generators are projected
     * @param sorted_indices [a_1,...,a_n] s.t. wpush(generator[a_1]) <= wpush(generator[a_2]) <= ..
     * @param projections sorted weighted pushes of generators
     */

    void
    ModulePresentation::project_generators(const DualPoint& slice, IndexVec& sorted_indices, RealVec& projections) const
    {
        size_t num_gens = generators_.size();

        RealVec gen_values;
        gen_values.reserve(num_gens);
        for(const auto& pos : generators_) {
            gen_values.push_back(slice.weighted_push(pos));
        }
        sorted_indices = get_sorted_indices(gen_values);
        projections.clear();
        projections.reserve(num_gens);
        for(auto i : sorted_indices) {
            projections.push_back(gen_values[i]);
        }
    }

    void ModulePresentation::project_relations(const DualPoint& slice, IndexVec& sorted_rel_indices,
            RealVec& projections) const
    {
        size_t num_rels = relations_.size();

        RealVec rel_values;
        rel_values.reserve(num_rels);
        for(const auto& rel : relations_) {
            rel_values.push_back(slice.weighted_push(rel.position_));
        }
        sorted_rel_indices = get_sorted_indices(rel_values);
        projections.clear();
        projections.reserve(num_rels);
        for(auto i : sorted_rel_indices) {
            projections.push_back(rel_values[i]);
        }
    }

    Diagram ModulePresentation::weighted_slice_diagram(const DualPoint& slice) const
    {
        IndexVec sorted_gen_indices, sorted_rel_indices;
        RealVec gen_projections, rel_projections;

        project_generators(slice, sorted_gen_indices, gen_projections);
        project_relations(slice, sorted_rel_indices, rel_projections);

        phat::boundary_matrix<> phat_matrix;

        phat_matrix.set_num_cols(relations_.size());

        for(Index i = 0; i < (Index) relations_.size(); i++) {
            IndexVec current_relation = relations_[sorted_rel_indices[i]].components_;
            for(auto& j : current_relation) {
                j = sorted_gen_indices[j];
            }
            std::sort(current_relation.begin(), current_relation.end());
            phat_matrix.set_dim(i, current_relation.size());
            phat_matrix.set_col(i, current_relation);
        }

        phat::persistence_pairs phat_persistence_pairs;
        phat::compute_persistence_pairs<phat::twist_reduction>(phat_persistence_pairs, phat_matrix);

        Diagram dgm;

        constexpr Real real_inf = std::numeric_limits<Real>::infinity();

        for(Index i = 0; i < (Index) phat_persistence_pairs.get_num_pairs(); i++) {
            std::pair<phat::index, phat::index> new_pair = phat_persistence_pairs.get_pair(i);
            bool is_finite_pair = new_pair.second != phat::k_infinity_index;
            Real birth = gen_projections.at(new_pair.first);
            Real death = is_finite_pair ? rel_projections.at(new_pair.second) : real_inf;
            if (birth != death) {
                dgm.emplace_back(birth, death);
            }
        }

        return dgm;
    }

    PointVec ModulePresentation::positions() const
    {
        return positions_;
    }

    Box ModulePresentation::bounding_box() const
    {
        return bounding_box_;
    }

}
