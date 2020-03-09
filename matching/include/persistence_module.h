#ifndef MATCHING_DISTANCE_PERSISTENCE_MODULE_H
#define MATCHING_DISTANCE_PERSISTENCE_MODULE_H

#include <cassert>
#include <vector>
#include <utility>
#include <string>
#include <numeric>
#include <algorithm>
#include <unordered_set>

#include "phat/boundary_matrix.h"
#include "phat/compute_persistence_pairs.h"

#include "common_util.h"
#include "dual_point.h"
#include "box.h"

namespace md {

    /* ModulePresentation contains only information needed for matching
     * distance computation over Z/2.
     * Generators are represented as points (critical values),
     * id i of generator g_i = its index in * vector generators_.
     *
     * Relations are represented by struct Relation, which has two members:
     *   position_ is a point at which relation appears,
     *   components_ contains indices of generators that sum to zero:
     *     if components_ = [i, ..., j], then the relation is g_i +...+ g_j = 0.
     *
     *  ModulePresentation has member positions_ that contains all
     *  distinct positions of generators and relations;
     *  this member simplifies computing local linear bound.
     */


    template<class Real>
    class ModulePresentation {
    public:

        using RealVec = std::vector<Real>;

        enum Format { rivet_firep };

        struct Relation {
            Point<Real> position_;
            IndexVec components_;

            Relation() {}
            Relation(const Point<Real>& _pos, const IndexVec& _components) : position_(_pos), components_(_components) {}

            Real get_x() const { return position_.x; }
            Real get_y() const { return position_.y; }
        };

        using RelVec = std::vector<Relation>;

        ModulePresentation() {}

        ModulePresentation(const PointVec<Real>& _generators, const RelVec& _relations);

        Diagram<Real> weighted_slice_diagram(const DualPoint<Real>& line) const;

        // translate all points by vector (a,a)
        void translate(Real a);

        // return minimal value of x- and y-coordinates
        Real minimal_coordinate() const { return std::min(min_x(), min_y()); }

        // return box that contains all positions of all simplices
        Box<Real> bounding_box() const;

        Real max_x() const { return max_x_; }

        Real max_y() const { return max_y_; }

        Real min_x() const { return min_x_; }

        Real min_y() const { return min_y_; }

        PointVec<Real> positions() const;

#ifndef MD_TEST_CODE
    private:
#endif

        PointVec<Real> generators_;
        std::vector<Relation> relations_;
        PointVec<Real> positions_;


        Real max_x_ { std::numeric_limits<Real>::max() };
        Real max_y_ { std::numeric_limits<Real>::max() };
        Real min_x_ { -std::numeric_limits<Real>::max() };
        Real min_y_ { -std::numeric_limits<Real>::max() };
        Box<Real> bounding_box_;

        void init_boundaries();

        void project_generators(const DualPoint<Real>& slice, IndexVec& sorted_indices, RealVec& projections) const;
        void project_relations(const DualPoint<Real>& slice, IndexVec& sorted_indices, RealVec& projections) const;

        void get_slice_projection_matrix(const DualPoint<Real>& slice, phat::boundary_matrix<>& phat_matrix,
                RealVec& gen_projections, RealVec& rel_projections) const;

    };
} // namespace md

#include "persistence_module.hpp"

#endif //MATCHING_DISTANCE_PERSISTENCE_MODULE_H
