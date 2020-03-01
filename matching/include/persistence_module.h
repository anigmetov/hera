#ifndef MATCHING_DISTANCE_PERSISTENCE_MODULE_H
#define MATCHING_DISTANCE_PERSISTENCE_MODULE_H

#include <cassert>
#include <vector>
#include <utility>
#include <string>

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


    class ModulePresentation {
    public:

        enum Format { rivet_firep };

        struct Relation {
            Point position_;
            IndexVec components_;

            Relation() {}
            Relation(const Point& _pos, const IndexVec& _components);

            Real get_x() const { return position_.x; }
            Real get_y() const { return position_.y; }
        };

        using RelVec = std::vector<Relation>;

        ModulePresentation() {}

        ModulePresentation(const PointVec& _generators, const RelVec& _relations);

        Diagram weighted_slice_diagram(const DualPoint& line) const;

        // translate all points by vector (a,a)
        void translate(Real a);

        // return minimal value of x- and y-coordinates
        Real minimal_coordinate() const { return std::min(min_x(), min_y()); }

        // return box that contains all positions of all simplices
        Box bounding_box() const;

        friend std::ostream& operator<<(std::ostream& os, const ModulePresentation& mp);

        Real max_x() const { return max_x_; }

        Real max_y() const { return max_y_; }

        Real min_x() const { return min_x_; }

        Real min_y() const { return min_y_; }

        PointVec positions() const;

    private:

        PointVec generators_;
        std::vector<Relation> relations_;
        PointVec positions_;


        Real max_x_ { std::numeric_limits<Real>::max() };
        Real max_y_ { std::numeric_limits<Real>::max() };
        Real min_x_ { -std::numeric_limits<Real>::max() };
        Real min_y_ { -std::numeric_limits<Real>::max() };
        Box bounding_box_;

        void init_boundaries();
        void project_generators(const DualPoint& slice, IndexVec& sorted_indices, RealVec& projections) const;
        void project_relations(const DualPoint& slice, IndexVec& sorted_indices, RealVec& projections) const;
    };
} // namespace md


#endif //MATCHING_DISTANCE_PERSISTENCE_MODULE_H
