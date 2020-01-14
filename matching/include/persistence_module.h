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


    class ModulePresentation {
    public:

        struct Relation {

            Point position_;
            IntPoint bigrade_;

            Relation() {}
            Relation(const Point& _pos) : position_(_pos) { }

            Real get_x() const { return position_.x; }
            Real get_y() const { return position_.y; }

            IndexVec components_;

        };

        ModulePresentation() { }

        ModulePresentation(const std::string& fname); // read from file

        Diagram weighted_slice_diagram(const DualPoint& line) const;

        // translate all points by vector (a,a)
        void translate(Real a);

        // return minimal value of x- and y-coordinates
        Real minimal_coordinate() const { return std::min(min_x(), min_x()); }

        // return box that contains all positions of all simplices
        Box bounding_box() const;

        friend std::ostream& operator<<(std::ostream& os, const ModulePresentation& mp);

        Real max_x() const { return max_x_; }

        Real max_y() const { return max_y_; }

        Real min_x() const { return min_x_; }

        Real min_y() const { return min_y_; }

        PointVec positions() const;

    private:

        Real max_x_ { std::numeric_limits<Real>::max() };
        Real max_y_ { std::numeric_limits<Real>::max() };
        Real min_x_ { -std::numeric_limits<Real>::max() };
        Real min_y_ { -std::numeric_limits<Real>::max() };
        Box bounding_box_;

        PointVec generators_;
        std::vector<Relation> relations_;

        PointVec positions_;

        void init_boundaries();
        void project_generators(const DualPoint& slice, IndexVec& sorted_indices, RealVec& projections) const;
        void project_relations(const DualPoint& slice, IndexVec& sorted_indices, RealVec& projections) const;
    };

//    template<typename R = double>
//    R matching_distance(const PersistenceModule&, const PersistenceModule&, R)
//    {
//        return 0.0;
//    };

} // namespace md



#endif //MATCHING_DISTANCE_PERSISTENCE_MODULE_H
