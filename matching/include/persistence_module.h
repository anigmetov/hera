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

    class Relation {
    public:
        Relation() {}

        Relation(const Point& _pos)
                :position_(_pos) { }

        Real get_x() const { return position_.x; }

        Real get_y() const { return position_.y; }

    private:
        Point position_;
    };

    class PersistenceModule {
    public:
        using Box = md::Box;

        PersistenceModule() { }

        PersistenceModule(const std::string& fname); // read from file

        Diagram slice_diagram(const DualPoint& line);

        Box bounding_box() const;

    private:
        std::vector<Point> generators_;
        std::vector<Relation> relations_;

    };

//    template<typename R = double>
//    R matching_distance(const PersistenceModule&, const PersistenceModule&, R)
//    {
//        return 0.0;
//    };

} // namespace md



#endif //MATCHING_DISTANCE_PERSISTENCE_MODULE_H
