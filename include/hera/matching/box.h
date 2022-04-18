#ifndef MATCHING_DISTANCE_BOX_H
#define MATCHING_DISTANCE_BOX_H

#include <cassert>
#include <limits>

#include "common_util.h"

namespace md {

    template<class Real_>
    struct Box {
    public:
            using Real = Real_;
    private:
        Point<Real> ll;
        Point<Real> ur;

    public:
        Box(Point<Real> ll = Point<Real>(), Point<Real> ur = Point<Real>())
                :ll(ll), ur(ur)
        {
        }

        Box(Point<Real> center, Real width, Real height) :
            ll(Point<Real>(center.x - 0.5 * width, center.y - 0.5 * height)),
            ur(Point<Real>(center.x + 0.5 * width, center.y + 0.5 * height))
        {
        }


        inline double width() const { return ur.x - ll.x; }

        inline double height() const { return ur.y - ll.y; }

        inline Point<Real> lower_left() const { return ll; }
        inline Point<Real> upper_right() const { return ur; }
        inline Point<Real> center() const { return Point<Real>((ll.x + ur.x) / 2, (ll.y + ur.y) / 2); }

        inline bool operator==(const Box& p)
        {
            return this->ll == p.ll && this->ur == p.ur;
        }

        std::vector<Box> refine() const;

        std::vector<Point<Real>> corners() const;

        void translate(Real a);
    };

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const Box<Real>& box);

} // namespace md

#include "box.hpp"

#endif //MATCHING_DISTANCE_BOX_H
