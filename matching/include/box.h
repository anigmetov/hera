#ifndef MATCHING_DISTANCE_BOX_H
#define MATCHING_DISTANCE_BOX_H

#include <cassert>
#include <limits>

#include "common_util.h"

namespace md {

    struct Box {
    private:
        Point ll;
        Point ur;

    public:
        Box(Point ll = Point(), Point ur = Point())
                :ll(ll), ur(ur)
        {
        }

        Box(Point center, Real width, Real height) :
            ll(Point(center.x - 0.5 * width, center.y - 0.5 * height)),
            ur(Point(center.x + 0.5 * width, center.y + 0.5 * height))
        {
        }


        inline double width() const { return ur.x - ll.x; }

        inline double height() const { return ur.y - ll.y; }

        inline Point lower_left() const { return ll; }
        inline Point upper_right() const { return ur; }
        inline Point center() const { return Point((ll.x + ur.x) / 2, (ll.y + ur.y) / 2); }

//        bool inside(Point& p) const { return ll.x <= p.x && ll.y <= p.y && ur.x >= p.x && ur.y >= p.y; }

        inline bool operator==(const Box& p)
        {
            return this->ll == p.ll && this->ur == p.ur;
        }

        std::vector<Box> refine() const;

        std::vector<Point> corners() const;

        void translate(Real a);

        // return minimal and maximal value of func
        // on the corners of the box
        template<typename F>
        std::pair<Real, Real> min_max_on_corners(const F& func) const;

        friend std::ostream& operator<<(std::ostream& os, const Box& box);
    };

    std::ostream& operator<<(std::ostream& os, const Box& box);
//    template<typename InputIterator>
//    Box compute_bounding_box(InputIterator simplices_begin, InputIterator simplices_end)
//    {
//        if (simplices_begin == simplices_end) {
//            return Box();
//        }
//        Box bb;
//        bb.ll = bb.ur = simplices_begin->pos;
//        for (InputIterator it = simplices_begin; it != simplices_end; it++) {
//            Point& pos = it->pos;
//            if (pos.x < bb.ll.x) {
//                bb.ll.x = pos.x;
//            }
//            if (pos.y < bb.ll.y) {
//                bb.ll.y = pos.y;
//            }
//            if (pos.x > bb.ur.x) {
//                bb.ur.x = pos.x;
//            }
//            if (pos.y > bb.ur.y) {
//                bb.ur.y = pos.y;
//            }
//        }
//        return bb;
//    }

    Box get_enclosing_box(const Box& box_a, const Box& box_b);

    template<typename F>
    std::pair<Real, Real> Box::min_max_on_corners(const F& func) const
    {
        std::pair<Real, Real> min_max { std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max() };
        for(Point p : corners()) {
            Real value = func(p);
            min_max.first = std::min(min_max.first, value);
            min_max.second = std::max(min_max.second, value);
        }
        return min_max;
    };
} // namespace md

#endif //MATCHING_DISTANCE_BOX_H
