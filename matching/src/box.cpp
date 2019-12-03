
#include "box.h"

namespace md {

    std::ostream& operator<<(std::ostream& os, const Box& box)
    {
        os << "Box(lower_left = " << box.lower_left() << ", upper_right = " << box.upper_right() << ")";
        return os;
    }

    Box get_enclosing_box(const Box& box_a, const Box& box_b)
    {
        Point lower_left(std::min(box_a.lower_left().x, box_b.lower_left().x),
                std::min(box_a.lower_left().y, box_b.lower_left().y));
        Point upper_right(std::max(box_a.upper_right().x, box_b.upper_right().x),
                std::max(box_a.upper_right().y, box_b.upper_right().y));
        return Box(lower_left, upper_right);
    }

    void Box::translate(md::Real a)
    {
        ll.x += a;
        ll.y += a;
        ur.x += a;
        ur.y += a;
    }

    std::vector<Box> Box::refine() const
    {
        std::vector<Box> result;

//        1 | 2
//        0 | 3

        Point new_ll = lower_left();
        Point new_ur = center();
        result.emplace_back(new_ll, new_ur);

        new_ll.y = center().y;
        new_ur.y = ur.y;
        result.emplace_back(new_ll, new_ur);

        new_ll = center();
        new_ur = upper_right();
        result.emplace_back(new_ll, new_ur);

        new_ll.y = ll.y;
        new_ur.y = center().y;
        result.emplace_back(new_ll, new_ur);

        return result;
    }

    std::vector<Point> Box::corners() const
    {
        return {ll, Point(ll.x, ur.y), ur, Point(ur.x, ll.y)};
    };


}
