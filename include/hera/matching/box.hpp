namespace md {

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const Box<Real>& box)
    {
        os << "Box(lower_left = " << box.lower_left() << ", upper_right = " << box.upper_right() << ")";
        return os;
    }

    template<class Real>
    void Box<Real>::translate(Real a)
    {
        ll.x += a;
        ll.y += a;
        ur.x += a;
        ur.y += a;
    }

    template<class Real>
    std::vector<Box<Real>> Box<Real>::refine() const
    {
        std::vector<Box<Real>> result;

//        1 | 2
//        0 | 3

        Point<Real> new_ll = lower_left();
        Point<Real> new_ur = center();
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

    template<class Real>
    std::vector<Point<Real>> Box<Real>::corners() const
    {
        return {ll, Point<Real>(ll.x, ur.y), ur, Point<Real>(ur.x, ll.y)};
    };

}
