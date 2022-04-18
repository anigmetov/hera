#include <vector>
#include <utility>
#include <cmath>
#include <ostream>
#include <limits>
#include <algorithm>

#include <common_util.h>

namespace md {

    template<class Real>
    Point<Real> operator+(const Point<Real>& u, const Point<Real>& v)
    {
        return Point<Real>(u.x + v.x, u.y + v.y);
    }

    template<class Real>
    Point<Real> operator-(const Point<Real>& u, const Point<Real>& v)
    {
        return Point<Real>(u.x - v.x, u.y - v.y);
    }

    template<class Real>
    Point<Real> least_upper_bound(const Point<Real>& u, const Point<Real>& v)
    {
        return Point<Real>(std::max(u.x, v.x), std::max(u.y, v.y));
    }

    template<class Real>
    Point<Real> greatest_lower_bound(const Point<Real>& u, const Point<Real>& v)
    {
        return Point<Real>(std::min(u.x, v.x), std::min(u.y, v.y));
    }

    template<class Real>
    Point<Real> max_point()
    {
        return Point<Real>(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::min());
    }

    template<class Real>
    Point<Real> min_point()
    {
        return Point<Real>(-std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::min());
    }

    template<class Real>
    std::ostream& operator<<(std::ostream& ostr, const Point<Real>& vec)
    {
        ostr << "(" << vec.x << ", " << vec.y << ")";
        return ostr;
    }

    template<class Real>
    Real l_infty_norm(const Point<Real>& v)
    {
        return std::max(std::abs(v.x), std::abs(v.y));
    }

    template<class Real>
    Real l_2_norm(const Point<Real>& v)
    {
        return v.norm();
    }

    template<class Real>
    Real l_2_dist(const Point<Real>& x, const Point<Real>& y)
    {
        return l_2_norm(x - y);
    }

    template<class Real>
    Real l_infty_dist(const Point<Real>& x, const Point<Real>& y)
    {
        return l_infty_norm(x - y);
    }

    template<class Real>
    void DiagramKeeper<Real>::add_point(int dim, Real birth, Real death)
    {
        data_[dim].emplace_back(birth, death);
    }

    template<class Real>
    Diagram<Real> DiagramKeeper<Real>::get_diagram(int dim) const
    {
        if (data_.count(dim) == 1)
            return data_.at(dim);
        else
            return Diagram<Real>();
    }
}
