#ifndef MATCHING_DISTANCE_COMMON_UTIL_H
#define MATCHING_DISTANCE_COMMON_UTIL_H

#include <cassert>
#include <vector>
#include <utility>
#include <cmath>
#include <ostream>
#include <sstream>
#include <string>
#include <map>
#include <functional>

#include "common_defs.h"
#include "phat/helpers/misc.h"

namespace md {

    using Index = phat::index;
    using IndexVec = std::vector<Index>;

    //static constexpr Real pi = M_PI;

    using Column = std::vector<Index>;

    template<class Real>
    struct Point {
        Real x;
        Real y;

        Point(Real x = 0, Real y = 0)
                :x(x), y(y) { }

        inline Real norm() const { return sqrt(x * x + y * y); }

        inline void normalize()
        {
            Real nor = norm();
            x /= nor;
            y /= nor;
        }

        inline void translate(Real a)
        {
            x += a;
            y += a;
        }

        inline bool operator==(const Point& v) const
        {
            return this->x == v.x && this->y == v.y;
        }

        // compare both coordinates, as needed in persistence
        // do not overload operator<, requirements are not satisfied
        inline bool is_less(const Point& other, bool strict = false) const
        {
            if (x <= other.x and y <= other.y) {
                if (strict) {
                    return x != other.x or y != other.y;
                }
                else {
                    return true;
                }
            }
            return false;
        }
    };


    template<class Real>
    using PointVec = std::vector<Point<Real>>;

    template<class Real>
    Point<Real> operator+(const Point<Real>& u, const Point<Real>& v);

    template<class Real>
    Point<Real> operator-(const Point<Real>& u, const Point<Real>& v);


    template<class Real>
    Point<Real> least_upper_bound(const Point<Real>& u, const Point<Real>& v);

    template<class Real>
    Point<Real> greatest_lower_bound(const Point<Real>& u, const Point<Real>& v);

    template<class Real>
    Point<Real> max_point();

    template<class Real>
    Point<Real> min_point();

    template<class Real>
    std::ostream& operator<<(std::ostream& ostr, const Point<Real>& vec);

    template<class Real>
    using DiagramPoint = std::pair<Real, Real>;

    template<class Real>
    using Diagram = std::vector<DiagramPoint<Real>>;


    // to keep diagrams in all dimensions
    // TODO: store in Hera format?
    template<class Real>
    class DiagramKeeper {
    public:

        DiagramKeeper() { };

        void add_point(int dim, Real birth, Real death);

        Diagram<Real> get_diagram(int dim) const;

        void clear() { data_.clear(); }

    private:
        std::map<int, Diagram<Real>> data_;
    };

    template<typename C>
    std::string container_to_string(const C& cont)
    {
        std::stringstream ss;
        ss << "[";
        int i = 0;
        for (const auto& x : cont) {
            i++;
            ss << x;
            if (i != (int) cont.size())
                ss << ", ";
        }
        ss << "]";
        return ss.str();
    }

    // return true, if s is empty or starts with # (commented out line)
    // whitespaces in the beginning of s are ignored
    inline bool ignore_line(const std::string& s)
    {
        for(auto c : s) {
            if (isspace(c))
                continue;
            return (c == '#');
        }
        return true;
    }


    // split string by delimeter
    template<typename Out>
    void split_by_delim(const std::string &s, char delim, Out result) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            *(result++) = item;
        }
    }

    inline std::vector<std::string> split_by_delim(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split_by_delim(s, delim, std::back_inserter(elems));
        return elems;
    }
}

namespace std {
    template<class Real>
    struct hash<md::Point<Real>>
    {
        std::size_t operator()(const md::Point<Real>& p) const
        {
            auto hx = std::hash<decltype(p.x)>()(p.x);
            auto hy = std::hash<decltype(p.y)>()(p.y);
            return hx ^ (hy << 1);
        }
    };
};

#include "common_util.hpp"


#endif //MATCHING_DISTANCE_COMMON_UTIL_H
