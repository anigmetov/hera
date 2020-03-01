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


    using Real = double;
    using RealVec = std::vector<Real>;
    using Index = phat::index;
    using IndexVec = std::vector<Index>;

    static constexpr Real pi = M_PI;

    using Column = std::vector<Index>;

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


    using PointVec = std::vector<Point>;

    Point operator+(const Point& u, const Point& v);

    Point operator-(const Point& u, const Point& v);

    Point least_upper_bound(const Point& u, const Point& v);

    Point greatest_lower_bound(const Point& u, const Point& v);

    Point max_point();

    Point min_point();

    std::ostream& operator<<(std::ostream& ostr, const Point& vec);

    Real L_infty(const Point& v);

    Real l_2_norm(const Point& v);

    Real l_2_dist(const Point& x, const Point& y);

    Real l_infty_dist(const Point& x, const Point& y);

    using Interval = std::pair<Real, Real>;

    // return minimal interval that contains both a and b
    inline Interval minimal_covering_interval(Interval a, Interval b)
    {
        return {std::min(a.first, b.first), std::max(a.second, b.second)};
    }

    // to keep diagrams in all dimensions
    // TODO: store in Hera format?
    class DiagramKeeper {
    public:
        using DiagramPoint = std::pair<Real, Real>;
        using Diagram = std::vector<DiagramPoint>;

        DiagramKeeper() { };

        void add_point(int dim, Real birth, Real death);

        Diagram get_diagram(int dim) const;

        void clear() { data_.clear(); }

    private:
        std::map<int, Diagram> data_;
    };

    using Diagram = std::vector<std::pair<Real, Real>>;

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

    int gcd(int a, int b);

    struct Rational {
        int numerator {0};
        int denominator {1};
        Rational() = default;
        Rational(int n, int d) : numerator(n / gcd(n, d)), denominator(d / gcd(n, d)) {}
        Rational(std::pair<int, int> p) : Rational(p.first, p.second) {}
        Rational(int n) : numerator(n), denominator(1) {}
        Real to_real() const { return (Real)numerator / (Real)denominator; }
        void reduce();
        Rational& operator+=(const Rational& rhs);
        Rational& operator-=(const Rational& rhs);
        Rational& operator*=(const Rational& rhs);
        Rational& operator/=(const Rational& rhs);
    };

    using namespace std::rel_ops;

    bool operator==(const Rational& a, const Rational& b);
    bool operator<(const Rational& a, const Rational& b);
    std::ostream& operator<<(std::ostream& os, const Rational& a);

    // arithmetic
    Rational operator+(Rational a, const Rational& b);
    Rational operator-(Rational a, const Rational& b);
    Rational operator*(Rational a, const Rational& b);
    Rational operator/(Rational a, const Rational& b);

    Rational reduce(Rational frac);

    Rational midpoint(Rational a, Rational b);

    // return true, if s is empty or starts with # (commented out line)
    // whitespaces in the beginning of s are ignored
    bool ignore_line(const std::string& s);

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
    template<>
    struct hash<md::Point>
    {
        std::size_t operator()(const md::Point& p) const
        {
            auto hx = std::hash<decltype(p.x)>()(p.x);
            auto hy = std::hash<decltype(p.y)>()(p.y);
            return hx ^ (hy << 1);
        }
    };
};


#endif //MATCHING_DISTANCE_COMMON_UTIL_H
