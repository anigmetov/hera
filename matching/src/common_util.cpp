#include <vector>
#include <utility>
#include <cmath>
#include <ostream>
#include <limits>
#include <algorithm>

#include <common_util.h>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

namespace md {


    int gcd(int a, int b)
    {
        assert(a != 0 or b != 0);
        // make b <= a
        std::tie(b, a) = std::minmax({ abs(a), abs(b) });
        if (b == 0)
            return a;
        while((a = a % b)) {
            std::swap(a, b);
        }
        return b;
    }

    int signum(int a)
    {
        if (a < 0)
            return -1;
        else if (a > 0)
            return 1;
        else
            return 0;
    }

    Rational reduce(Rational frac)
    {
        int d = gcd(frac.numerator, frac.denominator);
        frac.numerator /= d;
        frac.denominator /= d;
        return frac;
    }

    void Rational::reduce() { *this = md::reduce(*this); }


    Rational& Rational::operator*=(const md::Rational& rhs)
    {
        numerator *= rhs.numerator;
        denominator *= rhs.denominator;
        reduce();
        return *this;
    }

    Rational& Rational::operator/=(const md::Rational& rhs)
    {
        numerator *= rhs.denominator;
        denominator *= rhs.numerator;
        reduce();
        return *this;
    }

    Rational& Rational::operator+=(const md::Rational& rhs)
    {
        numerator = numerator * rhs.denominator + denominator * rhs.numerator;
        denominator *= rhs.denominator;
        reduce();
        return *this;
    }

    Rational& Rational::operator-=(const md::Rational& rhs)
    {
        numerator = numerator * rhs.denominator - denominator * rhs.numerator;
        denominator *= rhs.denominator;
        reduce();
        return *this;
    }


    Rational midpoint(Rational a, Rational b)
    {
        return reduce({a.numerator * b.denominator + a.denominator * b.numerator, 2 * a.denominator * b.denominator });
    }

    Rational operator+(Rational a, const Rational& b)
    {
        a += b;
        return a;
    }

    Rational operator-(Rational a, const Rational& b)
    {
        a -= b;
        return a;
    }

    Rational operator*(Rational a, const Rational& b)
    {
        a *= b;
        return a;
    }

    Rational operator/(Rational a, const Rational& b)
    {
        a /= b;
        return a;
    }

    bool is_less(Rational a, Rational b)
    {
        // compute a - b = a_1 / a_2 - b_1 / b_2
        long numer = a.numerator * b.denominator - a.denominator * b.numerator;
        long denom = a.denominator * b.denominator;
        assert(denom != 0);
        return signum(numer) * signum(denom) < 0;
    }

    bool operator==(const Rational& a, const Rational& b)
    {
        return std::tie(a.numerator, a.denominator) == std::tie(b.numerator, b.denominator);
    }

    bool operator<(const Rational& a, const Rational& b)
    {
        // do not remove signum - overflow
        long numer = a.numerator * b.denominator - a.denominator * b.numerator;
        long denom = a.denominator * b.denominator;
        assert(denom != 0);
//        spdlog::debug("a = {}, b = {}, numer = {}, denom = {}, result = {}", a, b, numer, denom, signum(numer) * signum(denom) <= 0);
        return signum(numer) * signum(denom) < 0;
     }

    bool is_leq(Rational a, Rational b)
    {
        // compute a - b = a_1 / a_2 - b_1 / b_2
        long numer = a.numerator * b.denominator - a.denominator * b.numerator;
        long denom = a.denominator * b.denominator;
        assert(denom != 0);
        return signum(numer) * signum(denom) <= 0;
    }

    bool is_greater(Rational a, Rational b)
    {
        return not is_leq(a, b);
    }

    bool is_geq(Rational a, Rational b)
    {
        return not is_less(a, b);
    }

    Point operator+(const Point& u, const Point& v)
    {
        return Point(u.x + v.x, u.y + v.y);
    }

    Point operator-(const Point& u, const Point& v)
    {
        return Point(u.x - v.x, u.y - v.y);
    }

    Point least_upper_bound(const Point& u, const Point& v)
    {
        return Point(std::max(u.x, v.x), std::max(u.y, v.y));
    }

    Point greatest_lower_bound(const Point& u, const Point& v)
    {
        return Point(std::min(u.x, v.x), std::min(u.y, v.y));
    }

    Point max_point()
    {
        return Point(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::min());
    }

    Point min_point()
    {
        return Point(-std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::min());
    }

    std::ostream& operator<<(std::ostream& ostr, const Point& vec)
    {
        ostr << "(" << vec.x << ", " << vec.y << ")";
        return ostr;
    }

    Real l_infty_norm(const Point& v)
    {
        return std::max(std::abs(v.x), std::abs(v.y));
    }

    Real l_2_norm(const Point& v)
    {
        return v.norm();
    }

    Real l_2_dist(const Point& x, const Point& y)
    {
        return l_2_norm(x - y);
    }

    Real l_infty_dist(const Point& x, const Point& y)
    {
        return l_infty_norm(x - y);
    }

    void DiagramKeeper::add_point(int dim, md::Real birth, md::Real death)
    {
        data_[dim].emplace_back(birth, death);
    }

    DiagramKeeper::Diagram DiagramKeeper::get_diagram(int dim) const
    {
        if (data_.count(dim) == 1)
            return data_.at(dim);
        else
            return DiagramKeeper::Diagram();
    }

    // return true, if line starts with #
    // or contains only spaces
    bool ignore_line(const std::string& s)
    {
        for(auto c : s) {
            if (isspace(c))
                continue;
            return (c == '#');
        }
        return true;
    }



    std::ostream& operator<<(std::ostream& os, const Rational& a)
    {
        os << a.numerator << " / " << a.denominator;
        return os;
    }
}
