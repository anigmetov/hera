#include <tuple>

#include "dual_point.h"

namespace md {

    std::ostream& operator<<(std::ostream& os, const AxisType& at)
    {
        if (at == AxisType::x_type)
            os << "x-type";
        else
            os << "y-type";
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const AngleType& at)
    {
        if (at == AngleType::flat)
            os << "flat";
        else
            os << "steep";
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const DualPoint& dp)
    {
        os << "Line(" << dp.axis_type() << ", ";
        os << dp.angle_type() << ", ";
        os << dp.lambda() << ", ";
        os << dp.mu() << ", equation: ";
        if (not dp.is_vertical()) {
            os << "y = " << dp.y_slope() << " x + " << dp.y_intercept();
        } else {
            os << "x = " << dp.x_intercept();
        }
        os << " )";
        return os;
    }

    bool DualPoint::operator<(const DualPoint& rhs) const
    {
        return std::tie(axis_type_, angle_type_, lambda_, mu_)
                < std::tie(rhs.axis_type_, rhs.angle_type_, rhs.lambda_, rhs.mu_);
    }

    DualPoint::DualPoint(AxisType axis_type, AngleType angle_type, Real lambda, Real mu)
            :
            axis_type_(axis_type),
            angle_type_(angle_type),
            lambda_(lambda),
            mu_(mu)
    {
        assert(sanity_check());
    }

    bool DualPoint::sanity_check() const
    {
        if (lambda_ < 0.0)
            throw std::runtime_error("Invalid line, negative lambda");
        if (lambda_ > 1.0)
            throw std::runtime_error("Invalid line, lambda > 1");
        if (mu_ < 0.0)
            throw std::runtime_error("Invalid line, negative mu");
        return true;
    }

    Real DualPoint::gamma() const
    {
        if (is_steep())
            return atan(Real(1.0) / lambda_);
        else
            return atan(lambda_);
    }

    DualPoint midpoint(DualPoint x, DualPoint y)
    {
        assert(x.angle_type() == y.angle_type() and x.axis_type() == y.axis_type());
        Real lambda_mid = (x.lambda() + y.lambda()) / 2;
        Real mu_mid = (x.mu() + y.mu()) / 2;
        return DualPoint(x.axis_type(), x.angle_type(), lambda_mid, mu_mid);

    }

    // return k in the line equation y = kx + b
    Real DualPoint::y_slope() const
    {
        if (is_flat())
            return lambda();
        else
            return Real(1.0) / lambda();
    }

    // return k in the line equation x = ky + b
    Real DualPoint::x_slope() const
    {
        if (is_flat())
            return Real(1.0) / lambda();
        else
            return lambda();
    }

    // return b in the line equation y = kx + b
    Real DualPoint::y_intercept() const
    {
        if (is_y_type()) {
            return mu();
        } else {
            // x = x_slope * y + mu = x_slope * (y + mu / x_slope)
            // x-intercept is -mu/x_slope = -mu * y_slope
            return -mu() * y_slope();
        }
    }

    // return k in the line equation x = ky + b
    Real DualPoint::x_intercept() const
    {
        if (is_x_type()) {
            return mu();
        } else {
            // y = y_slope * x + mu = y_slope (x + mu / y_slope)
            // x_intercept is -mu/y_slope = -mu * x_slope
            return -mu() * x_slope();
        }
    }

    Real DualPoint::x_from_y(Real y) const
    {
        if (is_horizontal())
            throw std::runtime_error("x_from_y called on horizontal line");
        else
            return x_slope() * y + x_intercept();
    }

    Real DualPoint::y_from_x(Real x) const
    {
        if (is_vertical())
            throw std::runtime_error("x_from_y called on horizontal line");
        else
            return y_slope() * x + y_intercept();
    }

    bool DualPoint::is_horizontal() const
    {
        return is_flat() and lambda() == 0;
    }

    bool DualPoint::is_vertical() const
    {
        return is_steep() and lambda() == 0;
    }
    
    bool DualPoint::contains(Point p) const
    {
        if (is_vertical())
            return p.x == x_from_y(p.y);
        else
            return p.y == y_from_x(p.x);
    }

    bool DualPoint::goes_below(Point p) const
    {
        if (is_vertical())
            return p.x <= x_from_y(p.y);
        else
            return p.y >= y_from_x(p.x);
    }

    bool DualPoint::goes_above(Point p) const
    {
        if (is_vertical())
            return p.x >= x_from_y(p.y);
        else
            return p.y <= y_from_x(p.x);
    }

    Point DualPoint::push(Point p) const
    {
        Point result;
        // if line is below p, we push horizontally
        bool horizontal_push = goes_below(p);
        if (is_x_type()) {
            if (is_flat()) {
                if (horizontal_push) {
                    result.x = p.y / lambda() + mu();
                    result.y = p.y;
                } else {
                    // vertical push
                    result.x = p.x;
                    result.y = lambda() * (p.x - mu());
                }
            } else {
                // steep
                if (horizontal_push) {
                    result.x = lambda() * p.y + mu();
                    result.y = p.y;
                } else {
                    // vertical push
                    result.x = p.x;
                    result.y = (p.x - mu()) / lambda();
                }
            }
        } else {
            // y-type
            if (is_flat()) {
                if (horizontal_push) {
                    result.x = (p.y - mu()) / lambda();
                    result.y = p.y;
                } else {
                    // vertical push
                    result.x = p.x;
                    result.y = lambda() * p.x + mu();
                }
            } else {
                // steep
                if (horizontal_push) {
                    result.x = (p.y - mu()) * lambda();
                    result.y = p.y;
                } else {
                    // vertical push
                    result.x = p.x;
                    result.y = p.x / lambda() + mu();
                }
            }
        }
        return result;
    }

    Real DualPoint::weighted_push(Point p) const
    {
        // if line is below p, we push horizontally
        bool horizontal_push = goes_below(p);
        if (is_x_type()) {
            if (is_flat()) {
                if (horizontal_push) {
                    return p.y;
                } else {
                    // vertical push
                    return lambda() * (p.x - mu());
                }
            } else {
                // steep
                if (horizontal_push) {
                    return lambda() * p.y;
                } else {
                    // vertical push
                    return (p.x - mu());
                }
            }
        } else {
            // y-type
            if (is_flat()) {
                if (horizontal_push) {
                    return p.y - mu();
                } else {
                    // vertical push
                    return lambda() * p.x;
                }
            } else {
                // steep
                if (horizontal_push) {
                    return lambda() * (p.y - mu());
                } else {
                    // vertical push
                    return p.x;
                }
            }
        }
    }

    bool DualPoint::operator==(const DualPoint& other) const
    {
        return axis_type() == other.axis_type() and
               angle_type() == other.angle_type() and
               mu() == other.mu() and
               lambda() == other.lambda();
    }

    Real DualPoint::weight() const
    {
        return lambda_ / sqrt(1 + lambda_ * lambda_);
    }
} // namespace md
