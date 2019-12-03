//
// Created by narn on 12.02.19.
//

#ifndef MATCHING_DISTANCE_DUAL_POINT_H
#define MATCHING_DISTANCE_DUAL_POINT_H

#include <vector>
#include <ostream>

#include "common_util.h"
#include "box.h"

namespace md {

    enum class AxisType {
        x_type, y_type
    };
    enum class AngleType {
        flat, steep
    };

    // class has two flags of AxisType and AngleType.
    // ATTENTION. == operator is not overloaded,
    // so, e.g., line y = x has 4 different non-equal representation.
    // we are unlikely to ever need this, because 4 cases are
    // always treated separately.
    class DualPoint {
    public:
        using Real = md::Real;

        DualPoint() = default;

        DualPoint(const DualPoint&) = default;

        DualPoint(DualPoint&&) = default;

        DualPoint& operator=(const DualPoint& other)& = default;

        DualPoint& operator=(DualPoint&& other) = default;

        DualPoint(AxisType axis_type, AngleType angle_type,  Real lambda, Real mu);

        Real lambda() const { return lambda_; }

        Real mu() const { return mu_; }

        // angle between line and x-axis
        Real gamma() const;

        bool is_steep() const { return angle_type_ == AngleType::steep; }

        bool is_flat() const { return angle_type_ == AngleType::flat; }

        bool is_x_type() const { return axis_type_ == AxisType::x_type; }

        bool is_y_type() const { return axis_type_ == AxisType::y_type; }

        friend std::ostream& operator<<(std::ostream& os, const DualPoint& dp);
        bool operator<(const DualPoint& rhs) const;

        AxisType axis_type() const { return axis_type_; }
        AngleType angle_type() const { return angle_type_; }

        // throw exception, if fields are invalid
        // return true otherwise
        bool sanity_check() const;

        Real weighted_push(Point p) const;
        Point push(Point p) const;

        bool is_horizontal() const;
        bool is_vertical() const;

        bool goes_below(Point p) const;
        bool goes_above(Point p) const;

        bool contains(Point p) const;

        Real x_slope() const;
        Real y_slope() const;

        Real x_intercept() const;
        Real y_intercept() const;

        Real x_from_y(Real y) const;
        Real y_from_x(Real x) const;

        Real weight() const;

        bool operator==(const DualPoint& other) const;

    private:
        AxisType axis_type_ {AxisType::y_type};
        AngleType angle_type_ {AngleType::flat};
        // both initial values are invalid: lambda must be between 0 and 1
        Real lambda_ {-1.0};
        Real mu_ {-1.0};
    };

    std::ostream& operator<<(std::ostream& os, const DualPoint& dp);

    DualPoint midpoint(DualPoint x, DualPoint y);
};

#endif //MATCHING_DISTANCE_DUAL_POINT_H
