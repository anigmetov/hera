#ifndef MATCHING_DISTANCE_DUAL_BOX_H
#define MATCHING_DISTANCE_DUAL_BOX_H

#include <ostream>
#include <limits>
#include <vector>

#include "common_util.h"
#include "dual_point.h"

namespace md {

    class DualBox {
    public:

        DualBox(DualPoint ll, DualPoint ur);

        DualBox() = default;
        DualBox(const DualBox&) = default;
        DualBox(DualBox&&) = default;

        DualBox& operator=(const DualBox& other) & = default;
        DualBox& operator=(DualBox&& other) = default;


        DualPoint center() const { return midpoint(lower_left_, upper_right_); }
        DualPoint lower_left() const { return lower_left_; }
        DualPoint upper_right() const { return upper_right_; }

        DualPoint lower_right() const;
        DualPoint upper_left() const;

        AxisType axis_type() const { return lower_left_.axis_type(); }
        AngleType angle_type() const { return lower_left_.angle_type(); }

        Real mu_min() const { return lower_left_.mu(); }
        Real mu_max() const { return upper_right_.mu(); }
        Real lambda_min() const { return lower_left_.lambda(); }
        Real lambda_max() const { return upper_right_.lambda(); }

        // return true, if all lines in dual_box are flat
        bool is_flat() const { return upper_right_.is_flat(); }
        bool is_steep() const { return lower_left_.is_steep(); }

        // return minimal and maximal value of func
        // on the corners of the box
        template<typename F>
        std::pair<Real, Real> min_max_on_corners(const F& func) const;

        template<typename F>
        Real max_abs_value(const F& func) const;


        std::vector<DualBox> refine() const;
        std::vector<DualPoint> corners() const;
        std::vector<DualPoint> critical_points(const Point& p) const;
        // sample n points from the box uniformly; for tests
        std::vector<DualPoint> random_points(int n) const;

        // return 2 dual points at the boundary
        // where push changes from horizontal to vertical
        std::vector<DualPoint> push_change_points(const Point& p) const;

        friend std::ostream& operator<<(std::ostream& os, const DualBox& db);

        // check that a has same sign, angles are all flat or all steep
        bool sanity_check() const;
        bool contains(const DualPoint& dp) const;

        bool operator==(const DualBox& other) const;

    private:
        DualPoint lower_left_;
        DualPoint upper_right_;
    };

    std::ostream& operator<<(std::ostream& os, const DualBox& db);

    template<typename F>
    std::pair<Real, Real> DualBox::min_max_on_corners(const F& func) const
    {
        std::pair<Real, Real> min_max { std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max() };
        for(auto p : corners()) {
            Real value = func(p);
            min_max.first = std::min(min_max.first, value);
            min_max.second = std::max(min_max.second, value);
        }
        return min_max;
    };


    template<typename F>
    Real DualBox::max_abs_value(const F& func) const
    {
        Real result = 0;
        for(auto p_1 : corners()) {
            for(auto p_2 : corners()) {
                Real value = fabs(func(p_1, p_2));
                result = std::max(value, result);
            }
        }
        return result;
    };

}

#endif //MATCHING_DISTANCE_DUAL_BOX_H
