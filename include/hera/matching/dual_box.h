#ifndef MATCHING_DISTANCE_DUAL_BOX_H
#define MATCHING_DISTANCE_DUAL_BOX_H

#include <ostream>
#include <limits>
#include <vector>
#include <random>

#include "common_util.h"
#include "dual_point.h"

namespace md {


    template<class Real>
    class DualBox {
    public:

        DualBox(DualPoint<Real> ll, DualPoint<Real> ur);

        DualBox() = default;
        DualBox(const DualBox&) = default;
        DualBox(DualBox&&) = default;

        DualBox& operator=(const DualBox& other) & = default;
        DualBox& operator=(DualBox&& other) = default;


        DualPoint<Real> center() const { return midpoint(lower_left_, upper_right_); }
        DualPoint<Real> lower_left() const { return lower_left_; }
        DualPoint<Real> upper_right() const { return upper_right_; }

        DualPoint<Real> lower_right() const;
        DualPoint<Real> upper_left() const;

        AxisType axis_type() const { return lower_left_.axis_type(); }
        AngleType angle_type() const { return lower_left_.angle_type(); }

        Real mu_min() const { return lower_left_.mu(); }
        Real mu_max() const { return upper_right_.mu(); }
        Real lambda_min() const { return lower_left_.lambda(); }
        Real lambda_max() const { return upper_right_.lambda(); }

        // return true, if all lines in dual_box are flat
        bool is_flat() const { return upper_right_.is_flat(); }
        bool is_steep() const { return lower_left_.is_steep(); }

        std::vector<DualBox> refine() const;
        std::vector<DualPoint<Real>> corners() const;
        std::vector<DualPoint<Real>> critical_points(const Point<Real>& p) const;
        // sample n points from the box uniformly; for tests
        std::vector<DualPoint<Real>> random_points(int n) const;

        // return 2 dual points at the boundary
        // where push changes from horizontal to vertical
        std::vector<DualPoint<Real>> push_change_points(const Point<Real>& p) const;

        // check that a has same sign, angles are all flat or all steep
        bool sanity_check() const;
        bool contains(const DualPoint<Real>& dp) const;

        bool operator==(const DualBox& other) const;

    private:
        DualPoint<Real> lower_left_;
        DualPoint<Real> upper_right_;
    };

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const DualBox<Real>& db)
    {
        os << "DualBox(" << db.lower_left() << ", " << db.upper_right() << ")";
        return os;
    }
}

#include "dual_box.hpp"

#endif //MATCHING_DISTANCE_DUAL_BOX_H
