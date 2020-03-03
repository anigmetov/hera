#ifndef MATCHING_DISTANCE_CELL_WITH_VALUE_H
#define MATCHING_DISTANCE_CELL_WITH_VALUE_H

#include <algorithm>

#include "common_defs.h"
#include "common_util.h"
#include "dual_box.h"

namespace md {

    enum class ValuePoint {
        center,
        lower_left,
        lower_right,
        upper_left,
        upper_right
    };

    inline std::ostream& operator<<(std::ostream& os, const ValuePoint& vp)
    {
        switch(vp) {
            case ValuePoint::upper_left :
                os << "upper_left";
                break;
            case ValuePoint::upper_right :
                os << "upper_right";
                break;
            case ValuePoint::lower_left :
                os << "lower_left";
                break;
            case ValuePoint::lower_right :
                os << "lower_right";
                break;
            case ValuePoint::center:
                os << "center";
                break;
            default:
                os << "FORGOTTEN ValuePoint";
        }
        return os;
    }

    const std::vector<ValuePoint> k_all_vps = {ValuePoint::center, ValuePoint::lower_left, ValuePoint::upper_left,
                       ValuePoint::upper_right, ValuePoint::lower_right};

    const std::vector<ValuePoint> k_corner_vps = {ValuePoint::lower_left, ValuePoint::upper_left,
                       ValuePoint::upper_right, ValuePoint::lower_right};

    // represents a cell in the dual space with the value
    // of the weighted bottleneck distance
    template<class Real_>
    class CellWithValue {
    public:
        using Real = Real_;

        CellWithValue() = default;

        CellWithValue(const CellWithValue&) = default;

        CellWithValue(CellWithValue&&) = default;

        CellWithValue& operator=(const CellWithValue& other)& = default;

        CellWithValue& operator=(CellWithValue&& other) = default;

        CellWithValue(const DualBox<Real>& b, int level)
                :dual_box_(b), level_(level) { }

        DualBox<Real> dual_box() const { return dual_box_; }

        DualPoint<Real> center() const { return dual_box_.center(); }

        Real value_at(ValuePoint vp) const;

        bool has_value_at(ValuePoint vp) const;

        DualPoint<Real> value_point(ValuePoint vp) const;

        int level() const { return level_; }

        void set_value_at(ValuePoint vp, Real new_value);

        bool has_corner_value() const;

        Real stored_upper_bound() const;

        Real max_corner_value() const;

        Real min_value() const;

        bool has_max_possible_value() const { return has_max_possible_value_; }

        std::vector<CellWithValue> get_refined_cells() const;

        void set_max_possible_value(Real new_upper_bound);

        int num_values() const;

#ifdef MD_DEBUG
        long long int id { 0 };

        static long long int max_id;

        std::vector<long long int> parent_ids;
#endif

    private:

        bool has_central_value() const { return central_value_ >= 0; }

        bool has_lower_left_value() const { return lower_left_value_ >= 0; }

        bool has_lower_right_value() const { return lower_right_value_ >= 0; }

        bool has_upper_left_value() const { return upper_left_value_ >= 0; }

        bool has_upper_right_value() const { return upper_right_value_ >= 0; }


        DualBox<Real> dual_box_;
        Real central_value_ {-1.0};
        Real lower_left_value_ {-1.0};
        Real lower_right_value_ {-1.0};
        Real upper_left_value_ {-1.0};
        Real upper_right_value_ {-1.0};

        Real max_possible_value_ {0.0};

        int level_ {0};

        bool has_max_possible_value_ {false};
    };

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const CellWithValue<Real>& cell);
} // namespace md

#include "cell_with_value.hpp"

#endif //MATCHING_DISTANCE_CELL_WITH_VALUE_H
