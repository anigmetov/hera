namespace md {

#ifdef MD_DEBUG
    long long int CellWithValue<Real>::max_id = 0;
#endif

    template<class Real>
    Real CellWithValue<Real>::value_at(ValuePoint vp) const
    {
        switch(vp) {
            case ValuePoint::upper_left :
                return upper_left_value_;
            case ValuePoint::upper_right :
                return upper_right_value_;
            case ValuePoint::lower_left :
                return lower_left_value_;
            case ValuePoint::lower_right :
                return lower_right_value_;
            case ValuePoint::center:
                return central_value_;
        }
        // to shut up compiler warning
        return 1.0 / 0.0;
    }

    template<class Real>
    bool CellWithValue<Real>::has_value_at(ValuePoint vp) const
    {
        switch(vp) {
            case ValuePoint::upper_left :
                return upper_left_value_ >= 0;
            case ValuePoint::upper_right :
                return upper_right_value_ >= 0;
            case ValuePoint::lower_left :
                return lower_left_value_ >= 0;
            case ValuePoint::lower_right :
                return lower_right_value_ >= 0;
            case ValuePoint::center:
                return central_value_ >= 0;
        }
        // to shut up compiler warning
        return 1.0 / 0.0;
    }

    template<class Real>
    DualPoint<Real> CellWithValue<Real>::value_point(md::ValuePoint vp) const
    {
        switch(vp) {
            case ValuePoint::upper_left :
                return dual_box().upper_left();
            case ValuePoint::upper_right :
                return dual_box().upper_right();
            case ValuePoint::lower_left :
                return dual_box().lower_left();
            case ValuePoint::lower_right :
                return dual_box().lower_right();
            case ValuePoint::center:
                return dual_box().center();
        }
        // to shut up compiler warning
        return DualPoint<Real>();
    }

    template<class Real>
    bool CellWithValue<Real>::has_corner_value() const
    {
        return has_lower_left_value() or has_lower_right_value() or has_upper_left_value()
                or has_upper_right_value();
    }

    template<class Real>
    Real CellWithValue<Real>::stored_upper_bound() const
    {
        assert(has_max_possible_value_);
        return max_possible_value_;
    }

    template<class Real>
    Real CellWithValue<Real>::max_corner_value() const
    {
        return std::max({lower_left_value_, lower_right_value_, upper_left_value_, upper_right_value_});
    }

    template<class Real>
    Real CellWithValue<Real>::min_value() const
    {
        Real result = std::numeric_limits<Real>::max();
        for(auto vp : k_all_vps) {
            if (not has_value_at(vp)) {
                continue;
            }
            result = std::min(result, value_at(vp));
        }
        return result;
    }

    template<class Real>
    std::vector<CellWithValue<Real>> CellWithValue<Real>::get_refined_cells() const
    {
        std::vector<CellWithValue<Real>> result;
        result.reserve(4);
        for(const auto& refined_box : dual_box_.refine()) {

            CellWithValue<Real> refined_cell(refined_box, level() + 1);

#ifdef MD_DEBUG
            refined_cell.parent_ids = parent_ids;
            refined_cell.parent_ids.push_back(id);
            refined_cell.id = ++max_id;
#endif

            if (refined_box.lower_left() == dual_box_.lower_left()) {
                // _|_
                // H|_

                refined_cell.set_value_at(ValuePoint::lower_left, lower_left_value_);
                refined_cell.set_value_at(ValuePoint::upper_right, central_value_);

            } else if (refined_box.upper_right() == dual_box_.upper_right()) {
                // _|H
                // _|_

                refined_cell.set_value_at(ValuePoint::lower_left, central_value_);
                refined_cell.set_value_at(ValuePoint::upper_right, upper_right_value_);

            } else if (refined_box.lower_right() == dual_box_.lower_right()) {
                // _|_
                // _|H

                refined_cell.set_value_at(ValuePoint::lower_right, lower_right_value_);
                refined_cell.set_value_at(ValuePoint::upper_left, central_value_);

            } else if (refined_box.upper_left() == dual_box_.upper_left()) {

                // H|_
                // _|_

                refined_cell.set_value_at(ValuePoint::lower_right, central_value_);
                refined_cell.set_value_at(ValuePoint::upper_left, upper_left_value_);
            }
            result.emplace_back(refined_cell);
        }
        return result;
    }

    template<class Real>
    void CellWithValue<Real>::set_value_at(ValuePoint vp, Real new_value)
    {
        if (has_value_at(vp)) {
            std::cerr << "CellWithValue<Real>: trying to re-assign value, vp = " << vp;
        }

        switch(vp) {
            case ValuePoint::upper_left :
                upper_left_value_ = new_value;
                break;
            case ValuePoint::upper_right :
                upper_right_value_ = new_value;
                break;
            case ValuePoint::lower_left :
                lower_left_value_ = new_value;
                break;
            case ValuePoint::lower_right :
                lower_right_value_ = new_value;
                break;
            case ValuePoint::center:
                central_value_ = new_value;
                break;
        }
    }

    template<class Real>
    int CellWithValue<Real>::num_values() const
    {
        int result = 0;
        for(ValuePoint vp : k_all_vps) {
            result += has_value_at(vp);
        }
        return result;
    }


    template<class Real>
    void CellWithValue<Real>::set_max_possible_value(Real new_upper_bound)
    {
        assert(new_upper_bound >= central_value_);
        assert(new_upper_bound >= lower_left_value_);
        assert(new_upper_bound >= lower_right_value_);
        assert(new_upper_bound >= upper_left_value_);
        assert(new_upper_bound >= upper_right_value_);
        has_max_possible_value_ = true;
        max_possible_value_ = new_upper_bound;
    }


    template<class Real>
    std::ostream& operator<<(std::ostream& os, const CellWithValue<Real>& cell)
    {
        os << "CellWithValue(box = " << cell.dual_box() << ", ";

#ifdef MD_DEBUG
        os << "id = " << cell.id;
        if (not cell.parent_ids.empty())
            os << ", parent_ids = " << container_to_string(cell.parent_ids) << ", ";
#endif

        for(ValuePoint vp : k_all_vps) {
            if (cell.has_value_at(vp)) {
                os << "value = " << cell.value_at(vp);
                os << ", at " << vp << " " << cell.value_point(vp);
            }
        }

        os << ", max_corner_value = ";
        if (cell.has_max_possible_value()) {
            os << cell.stored_upper_bound();
        } else {
            os << "-";
        }

        os << ", level = " << cell.level() << ")";
        return os;
    }

} // namespace md
