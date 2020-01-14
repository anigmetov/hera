#include <random>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

namespace spd = spdlog;

#include "dual_box.h"

namespace md {

    std::ostream& operator<<(std::ostream& os, const DualBox& db)
    {
        os << "DualBox(" << db.lower_left_ << ", " << db.upper_right_ << ")";
        return os;
    }

    DualBox::DualBox(DualPoint ll, DualPoint ur)
            :lower_left_(ll), upper_right_(ur)
    {
    }

    std::vector<DualPoint> DualBox::corners() const
    {
        return {lower_left_,
                DualPoint(axis_type(), angle_type(), lower_left_.lambda(), upper_right_.mu()),
                upper_right_,
                DualPoint(axis_type(), angle_type(), upper_right_.lambda(), lower_left_.mu())};
    }

    std::vector<DualPoint> DualBox::push_change_points(const Point& p) const
    {
        std::vector<DualPoint> result;
        result.reserve(2);

        bool is_y_type = lower_left_.is_y_type();
        bool is_flat = lower_left_.is_flat();

        auto mu_from_lambda = [p, is_y_type, is_flat](Real lambda) {
            bool is_x_type = not is_y_type, is_steep = not is_flat;
            if (is_y_type and is_flat) {
                return p.y - lambda * p.x;
            } else if (is_y_type and is_steep) {
                return p.y - p.x / lambda;
            } else if (is_x_type and is_flat) {
                return p.x - p.y / lambda;
            } else if (is_x_type and is_steep) {
                return p.x - lambda * p.y;
            }
            // to shut up compiler warning
            return static_cast<Real>(1.0 / 0.0);
        };

        auto lambda_from_mu = [p, is_y_type, is_flat](Real mu) {
            bool is_x_type = not is_y_type, is_steep = not is_flat;
            if (is_y_type and is_flat) {
                return (p.y - mu) / p.x;
            } else if (is_y_type and is_steep) {
                return p.x / (p.y - mu);
            } else if (is_x_type and is_flat) {
                return p.y / (p.x - mu);
            } else if (is_x_type and is_steep) {
                return (p.x - mu) / p.y;
            }
            // to shut up compiler warning
            return static_cast<Real>(1.0 / 0.0);
        };

        // all inequalities below are strict: equality means it is a corner
        // and critical_points() returns corners anyway

        Real mu_intersect_min = mu_from_lambda(lambda_min());

        if (mu_min() < mu_intersect_min && mu_intersect_min < mu_max())
            result.emplace_back(axis_type(), angle_type(), lambda_min(), mu_intersect_min);

        Real mu_intersect_max = mu_from_lambda(lambda_max());

        if (mu_max() < mu_intersect_max && mu_intersect_max < mu_max())
            result.emplace_back(axis_type(), angle_type(), lambda_max(), mu_intersect_max);

        Real lambda_intersect_min = lambda_from_mu(mu_min());

        if (lambda_min() < lambda_intersect_min && lambda_intersect_min < lambda_max())
            result.emplace_back(axis_type(), angle_type(), lambda_intersect_min, mu_min());

        Real lambda_intersect_max = lambda_from_mu(mu_max());
        if (lambda_min() < lambda_intersect_max && lambda_intersect_max < lambda_max())
            result.emplace_back(axis_type(), angle_type(), lambda_intersect_max, mu_max());

        assert(result.size() <= 2);

        if (result.size() > 2) {
            fmt::print("Error in push_change_points, p = {}, dual_box = {}, result = {}\n", p, *this,
                    container_to_string(result));
            throw std::runtime_error("push_change_points returned more than 2 points");
        }

        return result;
    }

    std::vector<DualPoint> DualBox::critical_points(const Point& /*p*/) const
    {
        // maximal difference is attained at corners
        return corners();
//        std::vector<DualPoint> result;
//        result.reserve(6);
//        for(auto dp : corners()) result.push_back(dp);
//        for(auto dp : push_change_points(p)) result.push_back(dp);
//        return result;
    }

    std::vector<DualPoint> DualBox::random_points(int n) const
    {
        assert(n >= 0);
        std::mt19937_64 gen(1);
        std::vector<DualPoint> result;
        result.reserve(n);
        std::uniform_real_distribution<Real> mu_distr(mu_min(), mu_max());
        std::uniform_real_distribution<Real> lambda_distr(lambda_min(), lambda_max());
        for(int i = 0; i < n; ++i) {
            result.emplace_back(axis_type(), angle_type(), lambda_distr(gen), mu_distr(gen));
        }
        return result;
    }

    bool DualBox::sanity_check() const
    {
        lower_left_.sanity_check();
        upper_right_.sanity_check();

        if (lower_left_.angle_type() != upper_right_.angle_type())
            throw std::runtime_error("angle types differ");

        if (lower_left_.axis_type() != upper_right_.axis_type())
            throw std::runtime_error("axis types differ");

        if (lower_left_.lambda() >= upper_right_.lambda())
            throw std::runtime_error("lambda of lower_left_ greater than lambda of upper_right ");

        if (lower_left_.mu() >= upper_right_.mu())
            throw std::runtime_error("mu of lower_left_ greater than mu of upper_right ");

        return true;
    }

    std::vector<DualBox> DualBox::refine() const
    {
        std::vector<DualBox> result;

        result.reserve(4);

        Real lambda_middle = (lower_left().lambda() + upper_right().lambda()) / 2.0;
        Real mu_middle = (lower_left().mu() + upper_right().mu()) / 2.0;

        DualPoint refinement_center(axis_type(), angle_type(), lambda_middle, mu_middle);

        result.emplace_back(lower_left_, refinement_center);

        result.emplace_back(DualPoint(axis_type(), angle_type(), lambda_middle, mu_min()),
                DualPoint(axis_type(), angle_type(), lambda_max(), mu_middle));

        result.emplace_back(refinement_center, upper_right_);

        result.emplace_back(DualPoint(axis_type(), angle_type(), lambda_min(), mu_middle),
                DualPoint(axis_type(), angle_type(), lambda_middle, mu_max()));
       return result;
    }

    bool DualBox::operator==(const DualBox& other) const
    {
        return lower_left() == other.lower_left() and
               upper_right() == other.upper_right();
    }

    bool DualBox::contains(const DualPoint& dp) const
    {
        return dp.angle_type() == angle_type() and dp.axis_type() == axis_type() and
               mu_max() >= dp.mu() and
               mu_min() <= dp.mu() and
               lambda_min() <= dp.lambda() and
               lambda_max() >= dp.lambda();
    }

    DualPoint DualBox::lower_right() const
    {
        return DualPoint(lower_left_.axis_type(), lower_left_.angle_type(), lambda_max(), mu_min());
    }

    DualPoint DualBox::upper_left() const
    {
        return DualPoint(lower_left_.axis_type(), lower_left_.angle_type(), lambda_min(), mu_max());
    }
}
