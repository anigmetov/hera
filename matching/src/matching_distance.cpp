#include <chrono>
#include <tuple>
#include <algorithm>

#include "common_defs.h"

#include "matching_distance.h"

namespace md {

    Real matching_distance(const Bifiltration& bif_a, const Bifiltration& bif_b,
            CalculationParams& params)
    {
        Real result;
        // compute distance only in one dimension
        if (params.dim != CalculationParams::ALL_DIMENSIONS) {
            BifiltrationProxy bifp_a(bif_a, params.dim);
            BifiltrationProxy bifp_b(bif_b, params.dim);
            DistanceCalculator<BifiltrationProxy> runner(bifp_a, bifp_b, params);
            result = runner.distance();
            params.n_hera_calls = runner.get_hera_calls_number();
        } else {
            // compute distance in all dimensions, return maximal
            result = -1;
            for(int dim = 0; dim < std::max(bif_a.maximal_dim(), bif_b.maximal_dim()); ++dim) {
                BifiltrationProxy bifp_a(bif_a, params.dim);
                BifiltrationProxy bifp_b(bif_a, params.dim);
                DistanceCalculator<BifiltrationProxy> runner(bifp_a, bifp_b, params);
                result = std::max(result, runner.distance());
                params.n_hera_calls += runner.get_hera_calls_number();
            }
        }
        return result;
    }


    Real matching_distance(const ModulePresentation& mod_a, const ModulePresentation& mod_b,
            CalculationParams& params)
    {
        DistanceCalculator<ModulePresentation> runner(mod_a, mod_b, params);
        Real result = runner.distance();
        params.n_hera_calls = runner.get_hera_calls_number();
        return result;
    }

    std::istream& operator>>(std::istream& is, BoundStrategy& s)
    {
        std::string ss;
        is >> ss;
        if (ss == "bruteforce") {
            s = BoundStrategy::bruteforce;
        } else if (ss == "local_grob") {
            s = BoundStrategy::local_dual_bound;
        } else if (ss == "local_combined") {
            s = BoundStrategy::local_combined;
        } else if (ss == "local_refined") {
            s = BoundStrategy::local_dual_bound_refined;
        } else if (ss == "local_for_each_point") {
            s = BoundStrategy::local_dual_bound_for_each_point;
        } else {
            throw std::runtime_error("UNKNOWN BOUND STRATEGY");
        }
        return is;
    }

    BoundStrategy bs_from_string(std::string s)
    {
        std::stringstream ss(s);
        BoundStrategy result;
        ss >> result;
        return result;
    }

    TraverseStrategy ts_from_string(std::string s)
    {
        std::stringstream ss(s);
        TraverseStrategy result;
        ss >> result;
        return result;
    }

    std::istream& operator>>(std::istream& is, TraverseStrategy& s)
    {
        std::string ss;
        is >> ss;
        if (ss == "DFS") {
            s = TraverseStrategy::depth_first;
        } else if (ss == "BFS") {
            s = TraverseStrategy::breadth_first;
        } else if (ss == "BFS-VAL") {
            s = TraverseStrategy::breadth_first_value;
        } else if (ss == "UB") {
            s = TraverseStrategy::upper_bound;
        } else {
            throw std::runtime_error("UNKNOWN TRAVERSE STRATEGY");
        }
        return is;
    }

    std::ostream& operator<<(std::ostream& os, const UbExperimentRecord& r)
    {
        os << r.time << "\t" << r.n_hera_calls << "\t" << r.error << "\t" << r.lower_bound << "\t" << r.upper_bound;
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const BoundStrategy& s)
    {
        switch(s) {
            case BoundStrategy::bruteforce :
                os << "bruteforce";
                break;
            case BoundStrategy::local_dual_bound :
                os << "local_grob";
                break;
            case BoundStrategy::local_combined :
                os << "local_combined";
                break;
            case BoundStrategy::local_dual_bound_refined :
                os << "local_refined";
                break;
            case BoundStrategy::local_dual_bound_for_each_point :
                os << "local_for_each_point";
                break;
            default:
                os << "FORGOTTEN BOUND STRATEGY";
        }
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const TraverseStrategy& s)
    {
        switch(s) {
            case TraverseStrategy::depth_first :
                os << "DFS";
                break;
            case TraverseStrategy::breadth_first :
                os << "BFS";
                break;
            case TraverseStrategy::breadth_first_value :
                os << "BFS-VAL";
                break;
            case TraverseStrategy::upper_bound :
                os << "UB";
                break;
            default:
                os << "FORGOTTEN TRAVERSE STRATEGY";
        }
        return os;
    }
}
