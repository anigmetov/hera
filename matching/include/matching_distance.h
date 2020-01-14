#pragma once

#include <vector>
#include <limits>
#include <utility>
#include <ostream>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

#include "common_defs.h"
#include "cell_with_value.h"
#include "box.h"
#include "dual_point.h"
#include "dual_box.h"
#include "persistence_module.h"
#include "bifiltration.h"
#include "bottleneck.h"

namespace spd = spdlog;

namespace md {

    using HeatMap = std::map<DualPoint, Real>;
    using HeatMaps = std::map<int, HeatMap>;

    enum class BoundStrategy {
        bruteforce,
        local_dual_bound,
        local_dual_bound_refined,
        local_dual_bound_for_each_point,
        local_combined
    };

    enum class TraverseStrategy {
        depth_first,
        breadth_first,
        breadth_first_value,
        upper_bound
    };

    std::ostream& operator<<(std::ostream& os, const BoundStrategy& s);

    std::ostream& operator<<(std::ostream& os, const TraverseStrategy& s);

    std::istream& operator>>(std::istream& is, BoundStrategy& s);

    std::istream& operator>>(std::istream& is, TraverseStrategy& s);

    BoundStrategy bs_from_string(std::string s);

    TraverseStrategy ts_from_string(std::string s);

    struct CalculationParams {
        static constexpr int ALL_DIMENSIONS = -1;

        Real hera_epsilon {0.001}; // relative error in hera call
        Real delta {0.1}; // relative error for matching distance
        int max_depth {6}; // maximal number of refinenemnts
        int initialization_depth {3};
        int dim {0}; // in which dim to calculate the distance; use ALL_DIMENSIONS to get max over all dims
        BoundStrategy bound_strategy {BoundStrategy::bruteforce};
        TraverseStrategy traverse_strategy {TraverseStrategy::breadth_first};
        bool tolerate_max_iter_exceeded {true};
        Real actual_error {std::numeric_limits<Real>::max()};
        int actual_max_depth {0};
        int n_hera_calls {0};  // for experiments only; is set in matching_distance function, input value is ignored

        // stop looping over points immediately, if current point's displacement is too large
        // to prune the cell
        // if true, cells are pruned immediately, and bounds may be unreliable
        // (we just return something large enough to prune the cell)
        bool stop_asap { true };

#ifdef PRINT_HEAT_MAP
        HeatMaps heat_maps;
#endif
    };

    template<class DiagramProvider>
    class DistanceCalculator {

//        using DiagramProvider = md::Bifiltration;
//        using DiagramProvider = md::ModulePresentation;
        using DualBox = md::DualBox;

        using CellValueVector = std::vector<CellWithValue>;

    public:
        DistanceCalculator(const DiagramProvider& a,
                const DiagramProvider& b,
                CalculationParams& params);

        Real distance();

        int get_hera_calls_number() const;
//    for tests - make everything public
//    private:
        DiagramProvider module_a_;
        DiagramProvider module_b_;

        CalculationParams& params_;

        int n_hera_calls_;
        std::map<int, int> n_hera_calls_per_level_;
        Real distance_;

        // if calculate_on_intermediate, then weighted distance
        // will be calculated on centers of each grid in between
        CellValueVector get_refined_grid(int init_depth, bool calculate_on_intermediate, bool calculate_on_last = true);

        CellValueVector get_initial_dual_grid(Real& lower_bound);

        void heatmap_in_dimension(int dim, int depth);

        Real get_max_x(int module) const;

        Real get_max_y(int module) const;

        void set_cell_central_value(CellWithValue& dual_cell);

        Real get_distance();

        Real get_distance_pq();

        // temporary, to try priority queue
        Real get_max_possible_value(const CellWithValue* first_cell_ptr, int n_cells);

        Real get_upper_bound(const CellWithValue& dual_cell, Real good_enough_upper_bound) const;

        Real get_single_dgm_bound(const CellWithValue& dual_cell, ValuePoint vp, int module,
                Real good_enough_value) const;

        // this bound depends only on dual box
        Real get_local_dual_bound(int module, const DualBox& dual_box) const;

        Real get_local_dual_bound(const DualBox& dual_box) const;

        // this bound depends only on dual box, is more accurate
        Real get_local_refined_bound(int module, const md::DualBox& dual_box) const;

        Real get_local_refined_bound(const md::DualBox& dual_box) const;

        Real get_good_enough_upper_bound(Real lower_bound) const;

        Real
        get_max_displacement_single_point(const CellWithValue& dual_cell, ValuePoint value_point, const Point& p) const;

        void check_upper_bound(const CellWithValue& dual_cell) const;

        Real distance_on_line(DualPoint line);
        Real distance_on_line_const(DualPoint line) const;

        Real current_error(Real lower_bound, Real upper_bound);
    };

    Real matching_distance(const Bifiltration& bif_a, const Bifiltration& bif_b, CalculationParams& params);

    Real matching_distance(const ModulePresentation& mod_a, const ModulePresentation& mod_b, CalculationParams& params);

    // for upper bound experiment
    struct UbExperimentRecord {
        Real error;
        Real lower_bound;
        Real upper_bound;
        CellWithValue cell;
        long long int time;
        long long int n_hera_calls;
    };

    std::ostream& operator<<(std::ostream& os, const UbExperimentRecord& r);

}

#include "matching_distance.hpp"
