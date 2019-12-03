#include "common_defs.h"

#include <iostream>
#include <string>
#include <cassert>
#include <experimental/filesystem>

#ifdef EXPERIMENTAL_TIMING
#include <chrono>
#endif

#include "opts/opts.h"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

//#include "persistence_module.h"
#include "bifiltration.h"
#include "box.h"
#include "matching_distance.h"

using namespace md;

namespace fs = std::experimental::filesystem;

void print_heat_map(const md::HeatMaps& hms, std::string fname, const CalculationParams& params)
{
#ifdef PRINT_HEAT_MAP
    spd::debug("Entered print_heat_map");
    std::set<Real> mu_vals, lambda_vals;
    auto hm_iter = hms.end();
    --hm_iter;
    int max_level = hm_iter->first;

    int level_cardinality = 4;
    for(int i = 0; i < params.initialization_depth; ++i) {
        level_cardinality *= 4;
    }
    for(int i = params.initialization_depth + 1; i <= max_level; ++i) {
        spd::debug("hms.at({}).size = {}, must be {}", i, hms.at(i).size(), level_cardinality);
        assert(static_cast<decltype(level_cardinality)>(hms.at(i).size()) == level_cardinality);
        level_cardinality *= 4;
    }

    std::map<std::pair<Real, Real>, Real> hm_x_flat, hm_x_steep, hm_y_flat, hm_y_steep;

    for(const auto& dual_point_value_pair : hms.at(max_level)) {
        const DualPoint& k = dual_point_value_pair.first;
        spd::debug("HM DP: {}", k);
        mu_vals.insert(k.mu());
        lambda_vals.insert(k.lambda());
    }

    std::vector<Real> lambda_vals_vec(lambda_vals.begin(), lambda_vals.end());
    std::vector<Real> mu_vals_vec(mu_vals.begin(), mu_vals.end());

    std::ofstream ofs {fname};
    if (not ofs.good()) {
        std::cerr << "Cannot write heat map to file " << fname << std::endl;
        throw std::runtime_error("Cannot open file for writing heat map");
    }

    std::vector<std::vector<Real>> heatmap_to_print(2 * mu_vals_vec.size(),
            std::vector<Real>(2 * lambda_vals_vec.size(), 0.0));

    for(auto axis_type : {AxisType::x_type, AxisType::y_type}) {
        bool is_x_type = axis_type == AxisType::x_type;
        for(auto angle_type : {AngleType::flat, AngleType::steep}) {
            bool is_flat = angle_type == AngleType::flat;

            int mu_idx_begin, mu_idx_end;

            if (is_x_type) {
                mu_idx_begin = mu_vals_vec.size() - 1;
                mu_idx_end = -1;
            } else {
                mu_idx_begin = 0;
                mu_idx_end = mu_vals_vec.size();
            }

            int lambda_idx_begin, lambda_idx_end;

            if (is_flat) {
                lambda_idx_begin = 0;
                lambda_idx_end = lambda_vals_vec.size();
            } else {
                lambda_idx_begin = lambda_vals_vec.size() - 1;
                lambda_idx_end = -1;
            }

            int mu_idx_final = is_x_type ? 0 : mu_vals_vec.size();

            for(int mu_idx = mu_idx_begin; mu_idx != mu_idx_end; (mu_idx_begin < mu_idx_end) ? mu_idx++ : mu_idx--) {
                Real mu = mu_vals_vec.at(mu_idx);

                if (mu == 0.0 and axis_type == AxisType::x_type)
                    continue;

                int lambda_idx_final = is_flat ? 0 : lambda_vals_vec.size();

                for(int lambda_idx = lambda_idx_begin;
                    lambda_idx != lambda_idx_end;
                    (lambda_idx_begin < lambda_idx_end) ? lambda_idx++ : lambda_idx--) {

                    Real lambda = lambda_vals_vec.at(lambda_idx);

                    if (lambda == 0.0 and angle_type == AngleType::flat)
                        continue;

                    DualPoint dp(axis_type, angle_type, lambda, mu);
                    Real dist_value = hms.at(max_level).at(dp);

                    heatmap_to_print.at(mu_idx_final).at(lambda_idx_final) = dist_value;

//                    fmt::print("HM, dp = {}, mu_idx_final = {}, lambda_idx_final = {}, value = {}\n", dp, mu_idx_final,
//                            lambda_idx_final, dist_value);

                    lambda_idx_final++;
                }
                mu_idx_final++;
            }
        }
    }

    for(size_t m_idx = 0; m_idx < heatmap_to_print.size(); ++m_idx) {
        for(size_t l_idx = 0; l_idx < heatmap_to_print.at(m_idx).size(); ++l_idx) {
            ofs << heatmap_to_print.at(m_idx).at(l_idx) << " ";
        }
        ofs << std::endl;
    }

    ofs.close();
    spd::debug("Exit print_heat_map");
#endif
}

int main(int argc, char** argv)
{
    spdlog::set_level(spdlog::level::info);
    //spdlog::set_pattern("[%L] %v");

    using opts::Option;
    using opts::PosOption;
    opts::Options ops;

    bool help = false;
    bool heatmap_only = false;
    bool no_stop_asap = false;
    CalculationParams params;

    std::string bounds_list_str = "local_combined";
    std::string traverse_list_str = "BFS";

    ops >> Option('m', "max-iterations", params.max_depth, "maximal number of iterations (refinements)")
        >> Option('e', "max-error", params.delta, "error threshold")
        >> Option('d', "dim", params.dim, "dim")
        >> Option('i', "initial-depth", params.initialization_depth, "initialization depth")
        >> Option("no-stop-asap", no_stop_asap,
                "don't stop looping over points, if cell cannot be pruned (asap is on by default)")
        >> Option("bounds", bounds_list_str, "bounds to use, separated by ,")
        >> Option("traverse", traverse_list_str, "traverse to use, separated by ,")
#ifdef PRINT_HEAT_MAP
        >> Option("heatmap-only", heatmap_only, "only save heatmap (bruteforce)")
#endif
        >> Option('h', "help", help, "show help message");

    std::string fname_a;
    std::string fname_b;

    if (!ops.parse(argc, argv) || help || !(ops >> PosOption(fname_a) >> PosOption(fname_b))) {
        std::cerr << "Usage: " << argv[0] << "bifiltration-file-1 bifiltration-file-2\n" << ops << std::endl;
        return 1;
    }

    params.stop_asap = not no_stop_asap;

    auto bounds_list = split_by_delim(bounds_list_str, ',');
    auto traverse_list = split_by_delim(traverse_list_str, ',');

    Bifiltration bif_a(fname_a, BifiltrationFormat::rene);
    Bifiltration bif_b(fname_b, BifiltrationFormat::rene);

    bif_a.sanity_check();
    bif_b.sanity_check();

    spd::info("Read bifiltrations {} {}", fname_a, fname_b);

    std::vector<BoundStrategy> bound_strategies;
    std::vector<TraverseStrategy> traverse_strategies;

    for(std::string s : bounds_list) {
        bound_strategies.push_back(bs_from_string(s));
    }

    for(std::string s : traverse_list) {
        traverse_strategies.push_back(ts_from_string(s));
    }

    for(auto bs : bound_strategies) {
        for(auto ts : traverse_strategies) {
            spd::info("Will test combination {} {}", bs, ts);
        }
    }

#ifdef EXPERIMENTAL_TIMING
    struct ExperimentResult {
        CalculationParams params {CalculationParams()};
        int n_hera_calls {0};
        double total_milliseconds_elapsed {0};
        double distance {0};
        double actual_error {std::numeric_limits<double>::max()};
        int actual_max_depth {0};

        int x_wins {0};
        int y_wins {0};
        int ad_wins {0};

        int seconds_elapsed() const
        {
            return static_cast<int>(total_milliseconds_elapsed / 1000);
        }

        double savings_ratio_old() const
        {
            long int max_possible_calls = 0;
            long int calls_on_level = 4;
            for(int i = 0; i <= actual_max_depth; ++i) {
                max_possible_calls += calls_on_level;
                calls_on_level *= 4;
            }
            return static_cast<double>(n_hera_calls) / static_cast<double>(max_possible_calls);
        }

        double savings_ratio() const
        {
            return static_cast<double>(n_hera_calls) / calls_on_actual_max_depth();
        }

        long long int calls_on_actual_max_depth() const
        {
            long long int result = 1;
            for(int i = 0; i < actual_max_depth; ++i) {
                result *= 4;
            }
            return result;
        }

        ExperimentResult() { }

        ExperimentResult(CalculationParams p, int nhc, double tme, double d)
                :
                params(p), n_hera_calls(nhc), total_milliseconds_elapsed(tme), distance(d) { }
    };

    const int n_repetitions = 1;

    if (heatmap_only) {
        bound_strategies.clear();
        bound_strategies.push_back(BoundStrategy::bruteforce);
        traverse_strategies.clear();
        traverse_strategies.push_back(TraverseStrategy::breadth_first);
    }

    std::map<std::tuple<BoundStrategy, TraverseStrategy>, ExperimentResult> results;
    for(BoundStrategy bound_strategy : bound_strategies) {
        for(TraverseStrategy traverse_strategy : traverse_strategies) {
            CalculationParams params_experiment;
            params_experiment.bound_strategy = bound_strategy;
            params_experiment.traverse_strategy = traverse_strategy;
            params_experiment.max_depth = params.max_depth;
            params_experiment.initialization_depth = params.initialization_depth;
            params_experiment.delta = params.delta;
            params_experiment.dim = params.dim;
            params_experiment.hera_epsilon = params.hera_epsilon;
            params_experiment.stop_asap = params.stop_asap;

            if (traverse_strategy == TraverseStrategy::depth_first and bound_strategy == BoundStrategy::bruteforce)
                continue;

            // if bruteforce, clamp max iterations number to 7,
            // save user-provided max_iters in user_max_iters and restore it later.
            // remember: params is passed by reference to return real relative error and heat map

            int user_max_iters = params.max_depth;
            if (bound_strategy == BoundStrategy::bruteforce and not heatmap_only) {
                params_experiment.max_depth = std::min(7, params.max_depth);
            }
            double total_milliseconds_elapsed = 0;
            int total_n_hera_calls = 0;
            Real dist;
            for(int i = 0; i < n_repetitions; ++i) {
                spd::debug("Processing bound_strategy {}, traverse_strategy {}, iteration = {}", bound_strategy,
                        traverse_strategy, i);
                auto t1 = std::chrono::high_resolution_clock().now();
                dist = matching_distance(bif_a, bif_b, params_experiment);
                auto t2 = std::chrono::high_resolution_clock().now();
                total_milliseconds_elapsed += std::chrono::duration_cast<std::chrono::milliseconds>(
                        t2 - t1).count();
                total_n_hera_calls += params_experiment.n_hera_calls;
            }

            auto key = std::make_tuple(bound_strategy, traverse_strategy);
            results[key].params = params_experiment;
            results[key].n_hera_calls = total_n_hera_calls / n_repetitions;
            results[key].total_milliseconds_elapsed = total_milliseconds_elapsed / n_repetitions;
            results[key].distance = dist;
            results[key].actual_error = params_experiment.actual_error;
            results[key].actual_max_depth = params_experiment.actual_max_depth;

            spd::info(
                    "Done (bound = {}, traverse = {}), n_hera_calls = {}, time = {} sec, d = {}, error = {}, savings = {}, max_depth = {}",
                    bound_strategy, traverse_strategy, results[key].n_hera_calls, results[key].seconds_elapsed(),
                    dist,
                    params_experiment.actual_error, results[key].savings_ratio(), results[key].actual_max_depth);

            if (bound_strategy == BoundStrategy::bruteforce) { params_experiment.max_depth = user_max_iters; }

#ifdef PRINT_HEAT_MAP
            if (bound_strategy == BoundStrategy::bruteforce) {
                fs::path fname_a_path {fname_a.c_str()};
                fs::path fname_b_path {fname_b.c_str()};
                fs::path fname_a_wo = fname_a_path.filename();
                fs::path fname_b_wo = fname_b_path.filename();
                std::string heat_map_fname = fmt::format("{0}_{1}_dim_{2}_weighted_values_xyp.txt", fname_a_wo.string(),
                        fname_b_wo.string(), params_experiment.dim);
                fs::path path_hm = fname_a_path.replace_filename(fs::path(heat_map_fname.c_str()));
                spd::debug("Saving heatmap to {}", heat_map_fname);
                print_heat_map(params_experiment.heat_maps, path_hm.string(), params);
            }
#endif
            spd::debug("Finished processing bound_strategy {}", bound_strategy);
        }
    }

//    std::cout << "File_1;File_2;Boundstrategy;TraverseStrategy;InitalDepth;NHeraCalls;SavingsRatio;Time;Distance;Error;PushStrategy;MaxDepth;CallsOnMaxDepth;Delta;Dimension" << std::endl;
    for(auto bs : bound_strategies) {
        for(auto ts : traverse_strategies) {
            auto key = std::make_tuple(bs, ts);

            fs::path fname_a_path {fname_a.c_str()};
            fs::path fname_b_path {fname_b.c_str()};
            fs::path fname_a_wo = fname_a_path.filename();
            fs::path fname_b_wo = fname_b_path.filename();

            std::cout << fname_a_wo.string() << ";" << fname_b_wo.string() << ";" << bs << ";" << ts << ";";
            std::cout << results[key].params.initialization_depth << ";";
            std::cout << results[key].n_hera_calls << ";"
                      << results[key].savings_ratio() << ";"
                      << results[key].total_milliseconds_elapsed << ";"
                      << results[key].distance << ";"
                      << results[key].actual_error << ";"
                      << "xyp" << ";"
                      << results[key].actual_max_depth << ";"
                      << results[key].calls_on_actual_max_depth() << ";"
                      << params.delta << ";"
                      << params.dim
                      << std::endl;
        }
    }
#else
    params.bound_strategy = bound_strategies.back();
    params.traverse_strategy = traverse_strategies.back();

    Real dist = matching_distance(bif_a, bif_b, params);
    std::cout << dist << std::endl;
#endif
    return 0;
}
