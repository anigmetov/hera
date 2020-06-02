namespace md {

    template<class R, class T>
    void DistanceCalculator<R, T>::check_upper_bound(const CellWithValue<R>& dual_cell) const
    {
        const int n_samples_lambda = 100;
        const int n_samples_mu = 100;
        DualBox<R> db = dual_cell.dual_box();
        R min_lambda = db.lambda_min();
        R max_lambda = db.lambda_max();
        R min_mu = db.mu_min();
        R max_mu = db.mu_max();

        R h_lambda = (max_lambda - min_lambda) / n_samples_lambda;
        R h_mu = (max_mu - min_mu) / n_samples_mu;
        for(int i = 1; i < n_samples_lambda; ++i) {
            for(int j = 1; j < n_samples_mu; ++j) {
                R lambda = min_lambda + i * h_lambda;
                R mu = min_mu + j * h_mu;
                DualPoint<R> l(db.axis_type(), db.angle_type(), lambda, mu);
                R other_result = distance_on_line_const(l);
                R diff = fabs(dual_cell.stored_upper_bound() - other_result);
                if (other_result > dual_cell.stored_upper_bound()) {
                    throw std::runtime_error("Wrong delta estimate");
                }
            }
        }
    }

    // for all lines l, l' inside dual box,
    // find the upper bound on the difference of weighted  pushes of p
    template<class R, class T>
    R
    DistanceCalculator<R, T>::get_max_displacement_single_point(const CellWithValue<R>& dual_cell, ValuePoint vp,
            const Point<R>& p) const
    {
        assert(p.x >= 0 && p.y >= 0);

#ifdef MD_DEBUG
        std::vector<long long int> debug_ids = {3, 13, 54, 218, 350, 382, 484, 795, 2040, 8415, 44076};
        bool debug = false; // std::find(debug_ids.begin(), debug_ids.end(), dual_cell.id) != debug_ids.end();
#endif
        DualPoint<R> line = dual_cell.value_point(vp);
        const R base_value = line.weighted_push(p);

        R result = 0.0;
        for(DualPoint<R> dp : dual_cell.dual_box().critical_points(p)) {
            R dp_value = dp.weighted_push(p);
            result = std::max(result, fabs(base_value - dp_value));
        }

#ifdef MD_DO_FULL_CHECK
        auto db = dual_cell.dual_box();
        std::uniform_real_distribution<R> dlambda(db.lambda_min(), db.lambda_max());
        std::uniform_real_distribution<R> dmu(db.mu_min(), db.mu_max());
        std::mt19937 gen(1);
        for(int i = 0; i < 1000; ++i) {
            R lambda = dlambda(gen);
            R mu = dmu(gen);
            DualPoint<R> dp_random { db.axis_type(), db.angle_type(), lambda, mu };
            R dp_value = dp_random.weighted_push(p);
            if (fabs(base_value - dp_value) > result) {
                throw std::runtime_error("error in get_max_displacement_single_value");
            }
        }
#endif

        return result;
    }

    template<class R, class T>
    typename DistanceCalculator<R, T>::CellValueVector DistanceCalculator<R, T>::get_initial_dual_grid(R& lower_bound)
    {
        CellValueVector result = get_refined_grid(params_.initialization_depth, false, true);

        lower_bound = -1;
        for(const auto& dc : result) {
            lower_bound = std::max(lower_bound, dc.max_corner_value());
        }

        assert(lower_bound >= 0);

        for(auto& dual_cell : result) {
            R good_enough_ub = get_good_enough_upper_bound(lower_bound);
            R max_value_on_cell = get_upper_bound(dual_cell, good_enough_ub);
            dual_cell.set_max_possible_value(max_value_on_cell);

#ifdef MD_DO_FULL_CHECK
            check_upper_bound(dual_cell);
#endif
        }

        return result;
    }

    template<class R, class T>
    typename DistanceCalculator<R, T>::CellValueVector
    DistanceCalculator<R, T>::get_refined_grid(int init_depth, bool calculate_on_intermediate, bool calculate_on_last)
    {
        const R y_max = std::max(module_a_.max_y(), module_b_.max_y());
        const R x_max = std::max(module_a_.max_x(), module_b_.max_x());

        const R lambda_min = 0;
        const R lambda_max = 1;

        const R mu_min = 0;

        DualBox<R> x_flat(DualPoint<R>(AxisType::x_type, AngleType::flat, lambda_min, mu_min),
                DualPoint<R>(AxisType::x_type, AngleType::flat, lambda_max, x_max));

        DualBox<R>  x_steep(DualPoint<R>(AxisType::x_type, AngleType::steep, lambda_min, mu_min),
                DualPoint<R>(AxisType::x_type, AngleType::steep, lambda_max, x_max));

        DualBox<R> y_flat(DualPoint<R>(AxisType::y_type, AngleType::flat, lambda_min, mu_min),
                DualPoint<R>(AxisType::y_type, AngleType::flat, lambda_max, y_max));

        DualBox<R> y_steep(DualPoint<R>(AxisType::y_type, AngleType::steep, lambda_min, mu_min),
                DualPoint<R>(AxisType::y_type, AngleType::steep, lambda_max, y_max));

        CellWithValue<R> x_flat_cell(x_flat, 0);
        CellWithValue<R> x_steep_cell(x_steep, 0);
        CellWithValue<R> y_flat_cell(y_flat, 0);
        CellWithValue<R> y_steep_cell(y_steep, 0);

        if (init_depth == 0) {
            DualPoint<R> diagonal_x_flat(AxisType::x_type, AngleType::flat, 1, 0);

            R diagonal_value = distance_on_line(diagonal_x_flat);
            n_hera_calls_per_level_[0]++;

            x_flat_cell.set_value_at(ValuePoint::lower_right, diagonal_value);
            y_flat_cell.set_value_at(ValuePoint::lower_right, diagonal_value);
            x_steep_cell.set_value_at(ValuePoint::lower_right, diagonal_value);
            y_steep_cell.set_value_at(ValuePoint::lower_right, diagonal_value);
        }

#ifdef MD_DEBUG
        x_flat_cell.id = 1;
        x_steep_cell.id = 2;
        y_flat_cell.id = 3;
        y_steep_cell.id = 4;
        CellWithValue<R>::max_id = 4;
#endif

        CellValueVector result {x_flat_cell, x_steep_cell, y_flat_cell, y_steep_cell};

        if (init_depth == 0) {
            return result;
        }

        CellValueVector refined_result;

        for(int i = 1; i <= init_depth; ++i) {
            refined_result.clear();
            for(const auto& dual_cell : result) {
                for(auto refined_cell : dual_cell.get_refined_cells()) {
                    // we calculate for init_dept - 1, not init_depth,
                    // because we want the cells to have value at a corner
                    if ((i == init_depth - 1 and calculate_on_last) or calculate_on_intermediate)
                        set_cell_central_value(refined_cell);
                    refined_result.push_back(refined_cell);
                }
            }
            result = std::move(refined_result);
        }
        return result;
    }

    template<class R, class T>
    DistanceCalculator<R, T>::DistanceCalculator(const T& a,
            const T& b,
            CalculationParams<R>& params)
            :
            module_a_(a),
            module_b_(b),
            params_(params)
    {
        // make all coordinates non-negative
        auto min_coord = std::min(module_a_.minimal_coordinate(),
                module_b_.minimal_coordinate());
        if (min_coord < 0) {
            module_a_.translate(-min_coord);
            module_b_.translate(-min_coord);
        }

        assert(std::min({module_a_.min_x(), module_b_.min_x(), module_a_.min_y(),
                         module_b_.min_y()}) >= 0);

    }

    template<class R, class T>
    R DistanceCalculator<R, T>::get_max_x(int module) const
    {
        return (module == 0) ? module_a_.max_x() : module_b_.max_x();
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::get_max_y(int module) const
    {
        return (module == 0) ? module_a_.max_y() : module_b_.max_y();
    }

    template<class R, class T>
    R
    DistanceCalculator<R, T>::get_local_refined_bound(const DualBox<R>& dual_box) const
    {
        return get_local_refined_bound(0, dual_box) + get_local_refined_bound(1, dual_box);
    }

    template<class R, class T>
    R
    DistanceCalculator<R, T>::get_local_refined_bound(int module, const DualBox<R>& dual_box) const
    {
        R d_lambda = dual_box.lambda_max() - dual_box.lambda_min();
        R d_mu = dual_box.mu_max() - dual_box.mu_min();
        R result;
        if (dual_box.axis_type() == AxisType::x_type) {
            if (dual_box.is_flat()) {
                result = dual_box.lambda_max() * d_mu + (get_max_x(module) - dual_box.mu_min()) * d_lambda;
            } else {
                result = d_mu + get_max_y(module) * d_lambda;
            }
        } else {
            // y-type
            if (dual_box.is_flat()) {
                result = d_mu + get_max_x(module) * d_lambda;
            } else {
                // steep
                result = dual_box.lambda_max() * d_mu + (get_max_y(module) - dual_box.mu_min()) * d_lambda;
            }
        }
        return result;
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::get_local_dual_bound(int module, const DualBox<R>& dual_box) const
    {
        R dlambda = dual_box.lambda_max() - dual_box.lambda_min();
        R dmu = dual_box.mu_max() - dual_box.mu_min();

        if (dual_box.is_flat()) {
            return get_max_x(module) * dlambda + dmu;
        } else {
            return get_max_y(module) * dlambda + dmu;
        }
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::get_local_dual_bound(const DualBox<R>& dual_box) const
    {
        return get_local_dual_bound(0, dual_box) + get_local_dual_bound(1, dual_box);
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::get_upper_bound(const CellWithValue<R>& dual_cell, R good_enough_ub) const
    {
        assert(good_enough_ub >= 0);

        switch(params_.bound_strategy) {
            case BoundStrategy::bruteforce:
                return std::numeric_limits<R>::max();

            case BoundStrategy::local_dual_bound:
                return dual_cell.min_value() + get_local_dual_bound(dual_cell.dual_box());

            case BoundStrategy::local_dual_bound_refined:
                return dual_cell.min_value() + get_local_refined_bound(dual_cell.dual_box());

            case BoundStrategy::local_combined: {
                R cheap_upper_bound = dual_cell.min_value() + get_local_refined_bound(dual_cell.dual_box());
                if (cheap_upper_bound < good_enough_ub) {
                    return cheap_upper_bound;
                } else {
                    [[fallthrough]];
                }
            }

            case BoundStrategy::local_dual_bound_for_each_point: {
                R result = std::numeric_limits<R>::max();
                for(ValuePoint vp : k_corner_vps) {
                    if (not dual_cell.has_value_at(vp)) {
                        continue;
                    }

                    R base_value = dual_cell.value_at(vp);
                    R bound_dgm_a = get_single_dgm_bound(dual_cell, vp, 0, good_enough_ub);

                    if (params_.stop_asap and bound_dgm_a + base_value >= good_enough_ub) {
                        // we want to return a valid upper bound, not just something that will prevent discarding the cell
                        // and we don't want to compute pushes for points in second bifiltration.
                        // so just return a constant time bound
                        return dual_cell.min_value() + get_local_refined_bound(dual_cell.dual_box());
                    }

                    R bound_dgm_b = get_single_dgm_bound(dual_cell, vp, 1,
                            std::max(R(0), good_enough_ub - bound_dgm_a));

                    result = std::min(result, base_value + bound_dgm_a + bound_dgm_b);

                    if (params_.stop_asap and result < good_enough_ub) {
                        break;
                    }
                }
                return result;
            }
        }
        // to suppress compiler warning
        return std::numeric_limits<R>::max();
    }

    // find maximal displacement of weighted points of m for all lines in dual_box
    template<class R, class T>
    R
    DistanceCalculator<R, T>::get_single_dgm_bound(const CellWithValue<R>& dual_cell,
            ValuePoint vp,
            int module,
            R good_enough_value) const
    {
        R result = 0;
        Point<R> max_point;

        const T& m = (module == 0) ? module_a_ : module_b_;
        for(const auto& position : m.positions()) {
            R x = get_max_displacement_single_point(dual_cell, vp, position);

            if (x > result) {
                result = x;
                max_point = position;
            }

            if (params_.stop_asap and result > good_enough_value) {
                // we want to return a valid upper bound,
                // now we just see it is worse than we need, but it may be even more
                // just return a valid upper bound
                result = get_local_refined_bound(dual_cell.dual_box());
                break;
            }
        }

        return result;
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::distance()
    {
        return get_distance_pq();
    }

    // calculate weighted bottleneneck distance between slices on line
    // increments hera calls counter
    template<class R, class T>
    R DistanceCalculator<R, T>::distance_on_line(DualPoint<R> line)
    {
        ++n_hera_calls_;
        R result = distance_on_line_const(line);
        return result;
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::distance_on_line_const(DualPoint<R> line) const
    {
        // TODO: think about this - how to call Hera
        auto dgm_a = module_a_.weighted_slice_diagram(line);
        auto dgm_b = module_b_.weighted_slice_diagram(line);
        R result;
        if (params_.hera_epsilon > static_cast<R>(0)) {
            result = hera::bottleneckDistApprox(dgm_a, dgm_b, params_.hera_epsilon) / ( params_.hera_epsilon + 1);
        } else {
            result = hera::bottleneckDistExact(dgm_a, dgm_b);
        }
        return result;
    }

    template<class R, class T>
    R DistanceCalculator<R, T>::get_good_enough_upper_bound(R lower_bound) const
    {
        R result;
        // in upper_bound strategy we only prune cells if they cannot improve the lower bound,
        // otherwise the experiment is supposed to run indefinitely
        if (params_.traverse_strategy == TraverseStrategy::upper_bound) {
            result = lower_bound;
        } else {
            result = (1.0 + params_.delta) * lower_bound;
        }
        return result;
    }

    // helper function
    // calculate weighted bt distance on cell center,
    // assign distance value to cell, keep it in heat_map, and return
    template<class R, class T>
    void DistanceCalculator<R, T>::set_cell_central_value(CellWithValue<R>& dual_cell)
    {
        DualPoint<R> central_line {dual_cell.center()};

        R new_value = distance_on_line(central_line);
        n_hera_calls_per_level_[dual_cell.level() + 1]++;
        dual_cell.set_value_at(ValuePoint::center, new_value);
        params_.actual_max_depth = std::max(params_.actual_max_depth, dual_cell.level() + 1);

#ifdef PRINT_HEAT_MAP
        if (params_.bound_strategy == BoundStrategy::bruteforce) {
            if (dual_cell.level() > params_.initialization_depth + 1
                    and params_.heat_maps[dual_cell.level()].count(dual_cell.center()) > 0) {
                auto existing = params_.heat_maps[dual_cell.level()].find(dual_cell.center());
            }
            assert(dual_cell.level() <= params_.initialization_depth + 1
                    or params_.heat_maps[dual_cell.level()].count(dual_cell.center()) == 0);
            params_.heat_maps[dual_cell.level()][dual_cell.center()] = new_value;
        }
#endif
    }

    // quick-and-dirty hack to efficiently traverse priority queue with dual cells
    // returns maximal possible value on all cells in queue
    // assumes that the underlying container is vector!
    // cell_ptr: pointer to the first element in queue
    // n_cells: queue size
    template<class R, class T>
    R DistanceCalculator<R, T>::get_max_possible_value(const CellWithValue<R>* cell_ptr, int n_cells)
    {
        R result = (n_cells > 0) ? cell_ptr->stored_upper_bound() : 0;
        for(int i = 0; i < n_cells; ++i, ++cell_ptr) {
            result = std::max(result, cell_ptr->stored_upper_bound());
        }
        return result;
    }

    // helper function:
    // return current error from lower and upper bounds
    // and save it in params_ (hence not const)
    template<class R, class T>
    R DistanceCalculator<R, T>::current_error(R lower_bound, R upper_bound)
    {
        R current_error = (lower_bound > 0.0) ? (upper_bound - lower_bound) / lower_bound
                                                 : std::numeric_limits<R>::max();
        params_.actual_error = current_error;

        return current_error;
    }

    // return matching distance
    // use priority queue to store dual cells
    // comparison function depends on the strategies in params_
    // ressets hera calls counter
    template<class R, class T>
    R DistanceCalculator<R, T>::get_distance_pq()
    {
        std::map<int, long> n_cells_considered;
        std::map<int, long> n_cells_pushed_into_queue;
        long int n_too_deep_cells = 0;
        std::map<int, long> n_cells_discarded;
        std::map<int, long> n_cells_pruned;

        std::chrono::high_resolution_clock timer;
        auto start_time = timer.now();

        n_hera_calls_ = 0;
        n_hera_calls_per_level_.clear();


        // if cell is too deep and is not pushed into queue,
        // we still need to take its max value into account;
        // the max over such cells is stored in max_result_on_too_fine_cells
        R upper_bound_on_deep_cells = -1;

        // user-defined less lambda function
        // to regulate priority queue depending on strategy
        auto dual_cell_less = [this](const CellWithValue<R>& a, const CellWithValue<R>& b) {

            int a_level = a.level();
            int b_level = b.level();
            R a_value = a.max_corner_value();
            R b_value = b.max_corner_value();
            R a_ub = a.stored_upper_bound();
            R b_ub = b.stored_upper_bound();
            if (this->params_.traverse_strategy == TraverseStrategy::upper_bound and
                    (not a.has_max_possible_value() or not b.has_max_possible_value())) {
                throw std::runtime_error("no upper bound on cell");
            }
            DualPoint<R> a_lower_left = a.dual_box().lower_left();
            DualPoint<R> b_lower_left = b.dual_box().lower_left();

            switch(this->params_.traverse_strategy) {
                // in both breadth_first searches we want coarser cells
                // to be processed first. Cells with smaller level must be larger,
                // hence the minus in front of level
                case TraverseStrategy::breadth_first:
                    return std::make_tuple(-a_level, a_lower_left)
                            < std::make_tuple(-b_level, b_lower_left);
                case TraverseStrategy::breadth_first_value:
                    return std::make_tuple(-a_level, a_value, a_lower_left)
                            < std::make_tuple(-b_level, b_value, b_lower_left);
                case TraverseStrategy::depth_first:
                    return std::make_tuple(a_value, a_level, a_lower_left)
                            < std::make_tuple(b_value, b_level, b_lower_left);
                case TraverseStrategy::upper_bound:
                    return std::make_tuple(a_ub, a_level, a_lower_left)
                            < std::make_tuple(b_ub, b_level, b_lower_left);
                default:
                    throw std::runtime_error("Forgotten case");
            }
        };

        std::priority_queue<CellWithValue<R>, CellValueVector, decltype(dual_cell_less)> dual_cells_queue(
                dual_cell_less);

        // weighted bt distance on the center of current cell
        R lower_bound = std::numeric_limits<R>::min();

        // init pq and lower bound
        for(auto& init_cell : get_initial_dual_grid(lower_bound)) {
            dual_cells_queue.push(init_cell);
        }

        R upper_bound = get_max_possible_value(&dual_cells_queue.top(), dual_cells_queue.size());

        std::vector<UbExperimentRecord> ub_experiment_results;

        while(not dual_cells_queue.empty()) {

            CellWithValue<R> dual_cell = dual_cells_queue.top();
            dual_cells_queue.pop();
            assert(dual_cell.has_corner_value()
                    and dual_cell.has_max_possible_value()
                    and dual_cell.max_corner_value() <= upper_bound);

            n_cells_considered[dual_cell.level()]++;

            bool discard_cell = false;

            if (not params_.stop_asap) {
                // if stop_asap is on, it is safer to never discard a cell
                if (params_.bound_strategy == BoundStrategy::bruteforce) {
                    discard_cell = false;
                } else if (params_.traverse_strategy == TraverseStrategy::upper_bound) {
                    discard_cell = (dual_cell.stored_upper_bound() <= lower_bound);
                } else {
                    discard_cell = (dual_cell.stored_upper_bound() <= (1.0 + params_.delta) * lower_bound);
                }
            }

            if (discard_cell) {
                n_cells_discarded[dual_cell.level()]++;
                continue;
            }

            // until now, dual_cell knows its value in one of its corners
            // new_value will be the weighted distance at its center
            set_cell_central_value(dual_cell);
            R new_value = dual_cell.value_at(ValuePoint::center);
            lower_bound = std::max(new_value, lower_bound);

            assert(upper_bound >= lower_bound);

            if (current_error(lower_bound, upper_bound) < params_.delta) {
                break;
            }

            // refine cell and push 4 smaller cells into queue
            for(auto refined_cell : dual_cell.get_refined_cells()) {

                if (refined_cell.num_values() == 0)
                    throw std::runtime_error("no value on cell");

                // if delta is smaller than good_enough_value, it allows to prune cell
                R good_enough_ub = get_good_enough_upper_bound(lower_bound);

                // upper bound of the parent holds for refined_cell
                // and can sometimes be smaller!
                R upper_bound_on_refined_cell = std::min(dual_cell.stored_upper_bound(),
                        get_upper_bound(refined_cell, good_enough_ub));

                refined_cell.set_max_possible_value(upper_bound_on_refined_cell);

#ifdef MD_DO_FULL_CHECK
                check_upper_bound(refined_cell);
#endif

                bool prune_cell = false;

                if (refined_cell.level() <= params_.max_depth) {
                    // cell might be added to queue; if it is not added, its maximal value can be safely ignored
                    if (params_.traverse_strategy == TraverseStrategy::upper_bound) {
                        prune_cell = (refined_cell.stored_upper_bound() <= lower_bound);
                    } else if (params_.bound_strategy != BoundStrategy::bruteforce) {
                        prune_cell = (refined_cell.stored_upper_bound() <= (1.0 + params_.delta) * lower_bound);
                    }
                    if (prune_cell)
                        n_cells_pruned[refined_cell.level()]++;
//                        prune_cell = (max_result_on_refined_cell <= lower_bound);
                } else {
                    // cell is too deep, it won't be added to queue
                    // we must memorize maximal value on this cell, because we won't see it anymore
                    prune_cell = true;
                    if (refined_cell.stored_upper_bound() > (1 + params_.delta) * lower_bound) {
                        n_too_deep_cells++;
                    }
                    upper_bound_on_deep_cells = std::max(upper_bound_on_deep_cells, refined_cell.stored_upper_bound());
                }

                if (not prune_cell) {
                    n_cells_pushed_into_queue[refined_cell.level()]++;
                    dual_cells_queue.push(refined_cell);
                }
            } // end loop over refined cells

            if (dual_cells_queue.empty())
                upper_bound = std::max(upper_bound, upper_bound_on_deep_cells);
            else
                upper_bound = std::max(upper_bound_on_deep_cells,
                        get_max_possible_value(&dual_cells_queue.top(), dual_cells_queue.size()));

            if (params_.traverse_strategy == TraverseStrategy::upper_bound) {
                upper_bound = dual_cells_queue.top().stored_upper_bound();

                if (get_hera_calls_number() < 20 || get_hera_calls_number() % 20 == 0) {
                    auto elapsed = timer.now() - start_time;
                    UbExperimentRecord ub_exp_record;

                    ub_exp_record.error = current_error(lower_bound, upper_bound);
                    ub_exp_record.lower_bound = lower_bound;
                    ub_exp_record.upper_bound = upper_bound;
                    ub_exp_record.cell = dual_cells_queue.top();
                    ub_exp_record.n_hera_calls = n_hera_calls_;
                    ub_exp_record.time = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();

#ifdef MD_DO_CHECKS
                    if (ub_experiment_results.size() > 0) {
                        auto prev = ub_experiment_results.back();
                        if (upper_bound > prev.upper_bound) {
                            throw std::runtime_error("die");
                        }

                        if (lower_bound < prev.lower_bound) {
                            throw std::runtime_error("die");
                        }
                    }
#endif

                    ub_experiment_results.emplace_back(ub_exp_record);

                    std::cerr << "[UB_EXPERIMENT]\t" <<  ub_exp_record << "\n";
                }
            }

            assert(upper_bound >= lower_bound);

            if (current_error(lower_bound, upper_bound) < params_.delta) {
                break;
            }
        }

        params_.actual_error = current_error(lower_bound, upper_bound);

        if (n_too_deep_cells > 0) {
            std::cerr << "Warning: error not guaranteed, there were too deep cells.";
            std::cerr << " Increase max_depth or delta, actual error = " << params_.actual_error << std::endl;
        }
        // otherwise actual_error in params can be larger than delta,
        // but this is OK

        if (params_.print_stats) {
            std::cout << "EXIT STATS, cells considered:\n";
            print_map(n_cells_considered);
            std::cout << "EXIT STATS, cells discarded:\n";
            print_map(n_cells_discarded);
            std::cout << "EXIT STATS, cells pruned:\n";
            print_map(n_cells_pruned);
            std::cout << "EXIT STATS, cells pushed:\n";
            print_map(n_cells_pushed_into_queue);
            std::cout << "EXIT STATS, hera calls:\n";
            print_map(n_hera_calls_per_level_);
            std::cout << "EXIT STATS, too deep cells with high value: " << n_too_deep_cells << "\n";
        }

        return lower_bound;
    }

    template<class R, class T>
    int DistanceCalculator<R, T>::get_hera_calls_number() const
    {
        return n_hera_calls_;
    }

    template<class R>
    R matching_distance(const Bifiltration<R>& bif_a, const Bifiltration<R>& bif_b,
            CalculationParams<R>& params)
    {
        R result;
        // compute distance only in one dimension
        if (params.dim != CalculationParams<R>::ALL_DIMENSIONS) {
            BifiltrationProxy<R> bifp_a(bif_a, params.dim);
            BifiltrationProxy<R> bifp_b(bif_b, params.dim);
            DistanceCalculator<R, BifiltrationProxy<R>> runner(bifp_a, bifp_b, params);
            result = runner.distance();
            params.n_hera_calls = runner.get_hera_calls_number();
        } else {
            // compute distance in all dimensions, return maximal
            result = -1;
            for(int dim = 0; dim < std::max(bif_a.maximal_dim(), bif_b.maximal_dim()); ++dim) {
                BifiltrationProxy<R> bifp_a(bif_a, params.dim);
                BifiltrationProxy<R> bifp_b(bif_a, params.dim);
                DistanceCalculator<R, BifiltrationProxy<R>> runner(bifp_a, bifp_b, params);
                result = std::max(result, runner.distance());
                params.n_hera_calls += runner.get_hera_calls_number();
            }
        }
        return result;
    }


    template<class R>
    R matching_distance(const ModulePresentation<R>& mod_a, const ModulePresentation<R>& mod_b,
            CalculationParams<R>& params)
    {
        DistanceCalculator<R, ModulePresentation<R>> runner(mod_a, mod_b, params);
        R result = runner.distance();
        params.n_hera_calls = runner.get_hera_calls_number();
        return result;
    }
} // namespace md
