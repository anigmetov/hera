namespace md {

    /**
     *
     * @param values vector of length n
     * @return [a_1,...,a_n] such that
     *           1) values[a_1] <= values[a_2] <= ... <= values[a_n]
     *           2) a_1,...,a_n is a permutation of 1,..,n
     */

    template<class T>
    IndexVec get_sorted_indices(const std::vector<T>& values)
    {
        IndexVec result(values.size());
        std::iota(result.begin(), result.end(), 0);
        std::sort(result.begin(), result.end(),
                [&values](size_t a, size_t b) { return values[a] < values[b]; });
        return result;
    }

    // helper function to initialize const member positions_ in ModulePresentation
    template<class Real>
    PointVec<Real> concat_gen_and_rel_positions(const PointVec<Real>& generators,
            const typename ModulePresentation<Real>::RelVec& relations)
    {
        std::unordered_set<Point<Real>> ps(generators.begin(), generators.end());
        for(const auto& rel : relations) {
            ps.insert(rel.position_);
        }
        return PointVec<Real>(ps.begin(), ps.end());
    }


    template<class Real>
    void ModulePresentation<Real>::init_boundaries()
    {
        max_x_ = -std::numeric_limits<Real>::max();
        max_y_ = -std::numeric_limits<Real>::max();
        min_x_ = std::numeric_limits<Real>::max();
        min_y_ = std::numeric_limits<Real>::max();

        for(const auto& gen : positions_) {
            min_x_ = std::min(gen.x, min_x_);
            min_y_ = std::min(gen.y, min_y_);
            max_x_ = std::max(gen.x, max_x_);
            max_y_ = std::max(gen.y, max_y_);
        }

        bounding_box_ = Box<Real>(Point<Real>(min_x_, min_y_), Point<Real>(max_x_, max_y_));
    }


    template<class Real>
    ModulePresentation<Real>::ModulePresentation(const PointVec<Real>& _generators, const RelVec& _relations) :
        generators_(_generators),
        relations_(_relations)
    {
        positions_ = concat_gen_and_rel_positions(generators_, relations_);
        init_boundaries();
    }

    template<class Real>
    void ModulePresentation<Real>::translate(Real a)
    {
        for(auto& g : generators_) {
            g.translate(a);
        }

        for(auto& r : relations_) {
            r.position_.translate(a);
        }

        positions_ = concat_gen_and_rel_positions(generators_, relations_);
        init_boundaries();
    }


    /**
     *
     * @param slice line on which generators are projected
     * @param sorted_indices [a_1,...,a_n] s.t. wpush(generator[a_1]) <= wpush(generator[a_2]) <= ..
     * @param projections sorted weighted pushes of generators
     */

    template<class Real>
    void ModulePresentation<Real>::project_generators(const DualPoint<Real>& slice,
            IndexVec& sorted_indices, RealVec& projections) const
    {
        size_t num_gens = generators_.size();

        RealVec gen_values;
        gen_values.reserve(num_gens);
        for(const auto& pos : generators_) {
            gen_values.push_back(slice.weighted_push(pos));
        }
        sorted_indices = get_sorted_indices(gen_values);
        projections.clear();
        projections.reserve(num_gens);
        for(auto i : sorted_indices) {
            projections.push_back(gen_values[i]);
        }
    }

    template<class Real>
    void ModulePresentation<Real>::project_relations(const DualPoint<Real>& slice, IndexVec& sorted_rel_indices,
            RealVec& projections) const
    {

        size_t num_rels = relations_.size();

        RealVec rel_values;
        rel_values.reserve(num_rels);
        for(const auto& rel : relations_) {
            rel_values.push_back(slice.weighted_push(rel.position_));
        }

        sorted_rel_indices = get_sorted_indices(rel_values);

        projections.clear();
        projections.reserve(num_rels);
        for(auto i : sorted_rel_indices) {
            projections.push_back(rel_values[i]);
        }
    }


    template<class Real>
    void ModulePresentation<Real>::get_slice_projection_matrix(const DualPoint<Real>& slice,
            phat::boundary_matrix<>& phat_matrix,
            RealVec& gen_projections, RealVec& rel_projections) const
    {
        IndexVec sorted_gen_indices, sorted_rel_indices;

        project_generators(slice, sorted_gen_indices, gen_projections);
        project_relations(slice, sorted_rel_indices, rel_projections);

        phat_matrix.set_num_cols(relations_.size());

        for(Index i = 0; i < (Index) relations_.size(); i++) {
            IndexVec current_relation = relations_[sorted_rel_indices[i]].components_;
            for(auto& j : current_relation) {
                j = sorted_gen_indices[j];
            }
            std::sort(current_relation.begin(), current_relation.end());
            // modules do not have dimension, set all to 0
            phat_matrix.set_dim(i, 0);
            phat_matrix.set_col(i, current_relation);
        }
    }


    template<class Real>
    Diagram<Real> ModulePresentation<Real>::weighted_slice_diagram(const DualPoint<Real>& slice) const
    {
        RealVec gen_projections, rel_projections;
        phat::boundary_matrix<> phat_matrix;

        get_slice_projection_matrix(slice, phat_matrix, gen_projections, rel_projections);

        phat::persistence_pairs phat_persistence_pairs;
        phat::compute_persistence_pairs<phat::twist_reduction>(phat_persistence_pairs, phat_matrix);

        Diagram<Real> dgm;

        constexpr Real real_inf = std::numeric_limits<Real>::infinity();

        for(Index i = 0; i < (Index) phat_persistence_pairs.get_num_pairs(); i++) {
            std::pair<phat::index, phat::index> new_pair = phat_persistence_pairs.get_pair(i);
            bool is_finite_pair = new_pair.second != phat::k_infinity_index;
            Real birth = gen_projections.at(new_pair.first);
            Real death = is_finite_pair ? rel_projections.at(new_pair.second) : real_inf;
            if (birth != death) {
                dgm.emplace_back(birth, death);
            }
        }

        return dgm;
    }

    template<class Real>
    PointVec<Real> ModulePresentation<Real>::positions() const
    {
        return positions_;
    }

    template<class Real>
    Box<Real> ModulePresentation<Real>::bounding_box() const
    {
        return bounding_box_;
    }

} // namespace md
