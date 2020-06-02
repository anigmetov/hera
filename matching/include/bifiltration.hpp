namespace md {

    template<class Real>
    void Bifiltration<Real>::init()
    {
        auto lower_left = max_point<Real>();
        auto upper_right = min_point<Real>();
        for(const auto& simplex : simplices_) {
            lower_left = greatest_lower_bound<>(lower_left, simplex.position());
            upper_right = least_upper_bound<>(upper_right, simplex.position());
            maximal_dim_ = std::max(maximal_dim_, simplex.dim());
        }
        bounding_box_ = Box<Real>(lower_left, upper_right);
    }

    template<class Real>
    Bifiltration<Real>::Bifiltration(const std::string& fname)
    {
        std::ifstream ifstr {fname.c_str()};
        if (!ifstr.good()) {
            std::string error_message = "Cannot open file " + fname;
            std::cerr << error_message << std::endl;
            throw std::runtime_error(error_message);
        }

        BifiltrationFormat input_format;

        std::string s;

        while(ignore_line(s)) {
            std::getline(ifstr, s);
        }

        if (s == "bifiltration") {
            input_format = BifiltrationFormat::rivet;
        } else if (s == "bifiltration_phat_like") {
            input_format = BifiltrationFormat::phat_like;
        } else {
            std::cerr << "Unknown format: '" << s << "' in file " << fname << std::endl;
            throw std::runtime_error("unknown bifiltration format");
        }

        switch(input_format) {
            case BifiltrationFormat::rivet :
                rivet_format_reader(ifstr);
                break;
            case BifiltrationFormat::phat_like :
                phat_like_format_reader(ifstr);
                break;
        }

        ifstr.close();

        init();
    }

    template<class Real>
    void Bifiltration<Real>::rivet_format_reader(std::ifstream& ifstr)
    {
        std::string s;
        // read axes names, ignore them
        std::getline(ifstr, s);
        std::getline(ifstr, s);

        Index index = 0;
        while(std::getline(ifstr, s)) {
            if (!ignore_line(s)) {
                simplices_.emplace_back(index++, s, BifiltrationFormat::rivet);
            }
        }
    }

    template<class Real>
    void Bifiltration<Real>::phat_like_format_reader(std::ifstream& ifstr)
    {
        // read stream line by line; do not use >> operator
        std::string s;
        std::getline(ifstr, s);

        // first line contains number of simplices
        long n_simplices = std::stol(s);

        // all other lines represent a simplex
        Index index = 0;
        while(index < n_simplices) {
            std::getline(ifstr, s);
            if (!ignore_line(s)) {
                simplices_.emplace_back(index++, s, BifiltrationFormat::phat_like);
            }
        }
    }

    template<class Real>
    void Bifiltration<Real>::scale(Real lambda)
    {
        for(auto& s : simplices_) {
            s.scale(lambda);
        }
        init();
    }

    template<class Real>
    void Bifiltration<Real>::sanity_check() const
    {
#ifdef DEBUG
        // check that boundary has correct number of simplices,
        // each bounding simplex has correct dim
        // and appears in the filtration before the simplex it bounds
        for(const auto& s : simplices_) {
            assert(s.dim() >= 0);
            assert(s.dim() == 0 or s.dim() + 1 == (int) s.boundary().size());
            for(auto bdry_idx : s.boundary()) {
                Simplex bdry_simplex = simplices()[bdry_idx];
                assert(bdry_simplex.dim() == s.dim() - 1);
                assert(bdry_simplex.position().is_less(s.position(), false));
            }
        }
#endif
    }

    template<class Real>
    Diagram<Real> Bifiltration<Real>::weighted_slice_diagram(const DualPoint<Real>& line, int dim) const
    {
        DiagramKeeper<Real> dgm;

        // make a copy for now; I want slice_diagram to be const
        std::vector<Simplex<Real>> simplices(simplices_);

//        std::vector<Simplex> simplices;
//        simplices.reserve(simplices_.size() / 2);
//        for(const auto& s : simplices_) {
//            if (s.dim() <= dim + 1 and s.dim() >= dim)
//                simplices.emplace_back(s);
//        }

        for(auto& simplex : simplices) {
            Real value = line.weighted_push(simplex.position());
            simplex.set_value(value);
        }

        std::sort(simplices.begin(), simplices.end(),
                [](const Simplex<Real>& a, const Simplex<Real>& b) { return a.value() < b.value(); });
        std::map<Index, Index> index_map;
        for(Index i = 0; i < (int) simplices.size(); i++) {
            index_map[simplices[i].id()] = i;
        }

        phat::boundary_matrix<> phat_matrix;
        phat_matrix.set_num_cols(simplices.size());
        std::vector<Index> bd_in_slice_filtration;
        for(Index i = 0; i < (int) simplices.size(); i++) {
            phat_matrix.set_dim(i, simplices[i].dim());
            bd_in_slice_filtration.clear();
            //std::cout << "new col" << i << std::endl;
            for(int j = 0; j < (int) simplices[i].boundary().size(); j++) {
                // F[i] contains the indices of its facet wrt to the
                // original filtration. We have to express it, however,
                // wrt to the filtration along the slice. That is why
                // we need the index_map
                //std::cout << "Found " << F[i].bd[j] << ", returning " << index_map[F[i].bd[j]] << std::endl;
                bd_in_slice_filtration.push_back(index_map[simplices[i].boundary()[j]]);
            }
            std::sort(bd_in_slice_filtration.begin(), bd_in_slice_filtration.end());
            phat_matrix.set_col(i, bd_in_slice_filtration);
        }
        phat::persistence_pairs phat_persistence_pairs;
        phat::compute_persistence_pairs<phat::twist_reduction>(phat_persistence_pairs, phat_matrix);

        dgm.clear();
        constexpr Real real_inf = std::numeric_limits<Real>::infinity();
        for(long i = 0; i < (long) phat_persistence_pairs.get_num_pairs(); i++) {
            std::pair<phat::index, phat::index> new_pair = phat_persistence_pairs.get_pair(i);
            bool is_finite_pair = new_pair.second != phat::k_infinity_index;
            Real birth = simplices.at(new_pair.first).value();
            Real death = is_finite_pair ? simplices.at(new_pair.second).value() : real_inf;
            int dim = simplices[new_pair.first].dim();
            assert(dim + 1 == simplices[new_pair.second].dim());
            if (birth != death) {
                dgm.add_point(dim, birth, death);
            }
        }

        return dgm.get_diagram(dim);
    }

    template<class Real>
    Box<Real> Bifiltration<Real>::bounding_box() const
    {
        return bounding_box_;
    }

    template<class Real>
    Real Bifiltration<Real>::minimal_coordinate() const
    {
        return std::min(bounding_box_.lower_left().x, bounding_box_.lower_left().y);
    }

    template<class Real>
    void Bifiltration<Real>::translate(Real a)
    {
        bounding_box_.translate(a);
        for(auto& simplex : simplices_) {
            simplex.translate(a);
        }
    }

    template<class Real>
    Real Bifiltration<Real>::max_x() const
    {
        if (simplices_.empty())
            return 1;
        auto me = std::max_element(simplices_.cbegin(), simplices_.cend(),
                [](const auto& s_a, const auto& s_b) { return s_a.position().x < s_b.position().x; });
        assert(me != simplices_.cend());
        return me->position().x;
    }

    template<class Real>
    Real Bifiltration<Real>::max_y() const
    {
        if (simplices_.empty())
            return 1;
        auto me = std::max_element(simplices_.cbegin(), simplices_.cend(),
                [](const auto& s_a, const auto& s_b) { return s_a.position().y < s_b.position().y; });
        assert(me != simplices_.cend());
        return me->position().y;
    }

    template<class Real>
    Real Bifiltration<Real>::min_x() const
    {
        if (simplices_.empty())
            return 0;
        auto me = std::min_element(simplices_.cbegin(), simplices_.cend(),
                [](const auto& s_a, const auto& s_b) { return s_a.position().x < s_b.position().x; });
        assert(me != simplices_.cend());
        return me->position().x;
    }

    template<class Real>
    Real Bifiltration<Real>::min_y() const
    {
        if (simplices_.empty())
            return 0;
        auto me = std::min_element(simplices_.cbegin(), simplices_.cend(),
                [](const auto& s_a, const auto& s_b) { return s_a.position().y < s_b.position().y; });
        assert(me != simplices_.cend());
        return me->position().y;
    }

    template<class Real>
    void Bifiltration<Real>::add_simplex(Index _id, Point<Real> birth, int _dim, const Column& _bdry)
    {
        simplices_.emplace_back(_id, birth, _dim, _bdry);
    }

    template<class Real>
    void Bifiltration<Real>::save(const std::string& filename, md::BifiltrationFormat format)
    {
        switch(format) {
            case BifiltrationFormat::rivet:
                throw std::runtime_error("Not implemented");
                break;
            case BifiltrationFormat::phat_like: {
                std::ofstream f(filename);
                if (not f.good()) {
                    std::cerr << "Bifiltration::save: cannot open file " << filename << std::endl;
                    throw std::runtime_error("Cannot open file for writing ");
                }
                f << simplices_.size() << "\n";

                for(const auto& s : simplices_) {
                    f << s.dim() << " " << s.position().x << " " << s.position().y << " ";
                    for(int b : s.boundary()) {
                        f << b << " ";
                    }
                    f << std::endl;
                }

            }
                break;
        }
    }

    template<class Real>
    void Bifiltration<Real>::postprocess_rivet_format()
    {
        std::map<Column, Index> facets_to_ids;

        // fill the map
        for(Index i = 0; i < (Index) simplices_.size(); ++i) {
            assert(simplices_[i].id() == i);
            facets_to_ids[simplices_[i].vertices_] = i;
        }

//        for(const auto& s : simplices_) {
//            facets_to_ids[s] = s.id();
//        }

        // main loop
        for(auto& s : simplices_) {
            assert(not s.vertices_.empty());
            assert(s.facet_indices_.empty());
            Column facet_indices;
            for(Index i = 0; i <= s.dim(); ++i) {
                Column facet;
                for(Index j : s.vertices_) {
                    if (j != i)
                        facet.push_back(j);
                }
                auto facet_index = facets_to_ids.at(facet);
                facet_indices.push_back(facet_index);
            }
            s.facet_indices_ = facet_indices;
        } // loop over simplices
    }

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const Bifiltration<Real>& bif)
    {
        os << "Bifiltration [" << std::endl;
        for(const auto& s : bif.simplices()) {
            os << s << std::endl;
        }
        os << "]" << std::endl;
        return os;
    }

    template<class Real>
    BifiltrationProxy<Real>::BifiltrationProxy(const Bifiltration<Real>& bif, int dim)
            :
            dim_(dim),
            bif_(bif)
    {
        cache_positions();
    }

    template<class Real>
    void BifiltrationProxy<Real>::cache_positions() const
    {
        cached_positions_.clear();
        for(const auto& simplex : bif_.simplices()) {
            if (simplex.dim() == dim_ or simplex.dim() == dim_ + 1)
                cached_positions_.push_back(simplex.position());
        }
    }

    template<class Real>
    PointVec<Real>
    BifiltrationProxy<Real>::positions() const
    {
        if (cached_positions_.empty()) {
            cache_positions();
        }
        return cached_positions_;
    }

    // translate all points by vector (a,a)
    template<class Real>
    void BifiltrationProxy<Real>::translate(Real a)
    {
        bif_.translate(a);
    }

    // return minimal value of x- and y-coordinates
    // among all simplices
    template<class Real>
    Real BifiltrationProxy<Real>::minimal_coordinate() const
    {
        return bif_.minimal_coordinate();
    }

    // return box that contains positions of all simplices
    template<class Real>
    Box<Real> BifiltrationProxy<Real>::bounding_box() const
    {
        return bif_.bounding_box();
    }

    template<class Real>
    Real BifiltrationProxy<Real>::max_x() const
    {
        return bif_.max_x();
    }

    template<class Real>
    Real BifiltrationProxy<Real>::max_y() const
    {
        return bif_.max_y();
    }

    template<class Real>
    Real BifiltrationProxy<Real>::min_x() const
    {
        return bif_.min_x();
    }

    template<class Real>
    Real BifiltrationProxy<Real>::min_y() const
    {
        return bif_.min_y();
    }


    template<class Real>
    Diagram<Real> BifiltrationProxy<Real>::weighted_slice_diagram(const DualPoint<Real>& slice) const
    {
        return bif_.weighted_slice_diagram(slice, dim_);
    }

}

