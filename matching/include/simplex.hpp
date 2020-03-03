namespace md {

    template<class Real>
    Simplex<Real>::Simplex(Index id, Point<Real> birth, int dim, const Column& bdry)
            :
            id_(id),
            pos_(birth),
            dim_(dim),
            facet_indices_(bdry) { }

    template<class Real>
    void Simplex<Real>::translate(Real a)
    {
        pos_.translate(a);
    }

    template<class Real>
    void Simplex<Real>::init_rivet(std::string s)
    {
        auto delim_pos = s.find_first_of(";");
        assert(delim_pos > 0);
        std::string vertices_str = s.substr(0, delim_pos);
        std::string pos_str = s.substr(delim_pos + 1);
        assert(not vertices_str.empty() and not pos_str.empty());
        // get vertices
        std::stringstream vertices_ss(vertices_str);
        int dim = 0;
        int vertex;
        while (vertices_ss >> vertex) {
            dim++;
            vertices_.push_back(vertex);
        }
        //
        std::sort(vertices_.begin(), vertices_.end());
        assert(dim > 0);

        std::stringstream pos_ss(pos_str);
        // TODO: get rid of 1-criticaltiy assumption
        pos_ss >> pos_.x >> pos_.y;
    }

    template<class Real>
    void Simplex<Real>::init_phat_like(std::string s)
    {
        facet_indices_.clear();
        std::stringstream ss(s);
        ss >> dim_ >> pos_.x >> pos_.y;
        if (dim_ > 0) {
            facet_indices_.reserve(dim_ + 1);
            for (int j = 0; j <= dim_; j++) {
                Index k;
                ss >> k;
                facet_indices_.push_back(k);
            }
        }
    }

    template<class Real>
    Simplex<Real>::Simplex(Index _id, std::string s, BifiltrationFormat input_format)
            :id_(_id)
    {
        switch (input_format) {
            case BifiltrationFormat::phat_like :
                init_phat_like(s);
                break;
            case BifiltrationFormat::rivet :
                init_rivet(s);
                break;
        }
    }

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const Simplex<Real>& x)
    {
        os << "Simplex<Real>(id = " << x.id() << ", dim = " << x.dim();
        os << ", boundary = " << container_to_string(x.boundary()) << ", pos = " << x.position() << ")";
        return os;
    }
}
