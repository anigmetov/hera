#include "simplex.h"

namespace md {

    std::ostream& operator<<(std::ostream& os, const AbstractSimplex& s)
    {
        os << "AbstractSimplex(id = " << s.id << ", vertices_ = " << container_to_string(s.vertices_) << ")";
        return os;
    }

    bool operator<(const AbstractSimplex& a, const AbstractSimplex& b)
    {
        return a.vertices_ < b.vertices_;
    }

    bool operator==(const AbstractSimplex& s1, const AbstractSimplex& s2)
    {
        return s1.vertices_ == s2.vertices_;
    }

    void AbstractSimplex::push_back(int v)
    {
        vertices_.push_back(v);
        std::sort(vertices_.begin(), vertices_.end());
    }

    AbstractSimplex::AbstractSimplex(std::vector<int> vertices, bool sort)
            :vertices_(vertices)
    {
        if (sort)
            std::sort(vertices_.begin(), vertices_.end());
    }

    std::vector<AbstractSimplex> AbstractSimplex::facets() const
    {
        std::vector<AbstractSimplex> result;
        for (int i = 0; i < static_cast<int>(vertices_.size()); ++i) {
            std::vector<int> facet_vertices;
            facet_vertices.reserve(dim());
            for (int j = 0; j < static_cast<int>(vertices_.size()); ++j) {
                if (j != i)
                    facet_vertices.push_back(vertices_[j]);
            }
            if (!facet_vertices.empty()) {
                result.emplace_back(facet_vertices, false);
            }
        }
        return result;
    }

    Simplex::Simplex(md::Index id, md::Point birth, int dim, const md::Column& bdry)
            :
            id_(id),
            pos_(birth),
            dim_(dim),
            facet_indices_(bdry) { }

    void Simplex::translate(Real a)
    {
        pos_.translate(a);
    }

    void Simplex::init_rivet(std::string s)
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

    void Simplex::init_phat_like(std::string s)
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

    Simplex::Simplex(Index _id, std::string s, BifiltrationFormat input_format)
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

    std::ostream& operator<<(std::ostream& os, const Simplex& x)
    {
        os << "Simplex(id = " << x.id() << ", dim = " << x.dim();
        os << ", boundary = " << container_to_string(x.boundary()) << ", pos = " << x.position() << ")";
        return os;
    }
}
