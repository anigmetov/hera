#ifndef MATCHING_DISTANCE_SIMPLEX_H
#define MATCHING_DISTANCE_SIMPLEX_H

#include <algorithm>
#include <vector>
#include <ostream>

#include "common_util.h"

namespace md {

    class Bifiltration;

    enum class BifiltrationFormat {
        rene, rivet
    };

    class AbstractSimplex {
    private:
        std::vector<int> vertices_;
    public:

        // this member is for convenience only;
        // abstract simplices are identified  by their set of vertices
        mutable int id {-1};

        decltype(auto) begin() { return vertices_.begin(); }

        decltype(auto) end() { return vertices_.end(); }

        decltype(auto) begin() const { return vertices_.begin(); }

        decltype(auto) end() const { return vertices_.end(); }

        decltype(auto) cbegin() const { return vertices_.cbegin(); }

        decltype(auto) cend() const { return vertices_.cend(); }

        int dim() const { return vertices_.size() - 1; }

        void push_back(int v);

        AbstractSimplex() { }

        AbstractSimplex(std::vector<int> vertices, bool sort = true);

        template<class Iter>
        AbstractSimplex(Iter beg_iter, Iter end_iter, bool sort = true)
                :
                vertices_(beg_iter, end_iter)
        {
            if (sort)
                std::sort(vertices_.begin(), end());
        }

        std::vector<AbstractSimplex> facets() const;

        friend std::ostream& operator<<(std::ostream& os, const AbstractSimplex& s);

        // compare by vertices_ only
        friend bool operator==(const AbstractSimplex& s1, const AbstractSimplex& s2);

        friend bool operator<(const AbstractSimplex&, const AbstractSimplex&);
    };

    std::ostream& operator<<(std::ostream& os, const AbstractSimplex& s);

    class Simplex {
    private:
        Index id_;
        Point pos_;
        int dim_;
        // in our format we use facet indices,
        // this is the fastest representation for homology
        // Rivet format fills vertices_ vector
        // Simplex alone cannot convert from one representation to the other,
        // conversion routines are in Bifiltration
        Column facet_indices_;
        Column vertices_;
        Real v {0.0}; // used when constructed a filtration for a slice
    public:
        Simplex(Index _id, std::string s, BifiltrationFormat input_format);

        Simplex(Index _id, Point birth, int _dim, const Column& _bdry);

        void init_rivet(std::string s);

        void init_rene(std::string s);

        Index id() const { return id_; }

        int dim() const { return dim_; }

        Column boundary() const { return facet_indices_; }

        Real value() const { return v; }

        // assumes 1-criticality
        Point position() const { return pos_; }

        void set_position(const Point& new_pos) { pos_ = new_pos; }

        void scale(Real lambda)
        {
            pos_.x *= lambda;
            pos_.y *= lambda;
        }

        void translate(Real a);

        void set_value(Real new_val) { v = new_val; }

        friend std::ostream& operator<<(std::ostream& os, const Simplex& s);

        friend Bifiltration;
    };

    std::ostream& operator<<(std::ostream& os, const Simplex& s);

}
#endif //MATCHING_DISTANCE_SIMPLEX_H
