#ifndef MATCHING_DISTANCE_BIFILTRATION_H
#define MATCHING_DISTANCE_BIFILTRATION_H

#include <string>
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "common_util.h"
#include "box.h"
#include "simplex.h"
#include "dual_point.h"
#include "phat/boundary_matrix.h"
#include "phat/compute_persistence_pairs.h"

#include "common_util.h"

namespace md {

    template<class Real>
    class Bifiltration {
    public:
        using SimplexVector = std::vector<Simplex<Real>>;

        Bifiltration() = default;

        Bifiltration(const Bifiltration&) = default;

        Bifiltration(Bifiltration&&) = default;

        Bifiltration& operator=(const Bifiltration& other)& = default;

        Bifiltration& operator=(Bifiltration&& other) = default;

        Bifiltration(const std::string& fname); // read from file

        template<class Iter>
        Bifiltration(Iter begin, Iter end)
                : simplices_(begin, end)
        {
            init();
        }

        Diagram<Real> weighted_slice_diagram(const DualPoint<Real>& line, int dim) const;

        SimplexVector simplices() const { return simplices_; }

        // translate all points by vector (a,a)
        void translate(Real a);

        // return minimal value of x- and y-coordinates
        // among all simplices
        Real minimal_coordinate() const;

        // return box that contains positions of all simplices
        Box<Real> bounding_box() const;

        void sanity_check() const;

        int maximal_dim() const { return maximal_dim_; }

        Real max_x() const;

        Real max_y() const;

        Real min_x() const;

        Real min_y() const;

        void add_simplex(Index _id, Point<Real> birth, int _dim, const Column& _bdry);

        void save(const std::string& filename, BifiltrationFormat format = BifiltrationFormat::rivet); // save to file

        void scale(Real lambda);

    private:
        SimplexVector simplices_;

        Box<Real> bounding_box_;
        int maximal_dim_ {-1};

        void init();

        void rivet_format_reader(std::ifstream&);

        void phat_like_format_reader(std::ifstream&);

        // in Rene format each simplex knows IDs of its boundary facets
        // postprocess_phat_like_format fills vector of IDs of boundary facets
        // in each simplex
        void postprocess_phat_like_format();

        // in Rivet format each simplex knows its vertices,
        // postprocess_rivet_format fills vector of IDs of boundary facets
        // in each simplex
        void postprocess_rivet_format();

    };

    template<class Real>
    std::ostream& operator<<(std::ostream& os, const Bifiltration<Real>& bif);

    template<class Real>
    class BifiltrationProxy {
    public:
        BifiltrationProxy(const Bifiltration<Real>& bif, int dim = 0);
        // return critical values of simplices that are important for current dimension (dim and dim+1)
        PointVec<Real> positions() const;
        // set current dimension
        int set_dim(int new_dim);

        // wrappers of Bifiltration
        int maximal_dim() const;
        void translate(Real a);
        Real minimal_coordinate() const;
        Box<Real> bounding_box() const;
        Real max_x() const;
        Real max_y() const;
        Real min_x() const;
        Real min_y() const;
        Diagram<Real> weighted_slice_diagram(const DualPoint<Real>& slice) const;

    private:
        int dim_ { 0 };
        mutable PointVec<Real> cached_positions_;
        Bifiltration<Real> bif_;

        void cache_positions() const;
    };
}

#include "bifiltration.hpp"

#endif //MATCHING_DISTANCE_BIFILTRATION_H
