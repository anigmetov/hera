#ifndef MATCHING_DISTANCE_BIFILTRATION_H
#define MATCHING_DISTANCE_BIFILTRATION_H

#include <string>
#include <ostream>

#include "common_util.h"
#include "box.h"
#include "simplex.h"
#include "dual_point.h"

namespace md {

    class Bifiltration {
    public:
        using Diagram = std::vector<std::pair<Real, Real>>;
        using Box = md::Box;
        using SimplexVector = std::vector<Simplex>;

        Bifiltration() = default;

        Bifiltration(const Bifiltration&) = default;

        Bifiltration(Bifiltration&&) = default;

        Bifiltration& operator=(const Bifiltration& other)& = default;

        Bifiltration& operator=(Bifiltration&& other) = default;

        Bifiltration(const std::string& fname,
                BifiltrationFormat input_format = BifiltrationFormat::rivet); // read from file


        template<class Iter>
        Bifiltration(Iter begin, Iter end)
                : simplices_(begin, end)
        {
            init();
        }

        Diagram weighted_slice_diagram(const DualPoint& line, int dim) const;

        SimplexVector simplices() const { return simplices_; }

        // translate all points by vector (a,a)
        void translate(Real a);

        // return minimal value of x- and y-coordinates
        // among all simplices
        Real minimal_coordinate() const;

        // return box that contains positions of all simplices
        Box bounding_box() const;

        void sanity_check() const;

        int maximal_dim() const { return maximal_dim_; }

        friend std::ostream& operator<<(std::ostream& os, const Bifiltration& bif);

        Real max_x() const;

        Real max_y() const;

        Real min_x() const;

        Real min_y() const;

        void add_simplex(Index _id, Point birth, int _dim, const Column& _bdry);

        void save(const std::string& filename, BifiltrationFormat format = BifiltrationFormat::rivet); // save to file

        void scale(Real lambda);

    private:
        SimplexVector simplices_;
        // axes names, for rivet bifiltration format only
        std::string parameter_1_name_ {"axis_1"};
        std::string parameter_2_name_ {"axis_2"};

        Box bounding_box_;
        int maximal_dim_ {-1};

        void init();

        void rivet_format_reader(std::ifstream&);

        void rene_format_reader(std::ifstream&);

        // in Rene format each simplex knows IDs of its boundary facets
        // postprocess_rene_format fills vector of IDs of boundary facets
        // in each simplex
        void postprocess_rene_format();

        // in Rivet format each simplex knows its vertices,
        // postprocess_rivet_format fills vector of IDs of boundary facets
        // in each simplex
        void postprocess_rivet_format();

    };

    std::ostream& operator<<(std::ostream& os, const Bifiltration& bif);

    class BifiltrationProxy {
    public:
        BifiltrationProxy(const Bifiltration& bif, int dim = 0);
        // return critical values of simplices that are important for current dimension (dim and dim+1)
        PointVec positions() const;
        // set current dimension
        int set_dim(int new_dim);

        // wrappers of Bifiltration
        int maximal_dim() const;
        void translate(Real a);
        Real minimal_coordinate() const;
        Box bounding_box() const;
        Real max_x() const;
        Real max_y() const;
        Real min_x() const;
        Real min_y() const;
        Diagram weighted_slice_diagram(const DualPoint& slice) const;

    private:
        int dim_ { 0 };
        mutable PointVec cached_positions_;
        Bifiltration bif_;

        void cache_positions() const;
    };
}



#endif //MATCHING_DISTANCE_BIFILTRATION_H

//// The value type of OutputIterator is Simplex_in_2D_filtration
//template<typename OutputIterator>
//void read_input(std::string filename, OutputIterator out)
//{
//    std::ifstream ifstr;
//    ifstr.open(filename.c_str());
//    long n;
//    ifstr >> n; // number of simplices is the first number in file
//
//    Index k;   // used in loop
//    for (int i = 0; i < n; i++) {
//        Simplex_in_2D_filtration next;
//        next.index = i;
//        ifstr >> next.dim >> next.pos.x >> next.pos.y;
//        if (next.dim > 0) {
//            for (int j = 0; j <= next.dim; j++) {
//                ifstr >> k;
//                next.bd.push_back(k);
//            }
//        }
//        *out++ = next;
//    }
//}
