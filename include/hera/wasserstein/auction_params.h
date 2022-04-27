#ifndef HERA_AUCTION_PARAMS_H
#define HERA_AUCTION_PARAMS_H

#include <unordered_map>
#include <limits>
#include <cassert>
#include <iostream>

#include <hera/common.h>

namespace hera {

template<class Real = double>
struct AuctionParams {
    Real wasserstein_power {1.0};
    Real delta {0.01}; // relative error
    Real internal_p {get_infinity<Real>()};
    Real initial_epsilon {0.0}; // 0.0 means maxVal / 4.0
    Real epsilon_common_ratio {5.0};
    int max_num_phases {std::numeric_limits<decltype(max_num_phases)>::max()};
    int max_bids_per_round {1};  // imitate Gauss-Seidel is default behaviour
    unsigned int dim {2}; // for pure geometric version only; ignored in persistence diagrams
    Real final_relative_error;  // out parameter - after auction terminates, contains the real relative error
    bool tolerate_max_iter_exceeded {false}; // whether auction should throw an exception on max. iterations exceeded
    bool return_matching {false}; // whether to return optimal matching along with cost
    bool match_inf_points {true}; // whether to add infinite points to matching; ignored, if return_matching is false

    std::unordered_map<int, int> matching_a_to_b_;
    std::unordered_map<int, int> matching_b_to_a_;

    void clear_matching()
    {
        matching_a_to_b_.clear();
        matching_b_to_a_.clear();
    }

    void add_to_matching(int a, int b)
    {
        assert(matching_a_to_b_.count(a) == 0 and matching_b_to_a_.count(b) == 0);
        matching_a_to_b_[a] = b;
        matching_b_to_a_[b] = a;
    }
};

} // namespace hera

#endif //HERA_AUCTION_PARAMS_H
