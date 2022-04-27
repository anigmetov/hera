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
    Real wasserstein_power {1};
    Real delta {0.01}; // relative error
    Real internal_p {get_infinity<Real>()};
    Real initial_epsilon {0}; // 0.0 means maxVal / 4.0
    Real epsilon_common_ratio {5};
    int max_num_phases {std::numeric_limits<decltype(max_num_phases)>::max()};
    int max_bids_per_round {1};  // imitate Gauss-Seidel is default behaviour
    unsigned int dim {2}; // for pure geometric version only; ignored in persistence diagrams
    bool tolerate_max_iter_exceeded {false}; // whether auction should throw an exception on max. iterations exceeded
    bool return_matching {false}; // whether to return optimal matching along with cost
    bool match_inf_points {true}; // whether to add infinite points to matching; ignored, if return_matching is false

};

} // namespace hera

#endif //HERA_AUCTION_PARAMS_H
