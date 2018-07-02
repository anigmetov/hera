/*

Copyright (c) 2015, M. Kerber, D. Morozov, A. Nigmetov
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
(Enhancements) to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to copyright holder,
without imposing a separate written license agreement for such Enhancements,
then you hereby grant the following license: a  non-exclusive, royalty-free
perpetual license to install, use, modify, prepare derivative works, incorporate
into other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

  */

#ifndef AUCTION_ORACLE_BRUTEFORCE_H
#define AUCTION_ORACLE_BRUTEFORCE_H

#include <map>
#include <memory>
#include <set>

#include "basic_defs_ws.h"
#include "auction_oracle_base.h"


namespace hera {
namespace ws {

template <class Real_ = double, class PointContainer_ = std::vector<DiagramPoint<Real_>>>
struct AuctionOracleBruteforce : AuctionOracleBase<Real_, PointContainer_> {

    using PointContainer    = PointContainer_;
    using Real              = Real_;

    using DiagramPointR     = typename ws::DiagramPoint<Real>;
    using DebugOptimalBidR  = typename ws::DebugOptimalBid<Real>;

    AuctionOracleBruteforce(const PointContainer& bidders, const PointContainer& items, const AuctionParams<Real>& params);
    // data members
    // temporarily make everything public
    size_t num_bidders_;
    size_t num_items_;
    std::vector<std::vector<Real>> cost_matrix_;
    Real max_val_;

    // methods
    void set_price(const IdxType items_idx, const Real new_priced);
    IdxValPair<Real> get_optimal_bid(const IdxType bidder_idx);
    void adjust_prices();
    void adjust_prices(const Real delta);

    // debug routines
    DebugOptimalBidR get_optimal_bid_debug(IdxType bidder_idx) const;
    void sanity_check();

    std::shared_ptr<spdlog::logger> console_logger;

    std::pair<Real, Real> get_minmax_price() const;

};

} // ws
} // hera


#include "auction_oracle_bruteforce.hpp"

#endif
