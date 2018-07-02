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
#ifndef AUCTION_ORACLE_BRUTEFORCE_HPP
#define AUCTION_ORACLE_BRUTEFORCE_HPP

#include <assert.h>
#include <algorithm>
#include <functional>
#include <iterator>

#include "def_debug_ws.h"
#include "basic_defs_ws.h"
#include "auction_oracle_bruteforce.h"

#ifdef FOR_R_TDA
#undef DEBUG_AUCTION
#endif

namespace hera {
namespace ws {


template<class Real_, class PointContainer_>
AuctionOracleBruteforce<Real_, PointContainer_>::AuctionOracleBruteforce(const PointContainer_& _bidders,
                                                                   const PointContainer_& _items,
                                                                   const AuctionParams<Real_>& params) :
    AuctionOracleBase<Real_, PointContainer_>(_bidders, _items, params),
    num_bidders_(_bidders.size()),
    num_items_(_items.size()),
    cost_matrix_(num_bidders_, std::vector<Real_>(num_items_, 0.0))
{
    max_val_ = 0.0;
    for(size_t bidder_idx = 0; bidder_idx < num_bidders_; ++bidder_idx) {
        for(size_t item_idx = 0; item_idx < num_items_; ++item_idx) {
            Real curr_val = this->get_value_for_bidder(bidder_idx, item_idx);
            cost_matrix_[bidder_idx][item_idx] = curr_val;
            max_val_ = std::max(max_val_, curr_val);
        }
    }

    console_logger = spdlog::get("console");
    if (not console_logger) {
        console_logger = spdlog::stdout_logger_st("console");
    }
}


template<class Real_, class PointContainer_>
typename AuctionOracleBruteforce<Real_, PointContainer_>::DebugOptimalBidR
AuctionOracleBruteforce<Real_, PointContainer_>::get_optimal_bid_debug(IdxType bidder_idx) const
{
    throw std::runtime_error("Not implemented");
}


template<class Real_, class PointContainer_>
IdxValPair<Real_> AuctionOracleBruteforce<Real_, PointContainer_>::get_optimal_bid(IdxType bidder_idx)
{

    size_t best_item_idx { k_invalid_index };
    //  size_t second_best_item_idx { k_invalid_index };
    Real best_item_value { std::numeric_limits<Real>::max() };
    Real second_best_item_value { std::numeric_limits<Real>::max() };

    for(size_t item_idx = 0; item_idx < num_items_; ++item_idx) {
        Real curr_val = cost_matrix_[bidder_idx][item_idx] + this->prices[item_idx];
        if (curr_val < best_item_value) {
            best_item_idx = item_idx;
            best_item_value = curr_val;
        }
    }

    for(size_t item_idx = 0; item_idx < num_items_; ++item_idx) {
        if (item_idx == best_item_idx)
            continue;
        Real curr_val = cost_matrix_[bidder_idx][item_idx] + this->prices[item_idx];
        if (curr_val < second_best_item_value) {
            // second_best_item_idx = item_idx;
            second_best_item_value = curr_val;
        }
    }

    IdxValPair<Real> result;

    assert( second_best_item_value >= best_item_value );

    result.first = best_item_idx;
    result.second = ( second_best_item_value - best_item_value ) + this->prices[best_item_idx] + this->epsilon;

    return result;
}


template<class Real_, class PointContainer_>
void AuctionOracleBruteforce<Real_, PointContainer_>::set_price(IdxType item_idx, Real new_price)
{
    assert(this->prices.size() == this->items.size());
    this->prices[item_idx] = new_price;
}


template<class Real_, class PointContainer_>
void AuctionOracleBruteforce<Real_, PointContainer_>::adjust_prices(Real delta)
{
    if (delta == 0.0)
        return;

    for(auto& p : this->prices) {
        p -= delta;
    }
}

template<class Real_, class PointContainer_>
void AuctionOracleBruteforce<Real_, PointContainer_>::adjust_prices()
{
    auto pr_begin = this->prices.begin();
    auto pr_end = this->prices.end();
    Real min_price = *(std::min_element(pr_begin, pr_end));
    adjust_prices(min_price);
}

template<class Real_, class PointContainer_>
std::pair<Real_, Real_> AuctionOracleBruteforce<Real_, PointContainer_>::get_minmax_price() const
{
    auto r = std::minmax_element(this->prices.begin(), this->prices.end());
    return std::make_pair(*r.first, *r.second);
}

template<class Real_, class PointContainer_>
void AuctionOracleBruteforce<Real_, PointContainer_>::sanity_check()
{
}


} // ws
} // hera

#endif
