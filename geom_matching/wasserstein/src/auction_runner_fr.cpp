/*

Copyright (c) 2016, M. Kerber, D. Morozov, A. Nigmetov
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


#include <assert.h>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iterator>
#include <chrono>

#include "def_debug_ws.h"
#include "auction_runner_fr.h"
#include "wasserstein.h"

#ifdef FOR_R_TDA
#include "Rcpp.h"
#endif

#define PRINT_DETAILED_TIMING

namespace geom_ws {

// *****************************
// AuctionRunnerFR
// *****************************

std::ostream& operator<<(std::ostream& o, const AuctionRunnerFR& f)
{
    o << "AuctionRunnerFR {\nnumBidders = " << f.numBidders << "\n";
    o << "numItems = " << f.numItems << "\n";

    auto bp = f.oracleReverse->getPrices();
    auto ip = f.oracleForward->getPrices();
    o << "bidders:\n";
    for(size_t b = 0; b < f.numBidders; ++b) {
        o << b << ": " << f.bidders[b] << ", price = " << bp[b] << "\n";
    }

    o << "items:\n";
    for(size_t i = 0; i < f.numItems; ++i) {
        o << i << ": " << f.items[i] << ", price = " << ip[i] << "\n";
    }

    o << "Assignment:\n";
    for(size_t b = 0; b < f.biddersToItems.size(); ++b) {
        auto i = f.biddersToItems[b];
        o << b << ", " << f.bidders[b] << " <--> " << i;
        if ( i >= 0 ) {
            o << ", " << f.items[i];
        }
        o << "\n";
    }

    o << "epsilon = " << f.oracleReverse->getEpsilon() << "\n";

    return o;
}

AuctionRunnerFR::AuctionRunnerFR(const std::vector<DiagramPoint>& A, const std::vector<DiagramPoint>& B, const double q, const double _delta, const double _internal_p, const double _initialEpsilon, const double _epsFactor) :
    bidders(A),
    items(B),
    numBidders(A.size()),
    numItems(A.size()),
    itemsToBidders(A.size(), -1),
    biddersToItems(A.size(), -1),
    wassersteinPower(q),
    delta(_delta),
    internal_p(_internal_p),
    initialEpsilon(_initialEpsilon),
    epsilonCommonRatio(_epsFactor == 0.0 ? 5.0 : _epsFactor)
{
    assert(initialEpsilon >= 0.0 );
    assert(epsilonCommonRatio >= 0.0 );
    assert(A.size() == B.size());
    oracleForward = std::unique_ptr<AuctionOracle>(new AuctionOracle(bidders, items, wassersteinPower, internal_p));
    weightAdjConst = oracleForward->maxVal;
    oracleReverse = std::unique_ptr<AuctionOracle>(new AuctionOracle(items, bidders, wassersteinPower, internal_p));
    std::cout << "Weight adj =  " << weightAdjConst << "\n";
}


void AuctionRunnerFR::assignItemToBidder(IdxType itemIdx, IdxType bidderIdx)
{
    numRounds++;
    numForwardRounds++;
    //sanityCheck();
    // only unassigned bidders should submit bids and get items
    assert(biddersToItems[bidderIdx] == -1);
    IdxType oldItemOwner = itemsToBidders[itemIdx];

    // set new owner
    biddersToItems[bidderIdx] = itemIdx;
    itemsToBidders[itemIdx] = bidderIdx;

    // remove bidder and item from the set of unassigned
    unassignedBidders.erase(bidderIdx);
    unassignedItems.erase(itemIdx);


    // old owner becomes unassigned
    if (oldItemOwner != -1) {
        biddersToItems[oldItemOwner] = -1;
        unassignedBidders.insert(oldItemOwner);
    }
}

void AuctionRunnerFR::assignBidderToItem(IdxType itemIdx, IdxType bidderIdx)
{
    numRounds++;
    numRevRounds++;
    //sanityCheck();
    // only unassigned items submit bids and get bidders in reverse rounds
    assert(itemsToBidders[itemIdx] == -1);
    IdxType oldBidderOwner = biddersToItems[bidderIdx];

    // set new owner
    biddersToItems[bidderIdx] = itemIdx;
    itemsToBidders[itemIdx] = bidderIdx;

    // remove bidder and item from the set of unassigned
    unassignedBidders.erase(bidderIdx);
    unassignedItems.erase(itemIdx);


    // old owner becomes unassigned
    if (oldBidderOwner != -1) {
        itemsToBidders[oldBidderOwner] = -1;
        unassignedItems.insert(oldBidderOwner);
    }
}


void AuctionRunnerFR::flushAssignment(void)
{
    for(auto& b2i : biddersToItems) {
        b2i = -1;
    }
    for(auto& i2b : itemsToBidders) {
        i2b = -1;
    }
    // we must flush assignment only after we got perfect matching
    assert(unassignedBidders.empty());
    // all bidders become unassigned
    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        unassignedBidders.insert(bidderIdx);
    }
    assert(unassignedBidders.size() == bidders.size());
    // all items become unassigned
    for(size_t itemIdx = 0; itemIdx < numItems; ++itemIdx) {
        unassignedItems.insert(itemIdx);
    }
    assert(unassignedItems.size() == items.size());
}

void AuctionRunnerFR::runAuction(void)
{
#ifdef PRINT_DETAILED_TIMING
    std::chrono::high_resolution_clock hrClock;
    std::chrono::time_point<std::chrono::high_resolution_clock> startMoment;
    startMoment = hrClock.now();
    std::vector<double> iterResults;
    std::vector<double> iterEstRelErrors;
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> iterTimes;
#endif
    // choose some initial epsilon
    if (initialEpsilon == 0.0) {
        double minMaxVal = std::min(oracleReverse->maxVal, oracleForward->maxVal);
        oracleForward->setEpsilon(minMaxVal / 4.0);
        oracleReverse->setEpsilon(minMaxVal / 4.0);
    } else {
        oracleForward->setEpsilon(initialEpsilon);
        oracleReverse->setEpsilon(initialEpsilon);
    }
    assert( oracleForward->getEpsilon() > 0 );
    assert( oracleReverse->getEpsilon() > 0 );

    int iterNum { 0 };
    bool notDone { false };
    double currentResult;
    do {
        flushAssignment();
        runAuctionPhase();
        iterNum++;
        currentResult = getDistanceToQthPowerInternal();
        double denominator = currentResult - numBidders * oracleForward->getEpsilon();
        currentResult = pow(currentResult, 1.0 / wassersteinPower);

#ifdef PRINT_DETAILED_TIMING
#ifndef FOR_R_TDA
        iterResults.push_back(currentResult);
        iterTimes.push_back(hrClock.now());
        std::cout << "Iteration " << iterNum << " finished. ";
        std::cout << "Current result is " << currentResult  << ", epsilon = " << oracleForward->getEpsilon() << std::endl;
        std::cout << "Number of rounds (cumulative): " << numRounds << std::endl;
        std::cout << "Number of forward rounds (cumulative): " << numForwardRounds << std::endl;
        std::cout << "Number of reverse rounds (cumulative): " << numRevRounds << std::endl;
#endif
#endif

        if ( denominator <= 0 ) {
            //std::cout << "Epsilon is too big." << std::endl;
            notDone = true;
        } else {
            denominator = pow(denominator, 1.0 / wassersteinPower);
            double numerator = currentResult - denominator;
#ifdef PRINT_DETAILED_TIMING
#ifndef FOR_R_TDA
            std::cout << " numerator: " << numerator << " denominator: " << denominator;
            std::cout << "; error bound: " << numerator / denominator << std::endl;
#endif
#endif
            // if relative error is greater than delta, continue
            notDone = ( numerator / denominator > delta );
        }

        // add the difference to bidder prices
        // so that after epsilon-scaling eps-CS is still satisfied
        auto oldEpsilon = oracleForward->getEpsilon();
        auto newEpsilon = oldEpsilon / epsilonCommonRatio;

        for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
            auto newPrice = 1.0000001 *(oldEpsilon - newEpsilon) + oracleReverse->getPrice(bidderIdx);
            oracleReverse->setPrice(bidderIdx, newPrice);
        }

        // decrease epsilon for the next iteration
        oracleForward->setEpsilon( oracleForward->getEpsilon() / epsilonCommonRatio );
        oracleReverse->setEpsilon( oracleReverse->getEpsilon() / epsilonCommonRatio );

        if (iterNum > maxIterNum) {
#ifndef FOR_R_TDA
            std::cerr << "Maximum iteration number exceeded, exiting. Current result is:";
            std::cerr << wassersteinDistance << std::endl;
#endif
            throw std::runtime_error("Maximum iteration number exceeded");
        }
    } while ( notDone );
}


void AuctionRunnerFR::runAuctionPhase()
{
    while(!unassignedItems.empty()) {
        //std::cout << *this;
        checkEpsilonCS();
        checkAssignmentConsistency();

        runForwardAuctionPhase();

        checkEpsilonCS();
        checkAssignmentConsistency();

        if (!unassignedBidders.empty()) {
            runReverseAuctionPhase();
        }
    }

    checkEpsilonCS();
}

void AuctionRunnerFR::runReverseAuctionPhase()
{
    const auto unassignedOrigQty = unassignedItems.size();
    const long int requiredNewAssignments = std::max(1LU, unassignedOrigQty / 10);

    //std::cout << "Entered runAuctionReversePhase" << std::endl;
    do {
        size_t itemIdx = *unassignedItems.begin();

        // get optimal bid
        auto optimalBid = oracleReverse->getOptimalBid(itemIdx);
        auto optimalBidderIdx = optimalBid.first;
        auto bidValue = optimalBid.second;

        assignBidderToItem(itemIdx, optimalBidderIdx);

        auto oldBidderPrice = oracleReverse->getPrice(optimalBidderIdx);
        auto oldItemPrice = oracleForward->getPrice(itemIdx);

        // set prices in both oracles
        oracleReverse->setPrice(optimalBidderIdx, bidValue);
        auto newItemPrice = getPairCost(optimalBidderIdx, itemIdx) - bidValue;
        oracleForward->setPrice(itemIdx, newItemPrice);

        assert(oldBidderPrice < bidValue);
        assert(oldItemPrice - newItemPrice >= -0.000001);

        if ( numRevRounds % 10000 == 0 ) {
            std::cout << "in Reverse ";
            std::cout << "numForwardRounds = " << numForwardRounds << ", numRevRounds = " << numRevRounds << ", total = " << numRounds << ", unassigned = " << unassignedItems.size() << std::endl;
        }


#ifdef FOR_R_TDA
        if ( numRounds % 10000 == 0 ) {
            Rcpp::checkUserInterrupt();
        }
#endif
    //} while (unassignedItems.size() >= unassignedOrigQty);
    } while (unassignedOrigQty - unassignedItems.size() < requiredNewAssignments);
    //} while (not unassignedBidders.empty());
    //std::cout << "runAuctionPhase finished" << std::endl;
}




void AuctionRunnerFR::runForwardAuctionPhase()
{
    const auto unassignedOrigQty = unassignedBidders.size();
    const long int requiredNewAssignments = std::max(1LU, unassignedOrigQty / 10);

    //std::cout << "Entered runAuctionForwardPhase" << std::endl;
    do {
        size_t bidderIdx = *unassignedBidders.begin();

        // get optimal bid
        auto optimalBid = oracleForward->getOptimalBid(bidderIdx);
        auto optimalItemIdx = optimalBid.first;
        auto bidValue = optimalBid.second;

        assignItemToBidder(optimalItemIdx, bidderIdx);



        auto oldBidderPrice = oracleReverse->getPrice(bidderIdx);
        auto oldItemPrice = oracleForward->getPrice(optimalItemIdx);

        // update prices in both oracles
        oracleForward->setPrice(optimalItemIdx, bidValue);
        auto newBidderPrice = getPairCost(bidderIdx, optimalItemIdx) - bidValue;
        oracleReverse->setPrice(bidderIdx, newBidderPrice);

        assert(oldItemPrice < bidValue);
        //if (oldBidderPrice <= newBidderPrice) {
            //std::cerr << "old bidder price = " << oldBidderPrice << ", new price = " << newBidderPrice << ", bidderIdx = " << bidderIdx << std::endl;
        //}
        assert(oldBidderPrice - newBidderPrice >= -0.00001);

        if ( numForwardRounds % 10000 == 0 ) {
            std::cout << "in Forward ";
            std::cout << "numForwardRounds = " << numForwardRounds << ", numRevRounds = " << numRevRounds << ", total = " << numRounds << ", unassigned = " << unassignedItems.size() << std::endl;
        }


#ifdef FOR_R_TDA
        if ( numRounds % 10000 == 0 ) {
            Rcpp::checkUserInterrupt();
        }
#endif
    //} while (unassignedBidders.size() >= unassignedOrigQty);
    } while (unassignedOrigQty - unassignedBidders.size() < requiredNewAssignments);
    //} while (not unassignedBidders.empty());
    //std::cout << "runAuctionPhase finished" << std::endl;
}

double AuctionRunnerFR::getDistanceToQthPowerInternal(void)
{
    sanityCheck();
    double result = 0.0;
    for(size_t bIdx = 0; bIdx < numBidders; ++bIdx) {
        auto pA = bidders[bIdx];
        assert( 0 <= biddersToItems[bIdx] and biddersToItems[bIdx] < static_cast<int>(items.size()) );
        auto pB = items[biddersToItems[bIdx]];
        // getPairCost returns minus the distance
        result += -getPairCost(pA, pB);
    }
    wassersteinCost = result;
    wassersteinDistance = pow(result, 1.0 / wassersteinPower);
    return result;
}

// return negated distance raised to wasssersteinPower
double AuctionRunnerFR::getPairCost(size_t bidderIdx, size_t itemIdx) const
{
    return getPairCost(bidders[bidderIdx], items[itemIdx]);
}

// return negated distance raised to wasssersteinPower
double AuctionRunnerFR::getPairCost(const DiagramPoint& pA, const DiagramPoint& pB) const
{
    return -pow(distLp(pA, pB, internal_p), wassersteinPower);
}

double AuctionRunnerFR::getWassersteinDistance(void)
{
    runAuction();
    return wassersteinDistance;
}

double AuctionRunnerFR::getWassersteinCost(void)
{
    runAuction();
    return wassersteinCost;
}



// Debug routines

void AuctionRunnerFR::printDebug(void)
{
#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    sanityCheck();
    std::cout << "**********************" << std::endl;
    std::cout << "Current assignment:" << std::endl;
    for(size_t idx = 0; idx < biddersToItems.size(); ++idx) {
        std::cout << idx << " <--> " << biddersToItems[idx] << std::endl;
    }

    std::cout << "Prices - forward: " << std::endl;
    for(const auto price : oracleForward->getPrices()) {
        std::cout << price << std::endl;
    }

    std::cout << "Prices - reverse: " << std::endl;
    for(const auto price : oracleReverse->getPrices()) {
        std::cout << price << std::endl;
    }
    std::cout << "**********************" << std::endl;
#endif
#endif
}

void AuctionRunnerFR::checkAssignmentConsistency()
{
#ifdef DEBUG_AUCTION
    size_t nAssignedBidders = 0;
    size_t nAssignedItems = 0;

    for(size_t b = 0; b < numBidders; ++b) {
        auto i = biddersToItems[b];
        if (i>= 0) {
            assert(unassignedBidders.count(b) == 0);
            assert((size_t)itemsToBidders[i] == b);
            nAssignedBidders++;
        } else {
            assert(unassignedBidders.count(b) == 1);
        }
    }

    for(size_t i = 0; i < numItems; ++i) {
        auto b = itemsToBidders[i];
        if (b>= 0) {
            assert(unassignedItems.count(i) == 0);
            assert((size_t)biddersToItems[b] == i);
            nAssignedItems++;
        } else {
            assert(unassignedItems.count(i) == 1);
        }
    }

    assert(nAssignedItems + unassignedItems.size() == numItems);
    assert(nAssignedBidders + unassignedBidders.size() == numBidders);
#endif
}

void AuctionRunnerFR::checkEpsilonCS()
{
#ifdef DEBUG_AUCTION
    std::vector<double> bPrices = oracleReverse->getPrices();
    std::vector<double> iPrices = oracleForward->getPrices();
    double eps = oracleReverse->getEpsilon();

    for(size_t b = 0; b < numBidders; ++b) {
        for(size_t i = 0; i < numItems; ++i) {
            if ((bidders[b].isDiagonal() && items[i].isNormal()) ||
                (bidders[b].isNormal() && items[i].isDiagonal()) )
                if (b != i)
                    continue;
            if ( bPrices[b] + iPrices[i] + eps < getPairCost(b, i) - 0.000001 ) {
                std::cerr << *this << std::endl;
                std::cerr << "eps-CS violated I, b = " << b << ", i = " << i << ", b price  = " << bPrices[b] << ", i Price = " << iPrices[i];
                std::cerr << ", cost = " << getPairCost(b, i) << ", eps = " << eps << std::endl;
                std::cerr << "numForwardRounds = " << numForwardRounds << ", numRevRounds = " << numRevRounds << std::endl;
                throw std::runtime_error("eps-CS violated I");
            }
        }
    }

    for(size_t b = 0; b < numBidders; ++b) {
        auto i = biddersToItems[b];
        if (i >= 0) {
            if ( bPrices[b] + iPrices[i] - getPairCost(b, i) > 0.000001 ) {
                std::cerr << "eps-CS violated II, b = " << b << ", i = " << i << ", b price  = " << bPrices[b] << ", i Price = " << iPrices[i] << std::endl;
                throw std::runtime_error("eps-CS violated II");
            }
        }
    }
#endif
}

void AuctionRunnerFR::sanityCheck(void)
{
#ifdef DEBUG_AUCTION
    if (unassignedItems.size() != unassignedBidders.size()) {
        throw std::runtime_error("Unassigned items and bidders mismatch");
    }

    if (biddersToItems.size() != numBidders) {
#ifndef FOR_R_TDA
        std::cerr << "Wrong size of biddersToItems, must be " << numBidders << ", is " << biddersToItems.size() << std::endl;
#endif
        throw std::runtime_error("Wrong size of biddersToItems");
    }

    if (itemsToBidders.size() != numBidders) {
#ifndef FOR_R_TDA
        std::cerr << "Wrong size of itemsToBidders, must be " << numBidders << ", is " << itemsToBidders.size() << std::endl;
#endif
        throw std::runtime_error("Wrong size of itemsToBidders");
    }

    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] >= 0) {

            if ( std::count(biddersToItems.begin(),
                        biddersToItems.end(),
                        biddersToItems[bidderIdx]) > 1 ) {
#ifndef FOR_R_TDA
                std::cerr << "Item " << biddersToItems[bidderIdx];
                std::cerr << " appears in biddersToItems more than once" << std::endl;
#endif
                throw std::runtime_error("Duplicate in biddersToItems");
            }

            if (itemsToBidders.at(biddersToItems[bidderIdx]) != static_cast<int>(bidderIdx)) {
#ifndef FOR_R_TDA
                std::cerr << "Inconsitency: bidderIdx = " << bidderIdx;
                std::cerr << ", itemIdx in biddersToItems = ";
                std::cerr << biddersToItems[bidderIdx];
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[biddersToItems[bidderIdx]] << std::endl;
#endif
                throw std::runtime_error("inconsistent mapping");
            }
        }
    }

    for(IdxType itemIdx = 0; itemIdx < static_cast<IdxType>(numBidders); ++itemIdx) {
        if ( itemsToBidders[itemIdx] >= 0) {

            // check for uniqueness
            if ( std::count(itemsToBidders.begin(),
                        itemsToBidders.end(),
                        itemsToBidders[itemIdx]) > 1 ) {
#ifndef FOR_R_TDA
                std::cerr << "Bidder " << itemsToBidders[itemIdx];
                std::cerr << " appears in itemsToBidders more than once" << std::endl;
#endif
                throw std::runtime_error("Duplicate in itemsToBidders");
            }
            // check for consistency
            if (biddersToItems.at(itemsToBidders[itemIdx]) != static_cast<int>(itemIdx)) {
#ifndef FOR_R_TDA
                std::cerr << "Inconsitency: itemIdx = " << itemIdx;
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[itemIdx];
                std::cerr << ", itemIdx in biddersToItems= ";
                std::cerr << biddersToItems[itemsToBidders[itemIdx]] << std::endl;
#endif
                throw std::runtime_error("inconsistent mapping");
            }
        }
    }
#endif
}

void AuctionRunnerFR::printMatching(void)
{
//#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    sanityCheck();
    for(size_t bIdx = 0; bIdx < biddersToItems.size(); ++bIdx) {
        if (biddersToItems[bIdx] >= 0) {
            auto pA = bidders[bIdx];
            auto pB = items[biddersToItems[bIdx]];
            std::cout <<  pA << " <-> " << pB << "+" << pow(distLp(pA, pB, internal_p), wassersteinPower) << std::endl;
        } else {
            assert(false);
        }
    }
#endif
//#endif
}



} // end of namespace geom_ws
