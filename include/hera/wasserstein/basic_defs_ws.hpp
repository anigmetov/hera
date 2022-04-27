/*

Copyright (c) 2015, M. Kerber, D. Morozov, A. Nigmetov
Copyright (c) 2018, G. Spreemann
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

#pragma once

#include <algorithm>
#include <cfloat>
#include <set>
#include <algorithm>
#include <istream>
#include <cstdint>
#include "basic_defs_ws.h"

#ifndef FOR_R_TDA
#include <iostream>
#endif

#include <sstream>

namespace hera {
static const int64_t DIPHA_MAGIC = 8067171840;
static const int64_t DIPHA_PERSISTENCE_DIAGRAM = 2;

namespace ws {






template<class Container>
inline std::string format_container_to_log(const Container& cont)
{
    std::stringstream result;
    result << "[";
    for(auto iter = cont.begin(); iter != cont.end(); ++iter) {
        result << *iter;
        if (std::next(iter) != cont.end()) {
            result << ", ";
        }
    }
    result << "]";
    return result.str();
}

template<class Container>
inline std::string format_pair_container_to_log(const Container& cont)
{
    std::stringstream result;
    result << "[";
    for(auto iter = cont.begin(); iter != cont.end(); ++iter) {
        result << "(" << iter->first << ", " << iter->second << ")";
        if (std::next(iter) != cont.end()) {
            result << ", ";
        }
    }
    result << "]";
    return result.str();
}


template<class Real, class IndexContainer>
inline std::string format_point_set_to_log(const IndexContainer& indices,
                                    const std::vector<DiagramPoint<Real>>& points)
{
    std::stringstream result;
    result << "[";
    for(auto iter = indices.begin(); iter != indices.end(); ++iter) {
        DiagramPoint<Real> p = points[*iter];
        result << "(" << p.getRealX() << ", " << p.getRealY() << ")";
        if (std::next(iter) != indices.end())
            result << ", ";
    }
    result << "]";
    return result.str();
}

template<class T>
inline std::string format_int(T i)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << i;
    return ss.str();
}


} // end of namespace ws


template <typename T> inline void reverse_endianness(T & x)
{
    uint8_t * p = reinterpret_cast<uint8_t *>(&x);
    std::reverse(p, p + sizeof(T));
}

template <typename T> inline T read_le(std::istream & s)
{
    T result;
    s.read(reinterpret_cast<char *>(&result), sizeof(T));
    #ifdef BIGENDIAN
    reverse_endianness(result);
    #endif
    return result;
}

} // hera