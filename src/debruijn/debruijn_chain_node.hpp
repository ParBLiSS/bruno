/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * debruijn_chain_node.hpp
 *
 *  Created on: June 24, 2016
 *      Author: tony pan
 *
 */

#ifndef DEBRUIJN_CHAIN_NODE_HPP_
#define DEBRUIJN_CHAIN_NODE_HPP_

#include "bliss-config.hpp"

#include <tuple>        // tuple and utility functions

#include "common/alphabets.hpp"

namespace bliss
{
  namespace debruijn
  {
    // data type for debruijn graph chain compaction
    // first Kmer is in edge, second Kmer is out edge.  int is distance from end node.
    // 0 means that this node is a terminus.  negative numbers indicate that the edge is pointing to a terminus
    template <typename KMER>
    using simple_biedge = ::std::tuple<KMER, KMER, int, int>;


    namespace chain
    {

      /**
       * @brief a node in a compacted chain.
       *                      The KMER field points to the lex smaller terminus.
       *                      int is distance to that terminus.
       *                      uint8_t is the character in same encoding as KMER along the chain at that position.
       *
       *                      TODO: make into kmer for the last entry.
       */
      template <typename KMER>
      using compacted_chain_node = ::std::tuple<KMER, int, uint8_t>;


    }/*namespace chain*/


    namespace transform {

      /// standard companion function (used by lex_less) to compact_simple_biedge to get reverse complement of edge.
      template <typename KMER>
      inline ::bliss::debruijn::simple_biedge<KMER>
      reverse_complement(::bliss::debruijn::simple_biedge<KMER> const & x) {
        return ::bliss::debruijn::simple_biedge<KMER>(std::get<1>(x).reverse_complement(),
                                                      std::get<0>(x).reverse_complement(),
                                                      std::get<3>(x), std::get<2>(x));
      }


      template <typename KMER>
      inline ::bliss::debruijn::chain::compacted_chain_node<KMER>
      reverse_complement(::bliss::debruijn::chain::compacted_chain_node<KMER> const & x) {
        return ::bliss::debruijn::chain::compacted_chain_node<KMER>(std::get<0>(x).reverse_complement(),
                                                                    std::get<1>(x),
                                                                    KMER::KmerAlphabet::TO_COMPLEMENT[::std::get<2>(x)]);
      }


    } // namespace transform

  }/*namespace debruijn*/
}/*namespace bliss*/




#endif /* DEBRUIJN_CHAIN_NODE_HPP_ */
