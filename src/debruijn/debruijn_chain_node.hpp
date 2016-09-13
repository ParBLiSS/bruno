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
#include "common/kmer.hpp"

namespace bliss
{
  namespace debruijn
  {
    // data type for debruijn graph chain compaction
    // first Kmer is in edge, second Kmer is out edge.  int is distance from end nodes.
    // 0 means that this node is a terminus.  negative numbers indicate that the edge is pointing to a terminus
    // used with Kmer key in a map, to mean in-Kmer-out bi edge.
    template <typename KMER>
    using simple_biedge = ::std::tuple<KMER, KMER, int, int>;
    //static_assert(sizeof(simple_biedge<::bliss::common::Kmer<31, ::bliss::common::DNA> >) == sizeof(::bliss::common::Kmer<31, ::bliss::common::DNA>) * 2 + 2 * sizeof(int), "size of simple biedge is not what's expected");


    namespace chain
    {

      /**
       * @brief a node in a compacted chain.
       *
       * @tparam KMER 		first KMER field is the current k-mer, canonical
       * @tparam KMER		second KMER field is the lex smaller 5' terminus of the chain
       * @tparam int		int is distance to that terminus.  +/0 indicates both are on same strand.  - means otherwise.
       *
       *      note that a vector of these are sufficient to represent a complete compacted chain.
       */
      template <typename KMER>
      using compacted_chain_node = ::std::tuple<KMER, KMER, int>;

      //static_assert(sizeof(compacted_chain_node<::bliss::common::Kmer<31, ::bliss::common::DNA> >) == sizeof(::bliss::common::Kmer<31, ::bliss::common::DNA>) * 2 + sizeof(int), "size of compacted chain node is not what's expected");
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


    } // namespace transform

  }/*namespace debruijn*/
}/*namespace bliss*/




#endif /* DEBRUIJN_CHAIN_NODE_HPP_ */
