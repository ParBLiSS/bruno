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
 * debruijn_common.hpp
 *
 *  Created on: June 24, 2016
 *      Author: tony pan
 *
 */

#ifndef DEBRUIJN_COMMON_HPP_
#define DEBRUIJN_COMMON_HPP_

#include "bliss-config.hpp"

// #include "utils/logging.h"
// #include "utils/transform_utils.hpp"

#include "debruijn/debruijn_biedge_loader.hpp"
#include "debruijn/debruijn_chain_operations.hpp"
#include "debruijn/debruijn_chain_node.hpp"

// this file should be included AFTER the types used by lex_less are included, but before the map using lex_less and MapParams

namespace bliss
{
  namespace debruijn
  {

    template <typename KMER>
    struct lex_less {
        inline KMER operator()(KMER const & x) const  {
          auto y = x.reverse_complement();
          return (x < y) ? x : y;
        }

        // standard operator to get canonical k-mer and value.  value is also reversed and complemented.
        template <typename VAL>
        inline ::std::pair<KMER, VAL > operator()(std::pair<KMER, VAL > const & x) const {
          auto y = x.first.reverse_complement();
          return (x.first < y) ? x :   // if already canonical, just return input
              std::pair<KMER, VAL >(y, ::bliss::debruijn::transform::reverse_complement(x.second) );
        }

        // standard operator to get canonical k-mer and value.  value is also reversed and complemented.
        template <typename VAL>
        inline ::std::pair<KMER, VAL > operator()(std::pair<const KMER, VAL > const & x) const {
          auto y = x.first.reverse_complement();
          return (x.first < y) ? x :   // if already canonical, just return input
              std::pair<KMER, VAL >(y, ::bliss::debruijn::transform::reverse_complement(x.second) );
        }

    };


    template <typename Kmer, template <typename> class DistHash = ::bliss::index::kmer::DistHashMurmur >
    using CanonicalDeBruijnHashMapParams = ::dsc::HashMapParams<
        Kmer,
        ::bliss::debruijn::lex_less,  // precanonalizer.  OPERATES ON VALUE AS WELL
         ::bliss::transform::identity,  // only one that makes sense given InputTransform
          DistHash,
          ::std::equal_to,
           ::bliss::transform::identity,
            ::bliss::index::kmer::StoreHashMurmur,
            ::std::equal_to
          >;

    template <typename Kmer >
    using CanonicalDeBruijnHashMapParamsMurmur = CanonicalDeBruijnHashMapParams<Kmer, ::bliss::index::kmer::DistHashMurmur >;

  }/*namespace debruijn*/
}/*namespace bliss*/




#endif /* DEBRUIJN_CHAIN_NODE_HPP_ */
