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
 * de_bruijn_operations.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DE_BRUIJN_OPERATIONS_HPP_
#define DE_BRUIJN_OPERATIONS_HPP_

#include <tuple>

namespace bliss {
  namespace de_bruijn {
    namespace operation {

      static constexpr char IN = -1;
      static constexpr char OUT = 1;

      namespace chain {

        // data type for debruijn graph chain compaction
        // first Kmer is in edge, second Kmer is out edge.  int is distance from end node.
        template <typename Kmer>
        using compaction_metadata = ::std::tuple<Kmer, Kmer, int, int>;

        // k-mer is the source k-mer ON SAME STRAND as the canonical key kmer that this is used with.
        // int is the indicator for whether the Kmer is an in or out edge from key.
        template <typename Kmer>
        using terminus_update_md = ::std::pair<Kmer, char>;

        template <typename Kmer>
        struct terminus_update {
            void operator()(compaction_metadata<Kmer> & x,
                            terminus_update_md<Kmer> const & y)  {
              if (y.second == IN) {  // sense
                if ((std::get<2>(x) == 0) || (y.first != std::get<0>(x))) {
                  std::cout << "y: " << y.first << ", y ?= x[in] " << std::get<2>(x) << ": " <<   std::get<0>(x) << ", x[out] " << std::get<3>(x) << ": " << std::get<1>(x) << std::endl;
                }

                assert(std::get<2>(x) != 0);   // just to make sure that the edge we're trying to update is not already "cut".
                assert(y.first == std::get<0>(x));   // just to make sure that we are talking about the same kmer.

                std::get<2>(x) = 0;   // mark as terminus
              } else {
                if ((std::get<3>(x) == 0) || (y.first != std::get<1>(x))) {
                  std::cout << "y: " << y.first << ", x[in] " << std::get<2>(x) << ": " <<   std::get<0>(x) << ", y ?= x[out] " << std::get<3>(x) << ": " << std::get<1>(x) <<  std::endl;
                }

                assert(std::get<3>(x) != 0);   // just to make sure that the edge we're trying to update is not already "cut".
                assert(y.first == std::get<1>(x));   // just to make sure that we are talking about the same kmer.
                std::get<3>(x) = 0;   // mark as terminus
              }
            }
        };


        template <typename KMER>
        struct lex_less {
            inline KMER operator()(KMER const & x) const  {
              auto y = x.reverse_complement();
              return (x < y) ? x : y;
            }
            inline KMER operator()(KMER const & x, KMER const & rc) const  {
              return (x < rc) ? x : rc;
            }

            inline ::std::pair<KMER, ::bliss::de_bruijn::operation::chain::terminus_update_md<KMER> >
            operator()(std::pair<KMER, ::bliss::de_bruijn::operation::chain::terminus_update_md<KMER> > const & x) const {
              auto y = x.first.reverse_complement();
              if (x.first < y)
            	  return x;   // if already canonical, just return input
              else {
            	  ::bliss::de_bruijn::operation::chain::terminus_update_md<KMER> z(x.second.first.reverse_complement(), -(x.second.second));
            	  return  std::pair<KMER, ::bliss::de_bruijn::operation::chain::terminus_update_md<KMER> >(
                          y, z );

              }
            }
        };

  	  template <typename Kmer >
  	  using CanonicalDeBruijnChainMapParams = ::dsc::HashMapParams<
  	      Kmer,
  	      ::bliss::de_bruijn::operation::chain::lex_less,  // precanonalizer.  operates on the value as well
  	       ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
  	        ::bliss::index::kmer::DistHashFarm,
  	        ::std::equal_to,
  	         ::bliss::kmer::transform::identity,
  	          ::bliss::index::kmer::StoreHashFarm,
  	          ::std::equal_to
  	        >;

      } // namespace chain

    } //namespace operation
  } //namespace de_bruijn
} //namespace bliss




#endif // DE_BRUIJN_OPERATIONS_HPP_
