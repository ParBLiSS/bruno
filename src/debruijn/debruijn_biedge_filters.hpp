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
 * debruijn_biedge_filters.hpp
 *
 *  Created on: Sept 26, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_BIEDGE_FILTERS_HPP_
#define DEBRUIJN_BIEDGE_FILTERS_HPP_

#include <vector>
#include <debruijn/debruijn_graph_node.hpp>

namespace bliss {
  namespace debruijn {
    namespace filter {
      namespace biedge {

        /**
         * @brief filter biedges
         * @note  while the parser has a predicate, frequency is distributed so we should not look up k+1 mer frequency on the fly one by one.
         * @param pred  unary operator with no return, that operates on the input iterator (must be mutable), and changes the input directly.
         */
        template <typename KmerBiedgeIter, typename CondTransform>
        void filter_biedges(KmerBiedgeIter start, KmerBiedgeIter end, CondTransform const &pred = CondTransform()) {
          std::for_each(start, end, pred);
        }
        /**
         * @brief filter biedges
         * @note  while the parser has a predicate, frequency is distributed so we should not look up k+1 mer frequency on the fly one by one.
         * @param pred  unary operator with return, that operates on the input iterator and returns the transformed output.
         */
        template <typename KmerBiedgeIter, typename CondTransform>
        void filter_copy_biedges(KmerBiedgeIter start, KmerBiedgeIter end, KmerBiedgeIter out, CondTransform const &pred = CondTransform()) {
          std::transform(start, end, out, pred);
        }


        // TODO:
        // filter by edge frequency:  input: vector of kmer-biedge, conditional transform function
        //    query by kmer for frequency/presence first in bulk
        //    then build local count hashtable in bulk
        //    finally apply conditional transform.


        /**
         * @brief filter a single biedge based on frequency of the edges.
         * @note  a similar filter for node frequency can be constructed as well.
         *        operates on pair of kmer + simple biedge.
         *
         *        not a traditional filter with bool output.  instead, this should work with for_each, and directly modify the compact_simple_biedge.
         */
        template <typename KmerType, typename CountIndexType>
        struct edge_frequency_inplace_filter {
            CountIndexType & counter;

            edge_frequency_inplace_filter(CountIndexType const & _counter) : counter(_counter) {};


            void operator()(std::pair<KmerType, ::bliss::debruijn::compact_simple_biedge> const & x) {
              ::bliss::debruijn::get_in_edge_k1mer()
            }
        };


        /**
         * @brief filter biedges
         * @note  while the parser has a predicate, frequency is distributed so we should not look up k+1 mer frequency on the fly one by one.
         *
         */
        template <typename KmerBiedgeIter, typename Predicate = ::fsc::identity>
        void filter_biedge(KmerBiedgeIter start, KmerBiedgeIter end, Predicate const &pred = Predicate()) {
          std::for_each(start, end, pred);
        }





      } // namespace graph
    } //namespace filter
  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_GRAPH_FILTERS_HPP_
