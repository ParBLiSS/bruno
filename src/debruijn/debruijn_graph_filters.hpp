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
 * debruijn_graph_filters.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_GRAPH_FILTERS_HPP_
#define DEBRUIJN_GRAPH_FILTERS_HPP_

#include <vector>

namespace bliss {
  namespace debruijn {
    namespace filter {
      namespace graph {

        //============= predicates for global filter operations.

        struct IsBranchPoint {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (t.second.get_in_edge_count() >= 2) || (t.second.get_out_edge_count() >= 2);
            }
        };

        struct IsIsolated {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (t.second.get_in_edge_count() == 0) && (t.second.get_out_edge_count() == 0);
            }
        };

        struct IsTerminus {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last)  const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const  {
              uint8_t in = t.second.get_in_edge_count();
              uint8_t out = t.second.get_out_edge_count();

              return ((in == 0) && (out == 1)) ||
                  ((in == 1) && (out == 0));
            }
        };

        struct IsChainInternalNode {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const  {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t)  const {
              return (t.second.get_in_edge_count() == 1) && (t.second.get_out_edge_count() == 1);
            }
        };


        struct IsChainNode {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last)  const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t)  const {
              uint8_t in = t.second.get_in_edge_count();
              uint8_t out = t.second.get_out_edge_count();

              return ((in == 0) && (out == 1)) ||
                  ((in == 1) && (out == 0)) ||
                  ((in == 1) && (out == 1));
            }
        };

      } // namespace graph
    } //namespace filter
  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_GRAPH_FILTERS_HPP_