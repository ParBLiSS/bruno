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
 * de_bruijn_filter.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DE_BRUIJN_FILTER_HPP_
#define DE_BRUIJN_FILTER_HPP_

#include <vector>

namespace bliss {
  namespace de_bruijn {
    namespace filter {

      //============= predicates for global filter operations.

      struct IsBranchPoint {
          /// does not filter by group of results.
          template <typename Iter>
          inline bool operator()(Iter first, Iter last) const {  return true; }

          template <typename Kmer, typename Edge>
          inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
            uint8_t in = t.second.get_in_edge_count();
            uint8_t out = t.second.get_out_edge_count();

            return (in >= 2) || (out >= 2);
          }
      };

      struct IsIsolated {
          /// does not filter by group of results.
          template <typename Iter>
          inline bool operator()(Iter first, Iter last) const {  return true; }

          template <typename Kmer, typename Edge>
          inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
            uint8_t in = t.second.get_in_edge_count();
            uint8_t out = t.second.get_out_edge_count();

            return (in == 0) && (out == 0);
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
            uint8_t in = t.second.get_in_edge_count();
            uint8_t out = t.second.get_out_edge_count();

            return (in == 1) && (out == 1);
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




    } //namespace filter
  } //namespace de_bruijn
} //namespace bliss




#endif // DE_BRUIJN_FILTER_HPP_
