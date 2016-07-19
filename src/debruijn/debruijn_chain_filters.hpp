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
 * debruijn_chain_filters.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_CHAIN_FILTERS_HPP_
#define DEBRUIJN_CHAIN_FILTERS_HPP_

#include <vector>

#include "debruijn/debruijn_chain_node.hpp"

namespace bliss {
  namespace debruijn {
    namespace filter {

      namespace chain {

        struct IsIsolated {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) == 0) && (std::get<3>(t.second) == 0);
            }
        };

        struct IsTerminus {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) == 0) ^ (std::get<3>(t.second) == 0);
            }
        };

        struct IsTerminusOrIsolated {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) == 0) || (std::get<3>(t.second) == 0);
            }
        };


        struct IsCanonicalTerminus {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              if (! ((std::get<2>(t.second) == 0) ^ (std::get<3>(t.second) == 0)) ) return false;

              if (std::get<2>(t.second) == 0) {
            	  return (t.first < std::get<1>(t.second).reverse_complement());
              } else if (std::get<3>(t.second) == 0) {
            	  return (t.first.reverse_complement() < std::get<0>(t.second));
              } else {
            	 return false;
              }
            }
        };

        struct IsCanonicalTerminusOrIsolated {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              if (! ((std::get<2>(t.second) == 0) || (std::get<3>(t.second) == 0)) ) return false;

              if ((std::get<2>(t.second) == 0) && (std::get<3>(t.second) == 0)) {
                return true;
              } else if (std::get<2>(t.second) == 0) {
                return (t.first < std::get<1>(t.second).reverse_complement());
              } else if (std::get<3>(t.second) == 0) {
                return (t.first.reverse_complement() < std::get<0>(t.second));
              } else {
               return false;
              }
            }
        };

        struct PointsToTermini {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both the node is isolated, (0, 0), both pointing to ends (-x,-y),
            // or is a terminus and pointing to other terminus (0, -y), then this is a node that's finished.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) <= 0) && (std::get<3>(t.second) <= 0);
            }
        };


        struct PointsToInternalNode {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) > 0) || (std::get<3>(t.second) > 0);
            }
        };

        // points to internal node that are not cycles
        struct IsCycleNode {
            int max_distance;

            IsCycleNode(size_t const iter) : max_distance(0x1 << iter) {};

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) == max_distance) &&
                  (std::get<3>(t.second) == max_distance);
            }
        };

        // points to internal node that are not cycles
        struct IsPalindrome {

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (t.first == t.first.reverse_complement());
            }
        };

        // points to internal node that are not cycles
        struct IsUncompactedNode {
            int max_distance;

            IsUncompactedNode(size_t const iter) : max_distance(0x1 << iter) {};

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return ((std::get<2>(t.second) > 0) || (std::get<3>(t.second) > 0)) && // at least 1 is positive
                  !((std::get<2>(t.second) == max_distance) && (std::get<3>(t.second) == max_distance));  // not both are at max_distance.
            }
        };


      }  // namespace dbg_chain

    } //namespace filter
  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_CHAIN_FILTERS_HPP_
