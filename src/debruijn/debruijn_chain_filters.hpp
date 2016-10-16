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

        struct IsInternal {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) != 0) && (std::get<3>(t.second) != 0);
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
                // if same as right rev comp, mark as canonical
            	  return (t.first <= std::get<1>(t.second).reverse_complement());
              } else if (std::get<3>(t.second) == 0) {
                // if revcomp same as left, choose the left.
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
              if ((std::get<2>(t.second) == 0) && (std::get<3>(t.second) == 0)) {
                return true;
              } else if (std::get<2>(t.second) == 0) {
                // if same as right rev comp, mark as canonical
                return (t.first <= std::get<1>(t.second).reverse_complement());
              } else if (std::get<3>(t.second) == 0) {
                // if revcomp same as left, choose the left.
                return (t.first.reverse_complement() < std::get<0>(t.second));
              } else {
               return false;
              }
            }
        };


        //=========== below are filters used during and after chain compaction.

        /// Used during compaction.  Indicate that a node has it's edges pointing to terminal nodes.
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


        /// Used during compaction.  Indicate that a node has it's edges pointing to non-terminal nodes.
        struct PointsToInternalNode {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) > 0) || (std::get<3>(t.second) > 0);
            }

            template <typename Iter>
                        inline bool operator()(Iter it) const {
            	return operator()(*it);
            }
        };

        /// A node that MAY BE part of a cycle.  used during chain compaction for stopping criteria
        ///    also used after compaction to remove true cycle nodes.
        struct IsCycleNode {
            int max_distance;

            IsCycleNode(size_t const iter) : max_distance(0x1 << iter) {
				// commented out because it should be valid for list ranking to have an isolated node with distance of 1 to itself.   a
            	// so preventing cycle check to run should be based on other criteria (i.e. no unfinished nodes).
//            	assert(iter > 0);  // has to have at least 1 iteration.
            };

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (std::get<2>(t.second) == max_distance) &&
                  (std::get<3>(t.second) == max_distance);
            }

            template <typename Iter>
                        inline bool operator()(Iter it) const {
            	return operator()(*it);
            }
        };
// not used.
//        // points to internal node that are not cycles
//        struct IsPalindrome {
//
//            /// does not filter by group of results.
//            template <typename Iter>
//            inline bool operator()(Iter first, Iter last) const {  return true; }
//
//            // if both one or both of the edge are not pointing to a terminus, then this is a node that's in progress.
//            template <typename Kmer, typename Edge>
//            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
//              return (t.first == t.first.reverse_complement());
//            }
//        };

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
