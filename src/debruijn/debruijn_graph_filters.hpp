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
#include "debruijn/debruijn_graph_node.hpp"
#include "debruijn/kmer_traits.hpp"

namespace bliss {
  namespace debruijn {
    namespace filter {
      namespace graph {

        //============= predicates for global filter operations.

        struct IsKPalindrome {
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const { return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              // if DNA, and k is odd, then impossible to be palindrome.
              using KM = typename std::remove_cv<Kmer>::type;
              return ::bliss::common::kmer::kmer_traits<KM>::is_rc_palindrome(t.first);

              // replaced
              // if (std::is_same<typename Kmer::KmerAlphabet, ::bliss::common::DNA>::value &&
              //     (Kmer::size & 0x1) == 1 ) return false;

              // // else check.
              // return (t.first == t.first.reverse_complement());
            }
        };

        struct IsKplus1Palindrome {

            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {return true;}


            template <typename Kmer, typename Alphabet, typename CountType, typename DUMMY>
//				typename ::std::enable_if<!std::is_same<CountType, bool>::value, int>::type = 0>
            inline bool operator()(::std::pair<Kmer,
            		::bliss::debruijn::graph::compact_multi_biedge<Alphabet, CountType, DUMMY> > const & t) const {

              // if DNA, and k is even, then impossible to be k+1 palindrome.
              if (std::is_same<typename Kmer::KmerAlphabet, ::bliss::common::DNA>::value &&
                  (Kmer::size & 0x1) == 0 ) return false;

              using KM = typename std::remove_cv<Kmer>::type;
              KM k = t.first;

// TODO use get_edge_frequency(0..3) instead.
	      bool km1_low_is_palindrome = ::bliss::common::kmer::kmer_traits<KM>::is_kminus1_rc_palindrome_low(k);
	      for (unsigned char c = 0; c < 4; ++c) {
		if (t.second.get_out_edge_frequency(c) > 0) {
			if (::bliss::common::kmer::kmer_traits<KM>::is_k1_k1pair_rc_palindrome(k, c, km1_low_is_palindrome)) return true;
		}
	      }

	      bool km1_high_is_palindrome = ::bliss::common::kmer::kmer_traits<KM>::is_kminus1_rc_palindrome_high(k);
              // in
              for (unsigned char c = 0; c < 4; ++c) {
		if (t.second.get_in_edge_frequency(c) > 0) {
			if (::bliss::common::kmer::kmer_traits<KM>::is_k1_k1pair_rc_palindrome(c, k, km1_high_is_palindrome)) return true;
		}
	      }


              return false;
            }

/*            template <typename Kmer, typename Alphabet, typename CountType, typename DUMMY,
			typename ::std::enable_if<std::is_same<CountType, bool>::value, int>::type = 0>
            inline bool operator()(::std::pair<Kmer,
            		::bliss::debruijn::graph::compact_multi_biedge<Alphabet, CountType, DUMMY> > const & t) const {

              // if DNA, and k is even, then impossible to be k+1 palindrome.
              if (std::is_same<typename Kmer::KmerAlphabet, ::bliss::common::DNA>::value &&
                  (Kmer::size & 0x1) == 0 ) return false;

              using KM = typename std::remove_cv<Kmer>::type;
              KM k = t.first;
              KM krevcomp = k.reverse_complement();


              // TODO: remember to change this back to the short cutting version.

              std::vector<KM> neighbors;
//              bool res = false;

//              std::cout << " kmer: " << k << " neighbors: " << t.second << std::endl;

              // out kmer is same as revcmp of kmer.
              t.second.get_out_neighbors(k, neighbors);
              for (auto x : neighbors) {
//                std::cout << "k+1 bool palindrome out " << x << " k " << k << (x == krevcomp ? " 1" : " 0") << std::endl;
            	  if (x == krevcomp) {
//                  res = true;
            	    return true;
            	  }
              }

              // in
              neighbors.clear();
              t.second.get_in_neighbors(k, neighbors);
              for (auto x : neighbors) {
//                std::cout << "k+1 bool palindrome in  " << x << " k " << k << (x == krevcomp ? " 1" : " 0") << std::endl;
            	  if (x == krevcomp) {
//                  res = true;
            	    return true;
            	  }
              }

//              return res;
              return false;
            }
  */      };

        struct IsBranchPoint {
            IsKPalindrome isKPalindrome;
            IsKplus1Palindrome isKplus1Palindrome;

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return ((t.second.get_in_edge_count() >= 2) ||
                  (t.second.get_out_edge_count() >= 2))  ||
                  (isKPalindrome(t)) ||
                  (isKplus1Palindrome(t));
            }
        };

// used only for debugging, and need to be updated with palindrome checks?
        struct IsIsolated {
            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t) const {
              return (t.second.get_in_edge_count() == 0) && (t.second.get_out_edge_count() == 0);
            }
        };
//
//        struct IsTerminus {
//            /// does not filter by group of results.
//            template <typename Iter>
//            inline bool operator()(Iter first, Iter last)  const {  return true; }
//
//            template <typename Kmer, typename Edge>
//            inline bool operator()(::std::pair<Kmer, Edge> const & t) const  {
//              uint8_t in = t.second.get_in_edge_count();
//              uint8_t out = t.second.get_out_edge_count();
//
//              return ((in == 0) && (out == 1)) ||
//                  ((in == 1) && (out == 0));
//            }
//        };

        struct IsChainInternalNode {
            IsKPalindrome isKPalindrome;
            IsKplus1Palindrome isKplus1Palindrome;

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last) const  {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t)  const {
              return (t.second.get_in_edge_count() == 1) &&
                  (t.second.get_out_edge_count() == 1) &&
                  (!isKPalindrome(t)) && (!isKplus1Palindrome(t));
            }
        };


        struct IsChainNode {
            IsKPalindrome isKPalindrome;
            IsKplus1Palindrome isKplus1Palindrome;

            /// does not filter by group of results.
            template <typename Iter>
            inline bool operator()(Iter first, Iter last)  const {  return true; }

            template <typename Kmer, typename Edge>
            inline bool operator()(::std::pair<Kmer, Edge> const & t)  const {
              return ((t.second.get_in_edge_count() < 2) &&
            		      (t.second.get_out_edge_count() < 2)) &&
                  (!isKPalindrome(t)) && (!isKplus1Palindrome(t));
            }
        };

      } // namespace graph
    } //namespace filter
  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_GRAPH_FILTERS_HPP_
