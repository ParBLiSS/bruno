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

#include <stdint.h>
#include <tuple>        // tuple and utility functions

//#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "debruijn/debruijn_graph_node.hpp"

namespace bliss
{
  namespace debruijn
  {
    // data type for debruijn graph chain compaction
    // first Kmer is in edge, second Kmer is out edge.  
	// uint is distance from end nodes, with 
	//		MSB=true  + d>0 indicates the edge is pointing to terminus.  this is end state for non-termini after chain compaction.
	//		MSB=false + d>0 indicates the edge is pointing to chain node.  this is intermediate state of non-termini during chain compaction.
	// 		MSB=false + 0 indicate that node is terminus and is pointing to a branch.  (determined from graph branches later)
	//		MSB=true  + 0 indicate that node is terminus and is pointing to self (daedend).  (not determined from graph branches later, so this should be default)
	// used with Kmer key in a map, to mean in-Kmer-out bi edge.
    template <typename KMER>
    using simple_biedge = ::std::tuple<KMER, KMER, uint, uint>;
    //static_assert(sizeof(simple_biedge<::bliss::common::Kmer<31, ::bliss::common::DNA> >) == sizeof(::bliss::common::Kmer<31, ::bliss::common::DNA>) * 2 + 2 * sizeof(int), "size of simple biedge is not what's expected");

	constexpr uint dist_mask = std::numeric_limits<uint>::max() >> 1;
	constexpr uint term_mask = ~(::bliss::debruijn::dist_mask);  // MSB=true and dist=0
	constexpr uint branch_neighbor = static_cast<uint>(0);

	inline constexpr uint get_chain_dist(uint const & x) {
		return x & dist_mask;
	}
	inline constexpr bool is_chain_terminal(uint const & x) {
		return (x & dist_mask) == 0;
	}
	inline constexpr bool points_to_branch(uint const & x) {
		return x == static_cast<uint>(0);
	}
	inline constexpr bool points_to_self(uint const & x) {  // DEADEND
		return x == term_mask;
	}
	inline constexpr bool points_to_chain_node(uint const & x) {  // MSB not set and d > 0
		return (x ^ term_mask) > term_mask;  // flip MSB, then compare to term_mask
	}
	inline constexpr bool points_to_terminal(uint const & x) { // for non-terminal chain nodes.  MSB set.
		return x > term_mask;
	}
	inline constexpr bool points_to_chain_node_or_self(uint const & x) {  // MSB not set and d > 0
		return x <= term_mask;  // flip MSB, then compare to term_mask
	}
	inline constexpr bool points_to_terminal_or_self(uint const & x) { // for non-terminal and terminal chain nodes.  MSB set.
		return (x ^ term_mask) <= term_mask;
	}
	inline uint mark_as_point_to_terminal(uint const & x) {  // unmark as self pointing, and keep dist.
		return x | term_mask;
	}


    template <typename KMER>
	struct to_simple_biedge {

    	template <typename EdgeEncoding, typename COUNT, typename DUMMY>
    	std::pair<KMER, simple_biedge<KMER> >
    	operator()(::std::pair<KMER, ::bliss::debruijn::graph::compact_multi_biedge<EdgeEncoding, COUNT, DUMMY> > const & t) {

			simple_biedge<KMER> node(KMER(), KMER(), ::bliss::debruijn::term_mask, ::bliss::debruijn::term_mask);   // default node

			// get the in neighbor
			std::vector<KMER> neighbors;
			t.second.get_in_neighbors(t.first, neighbors);
			assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
			if (neighbors.size() == 1) {
				std::get<0>(node) = neighbors[0];
				std::get<2>(node) = static_cast<uint>(1);
			} // else if there is no neighbor, leave kmer as blank and distance as 0
//			else {
//				std::cout << " node without IN edge: " << t << std::endl;
//			}

			// get the out neighbor
			neighbors.clear();
			t.second.get_out_neighbors(t.first, neighbors);
			assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
			if (neighbors.size() == 1) {
				std::get<1>(node) = neighbors[0];
				std::get<3>(node) = static_cast<uint>(1);
			} // else if there is no neighbor, leave kmer as blank and distance as 0
//			else {
//				std::cout << " node without OUT edge: " << t << std::endl;
//			}

    		return std::make_pair(t.first, node);
    	}

    	template <typename EdgeEncoding, typename COUNT, typename DUMMY>
    	std::pair<KMER, simple_biedge<KMER> >
    	operator()(::std::pair<const KMER, ::bliss::debruijn::graph::compact_multi_biedge<EdgeEncoding, COUNT, DUMMY> > const & t) {

			simple_biedge<KMER> node(KMER(), KMER(), ::bliss::debruijn::term_mask, ::bliss::debruijn::term_mask);   // default node

			// get the in neighbor
			std::vector<KMER> neighbors;
			t.second.get_in_neighbors(t.first, neighbors);
			assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
			if (neighbors.size() == 1) {
				std::get<0>(node) = neighbors[0];
				std::get<2>(node) = static_cast<uint>(1);
			}

			// get the out neighbor
			neighbors.clear();
			t.second.get_out_neighbors(t.first, neighbors);
			assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
			if (neighbors.size() == 1) {
				std::get<1>(node) = neighbors[0];
				std::get<3>(node) = static_cast<uint>(1);
			}

    		return std::make_pair(t.first, node);
    	}


    };

    namespace chain
    {

      /**
       * @brief a node in a compacted chain.
       *
       * @tparam KMER 		first KMER field is the current k-mer, canonical
       * @tparam KMER		second KMER field is the lex smaller 5' terminus of the chain
       * @tparam int		int is distance to that terminus.  +/0 indicates both are on same strand.  - means otherwise..  no -0 needed because we are on canonical strand.
       *
       *      note that a vector of these are sufficient to represent a complete compacted chain.
       */
      template <typename KMER>
      using listranked_chain_node = ::std::tuple<KMER, KMER, int>;

      //static_assert(sizeof(listranked_chain_node<::bliss::common::Kmer<31, ::bliss::common::DNA> >) == sizeof(::bliss::common::Kmer<31, ::bliss::common::DNA>) * 2 + sizeof(int), "size of compacted chain node is not what's expected");


      /**
       * @brief summarized chain.
       *
       * @tparam KMER		first KMER field is the kmer 5' to the the 5' terminus (branch or deadend)
       * @tparam KMER 		second KMER field is the chain representative, also the 5' terminus
       * @tparam KMER 		3' terminus.
       * @tparam KMER		3' to the 3' terminus.  (branch or deadend)
       * @tparam uint		chain length
       * @tparam uint8		5' self reference?  0 means has branch neighbor.
       * @tparam uint8		3' self reference?  1 means pointing to self (same definition as the dists).
	   * 
       */
	  template <typename KMER, typename COUNT = uint8_t>
	  using summarized_chain = ::std::tuple<KMER, KMER, KMER, KMER, uint, COUNT, COUNT>;

		template <typename KMER, typename COUNT = uint8_t>
		inline ::bliss::debruijn::chain::summarized_chain<KMER, COUNT> 
		merge_summarized_chains(::bliss::debruijn::chain::summarized_chain<KMER, COUNT> const & x,
			::bliss::debruijn::chain::summarized_chain<KMER, COUNT> const & y ) {
			// middle parts have to be identical.
			::bliss::debruijn::chain::summarized_chain<KMER, COUNT> out = x;

				assert((std::get<1>(out) == std::get<1>(y)) && 
						(std::get<2>(out) == std::get<2>(y)) && 
					   (std::get<4>(out) == std::get<4>(y)) &&
					   "the chain terminal kmer and distance do not match between the two halves.");
				// max of in and out edge freq, and corresponding kmers.
				if (std::get<5>(out) < std::get<5>(y)) {
					std::get<0>(out) = std::get<0>(y);
					std::get<5>(out) = std::get<5>(y);
				}
				if (std::get<6>(out) < std::get<6>(y)) {
					std::get<3>(out) = std::get<3>(y);
					std::get<6>(out) = std::get<6>(y);
				}
			return out;
		}

	  /// make a summarized chain.  assumes that the termini are at left of first and right of second, and first has kmer corresponding to chain representative
	  template <typename KMER, typename COUNT = uint8_t>
	  inline ::bliss::debruijn::chain::summarized_chain<KMER, COUNT> make_summarized_chain(KMER const & ff, simple_biedge<KMER> const & fs,
	  	 KMER const & sf, simple_biedge<KMER> const & ss ) {
			return std::make_tuple(
			std::get<0>(fs),  // prev at 5'
			ff,		// 5', also chain rep
			sf,		// 3', all same strand as chain rep.
			std::get<1>(ss),  // next at 3'
			::bliss::debruijn::get_chain_dist(std::get<2>(fs)) + ::bliss::debruijn::get_chain_dist(std::get<3>(fs)),	 // DIST between the termini.
			::bliss::debruijn::points_to_self(std::get<2>(fs)) ? static_cast<COUNT>(0) : static_cast<COUNT>(1),  // 5' terminal
			::bliss::debruijn::points_to_self(std::get<3>(ss)) ? static_cast<COUNT>(0) : static_cast<COUNT>(1) // 3' terminal
			);
	  }

	  template <typename KMER, typename COUNT = uint8_t>
	  inline ::bliss::debruijn::chain::summarized_chain<KMER, COUNT> make_summarized_chain(KMER const & ff, simple_biedge<KMER> const & fs,
	  	 int ) {
			return std::make_tuple(
			std::get<0>(fs),  // prev at 5'
			ff,		// 5', also chain rep
			std::get<1>(fs),		// 3', all same strand as chain rep.
			KMER(),  // next at 3'
			::bliss::debruijn::get_chain_dist(std::get<2>(fs)) + ::bliss::debruijn::get_chain_dist(std::get<3>(fs)),	 // DIST between the termini.
			::bliss::debruijn::points_to_self(std::get<2>(fs)) ? static_cast<COUNT>(0) : static_cast<COUNT>(1),  // 5' terminal
			static_cast<COUNT>(0) // 3' terminal
			);
	  }

  	  template <typename KMER, typename COUNT = uint8_t>
	  inline ::bliss::debruijn::chain::summarized_chain<KMER, COUNT> make_summarized_chain(int,
	  	 KMER const & sf, simple_biedge<KMER> const & ss  ) {
			return std::make_tuple(
			KMER(),  // prev at 5'
			std::get<0>(ss),		// 5', also chain rep
			sf,		// 3', all same strand as chain rep.
			std::get<1>(ss),  // next at 3'
			::bliss::debruijn::get_chain_dist(std::get<2>(ss)) + ::bliss::debruijn::get_chain_dist(std::get<3>(ss)),	 // DIST between the termini.
			static_cast<COUNT>(0),  // 5' terminal
			::bliss::debruijn::points_to_self(std::get<3>(ss)) ? static_cast<COUNT>(0) : static_cast<COUNT>(1) // 3' terminal
			);
	  }

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


      // reverse complement the chain summary.
      template <typename KMER, typename COUNT=uint8_t>
      inline ::bliss::debruijn::chain::summarized_chain<KMER, COUNT>
      reverse_complement(::bliss::debruijn::chain::summarized_chain<KMER, COUNT> const & x) {
        return ::bliss::debruijn::chain::summarized_chain<KMER, COUNT>(
			std::get<3>(x).reverse_complement(),
            std::get<2>(x).reverse_complement(),
			std::get<1>(x).reverse_complement(),
            std::get<0>(x).reverse_complement(),
			std::get<4>(x),
			std::get<6>(x), std::get<5>(x));
      }

    } // namespace transform

  }/*namespace debruijn*/
}/*namespace bliss*/




#endif /* DEBRUIJN_CHAIN_NODE_HPP_ */
