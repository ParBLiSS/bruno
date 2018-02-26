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
 * debruijn_compacted_graph.hpp
 *
 *  Created on: Oct 2, 2016
 *      Author: tony pan
 */

#ifndef DEBRUIJN_TOPO_OPERATIONS_HPP_
#define DEBRUIJN_TOPO_OPERATIONS_HPP_

//#include "index/kmer_index.hpp"
#include "utils/function_traits.hpp"

#include "iterators/algorithm.hpp"


#include "debruijn/debruijn_chain_graph.hpp"
#include "debruijn/debruijn_graph.hpp"


namespace bliss
{

namespace debruijn
{

namespace graph
{

	/**
	 * @brief 	chained debruijn graph has its straight chains list ranked from the chain representative kmers.
	 * @detials	just the chains.  makes the implementation simpler.
	 * 			operations on graphs with chains should be implemented
	 * 			as separate functions that takes as params a DBG and a CHAIN_GRAPH
	 *
	 * 			note that if graph were modified, e.g. edge removed, new straight chains may form
	 * 			and we may need to recompact.
	 * 				the list ranking algorithm used requires all chain nodes to participate.
	 * 				applying the same algorithm to all nodes requires
	 * 					a. branch nodes changing to chain nodes
	 * 					b. changed branch nodes update neighbors to indicate they are not
	 * 						termini anymore
	 * 					c. all termini and newly converted chain nodes perform list ranking
	 * 					d. all former internal chain nodes perform query to get new termini and distances
	 * 				the internal chain nodes of a compacted chain do not get updated during recompaction,
	 * 				since no other nodes point to them.  so they need to use query to resolve this.
	 *
	 * 			 using just the ends of the chain to perform compaction.
	 * 			 	then update internals after all done.
	 * 			 this is so that recompaction does not need to propagate value along entire chain one at a time.
	 *
	 * 			 problem:  could encounter cycles again.
	 * 			 	need a different way to detect cycles:
	 * 			 		are the distances all doubling?  not necessarily.  look at a cycle with 3 segments, x ,y ,z.
	 * 			 		after 1 iteration, we have x+y, y+z, z+x as the lengths.  next iteration
	 * 			 		x+2y+z, x+y+2z, 2x+y+z, 3 iterations:
	 * 			 		2x+3y+3z, 3x+2y+3z, 3x+3y+2z.
	 *
	 * 			 		total distance between nodes are x+y+z, 2x+2y+2z, 4x+4y+4z, 8x+8y+8z, etc.
	 * 			 		if cycle, then total distance double
	 *
	 * 			 		if total distance double, then cycle.
	 *
	 * 			 		total distance not double, then have non-cycle
	 *
	 * 			TODO: above.
	 *
	 * 			also note that because we update to max distance jumped, the intermediate nodes MAY end up not being updated properly, hence the query and local update step.
	 *				to prove that this is necessary, we need to show that there exists a "driver" node that propagates the terminal node id and distance to the target node
	 *				between 0 and n-1 during some iteration i.
	 *				update to node x is handle by driver node y at x-2^(i-1) or n - x + 2^(i-1)
	 *				connection between driver node y and x is set up during iterator i-1.
	 *				driver node y is ranked wrt to node 0 during iteration ceil(log(x-2^(i-1))).
	 *				recursive relationship, with base case distances = 0 or 1.
	 *				another way to look at it, each iterator allows the first and last 2^i nodes to complete 1 side of the edges.
	 *				during next iteration, another 2^i nodes from the ends are completed.
	 *				to reach other end, need ceil(log(n)) iterations.
	 *
	 *			above is prove that we do not need a query phase for list ranking
	 *
	 *			for recompacting, doubling approach requires log(n) iterations in the length of chain, even if chain is already compacted
	 *				recompact using end nodes as representatives, then update, should be faster.
	 *
	 *			FINAL, we want to be able to recompact multiple times (iteratively),
	 *				and we want to filter by data on chain nodes, e.g. frequency or quality score
	 *				suggesting that 1. we should keep the original dbg around, and chain nodes explicit.
	 *					(condensing may be memory efficient, but is much harder to get at data...)
	 *				separate call to compress.
	 *				when compressing, also calculate chain stats.
	 *
	 *			also keeping original dbg around allows for cutting chains -
	 *				can rebuild specific chains.
	 *
	 *			also, data distribution should be more even
	 *			don't have to deal with reverse complement of entire contig when recompacting, or when compacting read fragments.
	 *
	 *			trade off - not compact - ANYWAY TO COMPACT THIS MORE?
	 *
	 *			possible to use read fragments? - probably not.
	 *				read fragments may overlap.  we'd have to index all kmers in the fragment to be able to link
	 *				unless minhash is used.
	 *
	 *			need to define return types for variables.
	 *
	 * @note:  not subclassing.  this is just the chain part.
	 */


	// to do proper filtering, need to
	//   1. compute chain statistics - frequency min, max, mean, stdev.
	//	 2. then we can perform filtering
	//   3. then we can perform recompaction (using chain representatives)
	//   4. and new cycle detection and removal (sum of distances of all nodes is power of 2 (probably not strictly enough)


	//========= graph topological filtering.

	// approach:  1. find chain termini for both ends, sort, then filter out non-deadends.
	//			  2. find chain representatives for deadends, then find termini for these chains
	// 	return chain representative?  or termini?  or adjacent branches?
	//    subsequent:  get adjacent branches
	//					get termini
	//					get all chain nodes		  

	/**
	 * @brief return a vector of chain terminals from deadends.
	 */
	template <typename ChainGraph >
	std::vector<typename ChainGraph::mutable_value_type> find_deadend_termini(ChainGraph const & chains) {
		// Question: what is the "next" for a deadend terminal node?
			// in chain graph, AAAAAA with dist 0.
			// may be confused with true AAAAAA branch?
			// to real branch:  <TAAAAA, <XXXXXX, AAAAAA, X, 0> >
			// to deadend:  <YYYYYY, <XXXXXX, AAAAAA, X, 0> >
			// and reverse comp becomes TTTTTT?
			// 0 is indicating pointing to end, and deadend.
			//  may be better to consider as distance to "branch"

		
		// either search for terminal with no next, then search again for terminals with matching 
		// or search for all terminal then sort by 


		// distributed sort by chain rep, then send one to prev - 

	}
	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype count_deadends(std::pair<dbg, chain_graph> const & x, predicate) {
		
	}
	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype erase_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	}

	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype find_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	}
	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype count_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	}
	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype erase_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	}



	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype find_bubbles(std::pair<dbg, chain_graph> const & x, predicate) {

	}
	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype count_bubbles(std::pair<dbg, chain_graph> const & x, predicate) {

	}
	/**
	 * @brief return a vector of chain representatives, each a separate deadend.
	 */
	template <typename somedatatype, >
	somedatatype erase_bubbles(std::pair<dbg, chain_graph> const & x, predicate) {

	}
	// /**
	//  * @brief return a vector of chain representatives, each a separate deadend.
	//  */
	// template <typename somedatatype, >
	// somedatatype recursively_erase_bubbles(std::pair<dbg, chain_graph> const & x, predicate) {

	// }




} // ns: graph

} // ns: debruijn

} // ns: bliss

#endif // DEBRUIJN_TOPO_OPERATIONS_HPP_
