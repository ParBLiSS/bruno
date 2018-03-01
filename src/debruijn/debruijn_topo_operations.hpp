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
#include "debruijn/debruijn_chain_filters.hpp"

namespace bliss
{

namespace debruijn
{

namespace topology
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


	// approach:  1. find chain termini for both ends, sort, then filter out non-deadends.  <=
	//			  2. find chain representatives for deadends, then find termini for these chains.  the second find is more complex.
	// 	return chain representative?  or termini?  or adjacent branches?
	//    subsequent:  get adjacent branches
	//					get termini
	//					get all chain nodes?		  

	/**
	 * @brief return a vector of chain terminals from deadends.  isolated chain nodes are ignored as they are NOT deadends.
	 * 		
	 */
	template <typename ChainGraph >
	std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> >
	find_deadends(ChainGraph const & chains) {
		// Question: what is the "next" for a deadend terminal node?
			// in chain graph, AAAAAA with dist 0.
			// may be confused with true AAAAAA branch?  yes. (although in reality not too many of those)
			// to real branch:  <TAAAAA, <XXXXXX, AAAAAA, X, 0> >
			// to deadend:  <YYYYYY, <XXXXXX, AAAAAA, X, 0> >
			// and reverse comp becomes TTTTTT?
			// 0 is indicating pointing to end, and deadend.
			//  instead, dist should have type uint, with top bit indicating pointer to termini if dist is not 0, and pointing to self (no neighbor) if dist is 0.

		
		// either search for terminal with no next, then search again for terminals with matching - 
		//	this would require sorted structure, or hash table on a chainrep instead of node k-mer as key, for searching.
	    //   O(2M/P) for deadends,
		// or search for all terminal then distributed sort by chain rep, then send one to prev, then scan to filter out non-deadends
		//   O(2M/P) for getting termini, sort in O(2* (2M/P) log (2M/p) + distribute_cost), then scan O(2M/P).
		// search all termini then sort is SIMPLER.

		std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> > termini = chains.to_summarized_chains();

		// iterate over all termini to summarize deadend chains
		// then scan to locate deadends.  we must have a start and an end.  compact in place.  isolated are ignored.
		auto end = std::partition(termini.begin(), termini.end(),
			[](::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> const & x){
				return (std::get<5>(x) > 0) ^ (std::get<6>(x) > 0);  // at least one edge pointing to self (5. 6).  length can be 1.
		});
		termini.erase(end, termini.end());
		
		return termini;
	}
	
	template <typename ChainGraph >
	std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> >
	find_isolated_chains(ChainGraph const & chains) {
		// Question: what is the "next" for a deadend terminal node?
			// in chain graph, AAAAAA with dist 0.
			// may be confused with true AAAAAA branch?  yes. (although in reality not too many of those)
			// to real branch:  <TAAAAA, <XXXXXX, AAAAAA, X, 0> >
			// to deadend:  <YYYYYY, <XXXXXX, AAAAAA, X, 0> >
			// and reverse comp becomes TTTTTT?
			// 0 is indicating pointing to end, and deadend.
			//  instead, dist should have type uint, with top bit indicating pointer to termini if dist is not 0, and pointing to self (no neighbor) if dist is 0.

		
		// either search for terminal with no next, then search again for terminals with matching - 
		//	this would require sorted structure, or hash table on a chainrep instead of node k-mer as key, for searching.
	    //   O(2M/P) for deadends,
		// or search for all terminal then distributed sort by chain rep, then send one to prev, then scan to filter out non-deadends
		//   O(2M/P) for getting termini, sort in O(2* (2M/P) log (2M/p) + distribute_cost), then scan O(2M/P).
		// search all termini then sort is SIMPLER.

		std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> > termini = chains.to_summarized_chains();

		// iterate over all termini to summarize deadend chains
		// then scan to locate deadends.  we must have a start and an end.  compact in place.  isolated are ignored.
		auto end = std::partition(termini.begin(), termini.end(),
			[](::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> const & x){
				return (std::get<5>(x) > 0) && (std::get<6>(x) > 0);  // both edges are isolated.
		});
		termini.erase(end, termini.end());
		
		return termini;
	}


	
	template <typename ChainGraph >
	std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> >
	find_unit_chains(ChainGraph const & chains) {
		// Question: what is the "next" for a deadend terminal node?
			// in chain graph, AAAAAA with dist 0.
			// may be confused with true AAAAAA branch?  yes. (although in reality not too many of those)
			// to real branch:  <TAAAAA, <XXXXXX, AAAAAA, X, 0> >
			// to deadend:  <YYYYYY, <XXXXXX, AAAAAA, X, 0> >
			// and reverse comp becomes TTTTTT?
			// 0 is indicating pointing to end, and deadend.
			//  instead, dist should have type uint, with top bit indicating pointer to termini if dist is not 0, and pointing to self (no neighbor) if dist is 0.

		
		// either search for terminal with no next, then search again for terminals with matching - 
		//	this would require sorted structure, or hash table on a chainrep instead of node k-mer as key, for searching.
	    //   O(2M/P) for deadends,
		// or search for all terminal then distributed sort by chain rep, then send one to prev, then scan to filter out non-deadends
		//   O(2M/P) for getting termini, sort in O(2* (2M/P) log (2M/p) + distribute_cost), then scan O(2M/P).
		// search all termini then sort is SIMPLER.

		std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> > termini = chains.to_summarized_chains();

		// iterate over all termini to summarize deadend chains
		// then scan to locate deadends.  we must have a start and an end.  compact in place.  isolated are ignored.
		auto end = std::partition(termini.begin(), termini.end(),
			[](::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> const & x){
				return (std::get<4>(x) == 0) && (std::get<5>(x) == 0) && (std::get<6>(x) == 0);  // chains that sit between 2 branch nodes.
		});
		termini.erase(end, termini.end());
		
		return termini;
	}


	// /**
	//  * @brief 
	//  */
	// template <typename somedatatype, >
	// somedatatype count_deadends(std::pair<dbg, chain_graph> const & x, predicate) {
		
	// }
	// /**
	//  * @brief 
	//  */
	// template <typename somedatatype, >
	// somedatatype erase_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	// }

	// /**
	//  * @brief 
	//  */
	// template <typename somedatatype, >
	// somedatatype find_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	// }
	// /**
	//  * @brief 
	//  */
	// template <typename somedatatype, >
	// somedatatype count_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	// }
	// /**
	//  * @brief 
	//  */
	// template <typename somedatatype, >
	// somedatatype erase_deadends(std::pair<dbg, chain_graph> const & x, predicate) {

	// }


	/**
	 * @brief return a vector of chain terminals from bubbles. 
	 * 		note that bubble lengths will be at least k for alternate paths to form bubbles
	 */
	template <typename ChainGraph >
	std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> >
	find_bubbles(ChainGraph const & chains, ::mxx::comm const & comm) {
		// need to bring termini of the same chain together first in order to get the branches.
		// do this with sort by chain rep, then filter, then sort by the branches, and filter again.

		// don't believe that this can be done in 1 sort pass because left and right terminal each have incomplete info.
		// TODO: one possible alternative is when compacting, send the terminal's neighbors so that every node knows the branch.
		//   	this can't work, however, since the middle chain nodes would not have chain termini info for sorting/ordering later.

		// search for all terminal then distributed sort by chain rep, then send one to prev, then scan to filter out deadends and isolated.
		//   O(2M/P) for getting termini, sort in O(2* (2M/P) log (2M/p) + distribute_cost), then scan O(2M/P).

		std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> > termini = chains.to_summarized_chains();

		// iterate over all termini to summarize deadend chains
		// then scan to locate non-deadends.
		auto end = std::partition(termini.begin(), termini.end(),
			[](::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> const & x){
				return (std::get<4>(x) > 0) && (std::get<5>(x) == 0) && (std::get<6>(x) == 0);  
		});
		termini.erase(end, termini.end());

		// now that deadends and isolated chains have been removed, now sort by in edge of chain rep, followed by out edge.
		// evenly redistribute.
		mxx::distribute_inplace(termini, comm);
		
		// participate only if there is at least one entry.
		bool has_data = (termini.size() > 0);
		bool all_has_data = ::mxx::all_of(has_data, comm);
		mxx::comm subcomm = all_has_data ? comm.copy() : comm.split(has_data);

		// NOTE: two chains going out of 1 branch, but arrive at a distal branch at opposite strands, are NOT considered bubbles.

		if (has_data) {
			// CANONICALIZE summarized chains by the branch points (5' vs rc(3'))
			for (auto it = termini.begin(); it != termini.end(); ++it) {
				auto x = *it;
				if (std::get<0>(x) > std::get<3>(x).reverse_complement()) {
					*it = ::bliss::debruijn::transform::reverse_complement(x);
				}
			}

			// now sort by 5' branch, then by 3' branch
			mxx::sort(termini.begin(), termini.end(),
				[](::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> const & lhs,
					::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> const & rhs){
					// get chain representatives, then compare.
					return (std::get<0>(lhs) < std::get<0>(rhs)) || 
							((std::get<0>(lhs) == std::get<0>(rhs)) && (std::get<3>(lhs) < std::get<3>(rhs)));
			}, subcomm);


			// now that sorting is done, shift so all entries for the same branch are together.   note that there are at most 4 chains per source branch.
			std::pair<typename ChainGraph::kmer_type, int> splitter = 
				std::make_pair(std::get<0>(termini.back()), subcomm.size() - subcomm.rank());  // use 5' branch as splitter.
			// now do exclusive scan with max for key and min for rank.
			splitter = mxx::exscan(splitter, [](std::pair<typename ChainGraph::kmer_type, int> const & x, std::pair<typename ChainGraph::kmer_type, int> const & y){
				return (x.first == y.first) ? ((x.second > y.second) ? x : y) : 
					((x.first > y.first) ? x : y);  // if kmer same, return smaller comm rank.  else, return larger kmer.
			}, subcomm);
			splitter.second = subcomm.size() - splitter.second;  // convert back to rank
			// now set up alltoallv
			std::vector<size_t> send_counts(subcomm.size(), 0);
			// calculate the number to send.
			size_t i = 0;
			while (std::get<0>(termini[i]) == splitter.first) ++i;
			send_counts[splitter.second] = i;  // all else is 0.	
			// now move data using alltoallv
			std::vector<::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> > extras = 
				mxx::all2allv(termini.data(), send_counts, subcomm);
			termini.insert(termini.end(), extras.begin(), extras.end());  // order should be maintained.

			// then scan and keep only bubble entries, regardless of chain lengths.
			// start with i-th entry.
			size_t insert_at = 0;
			i = send_counts[splitter.second];  // these have been sent away.
			::bliss::debruijn::chain::summarized_chain<typename ChainGraph::kmer_type> first = termini[i];
			size_t first_i = i;
			++i;
			for (; i < termini.size(); ++i) {
				if ((std::get<0>(first) == std::get<0>(termini[i])) &&
					(std::get<3>(first) == std::get<3>(termini[i]))) {
					continue;
				}  // search until different.

				if ((i - first_i) > 1) {  // has at least 1 pair.
					// copy over.
					for (; first_i < i; ++first_i, ++insert_at) {
						termini[insert_at] = termini[first_i];
					}
				}

				first = termini[i];
				first_i = i;
			}
			// handle the last part.
			if ((i - first_i) > 1) {  // has at least 1 pair.
				// copy over.
				for (; first_i < i; ++first_i, ++insert_at) {
					termini[insert_at] = termini[first_i];
				}
			}
			// clear out the remainder.
			termini.erase(termini.begin() + insert_at, termini.end());
		}

		return termini;
	}
	// /**
	//  * @brief return a vector of chain representatives, each a separate deadend.
	//  */
	// template <typename somedatatype, >
	// somedatatype count_bubbles(std::pair<dbg, chain_graph> const & x, predicate) {

	// }
	// /**
	//  * @brief return a vector of chain representatives, each a separate deadend.
	//  */
	// template <typename somedatatype, >
	// somedatatype erase_bubbles(std::pair<dbg, chain_graph> const & x, predicate) {

	// }
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
