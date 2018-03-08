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
 * debruijn_chain_graph.hpp
 *
 *  Created on: Oct 2, 2016
 *      Author: tony pan
 */

#ifndef DEBRUIJN_CHAIN_GRAPH_HPP_
#define DEBRUIJN_CHAIN_GRAPH_HPP_

//#include "index/kmer_index.hpp"
#include "utils/function_traits.hpp"

#include "iterators/algorithm.hpp"


#include "debruijn/debruijn_chain_node.hpp"
//#include "debruijn/debruijn_chain_filters.hpp"
#include "debruijn/debruijn_chain_operations.hpp"
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
	template <typename KmerType, template <typename> class DistHash = ::bliss::index::kmer::DistHashMurmur>
	class debruijn_chain_graph {

	public:
		template <typename KMER>
		using map_params_template = ::bliss::debruijn::CanonicalDeBruijnHashMapParams<KMER, DistHash>;
		using map_params_type =
				map_params_template<KmerType>;

		using map_type = ::dsc::densehash_map<KmerType,
				::bliss::debruijn::simple_biedge<KmerType>,
				map_params_template,
				 ::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

		using kmer_type = typename map_type::key_type;
		using edge_type = typename map_type::mapped_type;
		using value_type = typename map_type::value_type;
		using mutable_value_type = ::std::pair<kmer_type, edge_type>;
		using iterator = typename map_type::iterator;

	protected:
		// full debruijn graph is stored elsewhere.
		// function that requires info from both should have both as parameters

		// chain map - populated after calling listrank.
		map_type map;
		map_type cycle_nodes;  // for storing the cycle nodes.
		map_type isolated;  // for storing the cycle nodes.

		// communicator
		mxx::comm comm;

		bool listranked;
		bool modified;

	public:
		debruijn_chain_graph(mxx::comm const & _comm) :
			map(_comm), cycle_nodes(_comm), isolated(_comm), comm(_comm.copy()), listranked(false), modified(false) {

			map.clear();
		}

		/// can only construct from debruijn graph
		template <typename C>
		debruijn_chain_graph(
				::bliss::debruijn::graph::debruijn_graph<KmerType, C> const & dbg,
				mxx::comm const & _comm)
			: debruijn_chain_graph(_comm) {

			// build the chain graph now.
			make_chain_map(dbg);

			// and compact.
			list_rank();
		}

		virtual ~debruijn_chain_graph() {};

		/// get internal chain map
		map_type & get_map() { return map; }
		map_type const & get_map() const { return map; }
		map_type & get_cycle_map() { return cycle_nodes; }
		map_type const & get_cycle_map() const { return cycle_nodes; }
		map_type & get_isolated_map() { return isolated; }
		map_type const & get_isolated_map() const { return isolated; }

		// INSERT OPERATIONS!
		/// insert chain nodes.  version for chain map native edge type.
		size_t insert(std::vector<mutable_value_type> & nodes,
				bool is_local = false) {

			// do not local reserve since repeats may be high and the map is a reduction map.
			//map.local_reserve(node.size());
			// take the realloc hit.

			bool all_local = mxx::all_of(is_local, comm);
			modified = true;

			if (all_local) {
				return map.local_insert(nodes.begin(), nodes.end());
			} else {
				return map.insert(nodes);  // includes distribution.
			}
		}

		/**
		 * @brief	insert the selected node into index.  predicate returns true means the element should be inserted.
		 * @details	example of predicate:  using another index along with a subpredicate, evaluate 1. if the element is in that index, and 2. is the subpredicate satisfied.
		 * 					we can insert e.g. just branch nodes that are in the graph, into the node.
		 */
		template <typename Predicate>
		size_t insert(std::vector<mutable_value_type> & nodes,
				Predicate const & pred,
				bool is_local = false) {

			// doing it this way to ensure that we get the intersection and not a complement
			bool all_local = mxx::all_of(is_local, comm);
			modified = true;

			if (all_local) {
				return map.local_insert(nodes.begin(), nodes.end(), pred);
			} else {
				return map.insert(nodes, pred);  // includes distribution.
			}
		}

		/// insert chain nodes.  version for dbg native edge type.
		template <typename C>
		size_t insert(
				std::vector<typename ::bliss::debruijn::graph::debruijn_graph<KmerType, C>::mutable_value_type> const & nodes,
				bool is_local = false) {

			// do not local reserve since repeats may be high and the map is a reduction map.
			//map.local_reserve(node.size());
			// take the realloc hit.

			// convert first.
			::bliss::debruijn::to_simple_biedge<KmerType> to_biedge;
			::bliss::debruijn::filter::graph::IsChainNode is_chain;

			bool all_local = mxx::all_of(is_local, comm);
			size_t count = 0;
			if (all_local) {
				for (auto t : nodes) {
					if (is_chain(t)) {
						map.get_local_container.insert(to_biedge(t));
						++count;
					}
				}
			} else {

				// communication is involved, so need to copy.
				std::vector<mutable_value_type> temp;
				::fsc::back_emplace_iterator<std::vector<mutable_value_type> > emplacer(temp);
				bliss::iterator::algorithm::transform_if(
						nodes.begin(), nodes.end(), emplacer,
						to_biedge, is_chain);

				count = map.insert(temp);  // includes distribution.
			}
			modified = true;

			return count;
		}

		/**
		 * @brief	insert the selected node into index.  predicate returns true means the element should be inserted.
		 * @details	example of predicate:  using another index along with a subpredicate, evaluate 1. if the element is in that index, and 2. is the subpredicate satisfied.
		 * 					we can insert e.g. just branch nodes that are in the graph, into the node.
		 */
		template <typename C, typename Predicate>
		size_t insert(
				std::vector<typename ::bliss::debruijn::graph::debruijn_graph<KmerType, C>::mutable_value_type> const & nodes,
				Predicate const & pred,
				bool is_local = false) {

			::bliss::debruijn::to_simple_biedge<KmerType> to_biedge;
			::bliss::debruijn::filter::graph::IsChainNode is_chain;

			std::vector<mutable_value_type> temp;
			::fsc::back_emplace_iterator<std::vector<mutable_value_type> > emplacer(temp);
			bliss::iterator::algorithm::transform_if(
					nodes.begin(), nodes.end(), emplacer,
					to_biedge, is_chain);

			// now hand off.
			return this->insert(temp, pred, is_local);
			// predicate needs to operate on range, so need to transform in batch.
		}



		// use the a2a collective find since this is not a multimap.
		std::vector<mutable_value_type> find(std::vector<KmerType> &query) const {
			return map.find(query);
		}

		template <typename Predicate>
		std::vector<mutable_value_type> find_if(std::vector<KmerType> &query, Predicate const &pred) const {
			return map.find(query, false, pred);
		}

		template <typename Predicate>
		std::vector<mutable_value_type> find_if(Predicate const &pred) const {
			return map.find(pred);
		}

		template <typename Predicate>
		std::vector<iterator> find_if_iterators(Predicate const &pred) const {
			return map.find_iterators(pred);
		}

		//==== SET of find function that performs transform of output during query processing.
		// use the a2a collective find since this is not a multimap.
		template <typename Transform>
		std::vector<typename ::bliss::functional::function_traits<Transform, mutable_value_type>::return_type >
		find_and_transform(std::vector<KmerType> &query, Transform const & trans) const {
			return map.find_transform(query, trans);
		}

		template <typename Predicate, typename Transform>
		std::vector<typename ::bliss::functional::function_traits<Transform, mutable_value_type>::return_type >
		find_if_and_transform(std::vector<KmerType> &query, Predicate const &pred, Transform const & trans) const {
			return map.find_transform(query, false, pred, trans);
		}

		template <typename Predicate, typename Transform>
		std::vector<typename ::bliss::functional::function_traits<Transform, mutable_value_type>::return_type >
		find_if_and_transform(Predicate const &pred, Transform const & trans) const {
			return map.find_transform(pred, trans);
		}


		std::vector< std::pair<KmerType, size_t> > count(std::vector<KmerType> &query) const {
			return map.count(query);
		}

		template <typename Predicate>
		std::vector<std::pair<KmerType, size_t> > count_if(std::vector<KmerType> &query, Predicate const &pred) const {
			return map.count(query, false, pred);
		}

		template <typename Predicate>
		std::vector<std::pair<KmerType, size_t> > count_if(Predicate const &pred) const {
			return map.count(pred);
		}


		void erase(std::vector<KmerType> &query) {
			map.erase(query);
		}

		template <typename Predicate>
		void erase_if(std::vector<KmerType> &query, Predicate const &pred) {
			map.erase(query, false, pred);
		}

		template <typename Predicate>
		void erase_if(Predicate const &pred) {
			map.erase(pred);
		}

		template <typename update_type, typename update_functor_type>
		size_t update(std::vector<update_type> & val, bool const & sorted,
				update_functor_type const & updater) {
			modified = true;
			return map.update(val, sorted, updater);
		}

	   typename map_type::const_iterator cbegin() const
	   {
		 return map.get_local_container().cbegin();
	   }

	   typename map_type::const_iterator cend() const {
		 return map.get_local_container().cend();
	   }

	   void clear() {
		   map.clear();
		   listranked = false;
		   modified = false;
	   }

	   size_t size() const {
		 return map.size();
	   }

	   size_t local_size() const {
		 return map.local_size();
	   }


	   //======  convenience functions

	   /// return local branch nodes
	   std::vector<mutable_value_type> get_internal_nodes() const {
		   return this->find_if(::bliss::debruijn::filter::chain::IsInternal());
	   }
	   /// return local branch node's kmer
	   std::vector<KmerType> get_internal_node_kmers() const {
		   return this->find_if_and_transform(::bliss::debruijn::filter::chain::IsInternal(),
				   [](mutable_value_type const & x){
			   return x.first;
		   });
	   }

	   /// return local branch nodes
	   std::vector<mutable_value_type> get_terminal_nodes() const {
		   return this->find_if(::bliss::debruijn::filter::chain::IsTerminusOrIsolated());
	   }
	   /// return local branch node's kmer
	   std::vector<KmerType> get_terminal_node_kmers() const {
		   return this->find_if_and_transform(::bliss::debruijn::filter::chain::IsTerminusOrIsolated(),
				   [](mutable_value_type const & x){
			   return x.first;
		   });
	   }

	   /// return local branch nodes
	   std::vector<mutable_value_type> get_chain_representatives() const {
		   return this->find_if(::bliss::debruijn::filter::chain::IsCanonicalTerminusOrIsolated());
	   }
	   /// return local branch node's kmer
	   std::vector<KmerType> get_chain_representative_kmers() const {
		   return this->find_if_and_transform(::bliss::debruijn::filter::chain::IsCanonicalTerminusOrIsolated(),
				   [](mutable_value_type const & x){
			   return x.first;
		   });
	   }

	   std::vector<mutable_value_type> get_cycle_nodes() const {
		   assert(this->listranked);

		   std::vector<mutable_value_type> results;
		   cycle_nodes.to_vector(results);
		   return results;
	   }
	   std::vector<KmerType> get_cycle_node_kmers() const {
		   assert(this->listranked);

		   std::vector<KmerType> results;
		   ::fsc::back_emplace_iterator<std::vector<KmerType> > emplacer(results);

		   std::transform(cycle_nodes.get_local_container().cbegin(),
				   cycle_nodes.get_local_container().cend(), emplacer,
				   [](value_type const & x){
			   return x.first;
		   });
		   return results;
	   }

	   std::vector<mutable_value_type> get_isolated() const {
		   assert(this->listranked);

		   std::vector<mutable_value_type> results;
		   isolated.to_vector(results);
		   return results;
	   }
	   std::vector<KmerType> get_isolated_kmers() const {
		   assert(this->listranked);

		   std::vector<KmerType> results;
		   ::fsc::back_emplace_iterator<std::vector<KmerType> > emplacer(results);

		   std::transform(isolated.get_local_container().cbegin(),
				   isolated.get_local_container().cend(), emplacer,
				   [](value_type const & x){
			   return x.first;
		   });
		   return results;
	   }



		template <typename C>
		void setup_chain_termini(::bliss::debruijn::graph::debruijn_graph<KmerType, C, DistHash> const & idx) {
			// allocate input
			//==  BATCHED version
			std::vector<std::pair<KmerType, bliss::debruijn::operation::chain::terminus_update_md<KmerType> > > all_neighbors;
			
			std::vector<KmerType> neighbors;
			neighbors.reserve(6); // ACGTN.

			// estimate the largest amount of memory to use.
			unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);

			// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
			size_t step = (free_mem / (8 * sizeof(std::pair<KmerType, bliss::debruijn::operation::chain::terminus_update_md<KmerType> > )));  // number of elements that can be held in freemem
			step = std::max(std::min(step, idx.local_size()), static_cast<size_t>(1));

			if (comm.rank() == 0) std::cout << "estimate num chain terminal updates=" << step << ", value_type size=" <<
					sizeof(std::pair<KmerType, bliss::debruijn::operation::chain::terminus_update_md<KmerType> > ) << " bytes" << std::endl;


			// 8 possible edges, so split the thing into 8 parts to bound the memory usage to about N.  lower in practice
			//size_t step = idx.size() / this->comm.size() / 8;
			//if (step == 0) step = this->comm.size();
			// size_t step = 1000000;
			all_neighbors.reserve(step);   // do in steps of 1000000
			size_t nsteps = (idx.local_size() + step - 1) / step;
			nsteps = mxx::allreduce(nsteps, [](size_t const & x, size_t const & y){
				return std::max(x, y);
			}, comm);

			auto iter = idx.cbegin();
			auto end = idx.cend();
			::bliss::debruijn::filter::graph::IsBranchPoint is_branch;
			::bliss::debruijn::operation::chain::terminus_update<KmerType> updater;
			size_t count = 0;

			for (size_t s = 0; s < nsteps; ++s) {

				all_neighbors.clear();

				// TODO: CONVERT TO USING filter_transform_iterator?

				// extract neighbors of branches
				for (size_t i = 0; (i < step) && (iter != end); ++i, ++iter) {
					// compute the chain rep
					if (!is_branch(*iter)) continue;
					auto t = *iter;

					neighbors.clear();
					t.second.get_out_neighbors(t.first, neighbors);
					for (auto n : neighbors) {  // OUT neighbor is used as IN edge
						// insert as is.  let lex_less handle flipping it.
						all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::debruijn::operation::IN));
						
						// NOT NEEDED - n is same as t.first, which is not a chain node so should not need to be updated.  // HANDLE k1 PALINDROME
						// if (n.reverse_complement() == t.first) 
						// 	all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::debruijn::operation::OUT));
					}

					neighbors.clear();
					t.second.get_in_neighbors(t.first, neighbors);
					for (auto n : neighbors) {  // IN neighbor is used as OUT edge
						// insert as is.  let lex_less handle flipping it.
						all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::debruijn::operation::OUT));

						// NOT NEEDED - n is same as t.first, which is not a chain node so should not need to be updated. // HANDLE k1 PALINDROME
						// if (n.reverse_complement() == t.first) 
						// 	all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::debruijn::operation::IN));
					}
				}

				// update the internal nodes.  only ones in chainmap are affected (no insertion)
				count += map.update(all_neighbors, false, updater );  // collective comm
			}

//				// ======= VERSION THAT IS NOT BATCHED.
//				//=== find branching nodes. local computation.
//				BL_BENCH_START(chain);
//				auto nodes = idx.find_if_iterators(::bliss::debruijn::filter::graph::IsBranchPoint());
//				BL_BENCH_COLLECTIVE_END(chain, "get_branches", nodes.size(), comm);
//
//				//=== mark neighbors of branch points.
//				BL_BENCH_START(chain);
//				std::vector<std::pair<KmerType, bliss::debruijn::operation::chain::terminus_update_md<KmerType> > > all_neighbors;
//				all_neighbors.reserve(nodes.size() * 2);  // branch has at least 2 edges
//
//				for (auto t : nodes) {
//					neighbors.clear();
//					(*t).second.get_out_neighbors((*t).first, neighbors);
//					for (auto n : neighbors) {
//						// insert as is.  let lex_less handle flipping it.
//						all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>((*t).first, bliss::debruijn::operation::IN));
//					}
//
//					neighbors.clear();
//					(*t).second.get_in_neighbors((*t).first, neighbors);
//					for (auto n : neighbors) {
//						// insert as is.  let lex_less handle flipping it.
//						all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>((*t).first, bliss::debruijn::operation::OUT));
//					}
//				}
//				BL_BENCH_COLLECTIVE_END(chain, "branch_neighbors", all_neighbors.size(), comm);
//
//
//				// mark chain termini that are adjacent to branch points, so we can mark them in the chainmap.
//				BL_BENCH_START(chain);
//				::bliss::debruijn::operation::chain::terminus_update<KmerType> updater;
//				size_t count = map.update(all_neighbors, false, updater );
//				BL_BENCH_COLLECTIVE_END(chain, "update_termini", count, comm);


		}



	   /**
	    * @brief extract chains from a debruijn graph
	    * @note	 existing chain map is cleared.
	    */
		template <typename C>
		void extract_chains(::bliss::debruijn::graph::debruijn_graph<KmerType, C, DistHash> const & idx) {

			size_t total_idx_size = idx.size();
			if (total_idx_size == 0) return;

			BL_BENCH_INIT(chain);

	          if (comm.rank() == 0) {
	        	  std::cout << "SIZES simple biedge size: " << sizeof(::bliss::debruijn::simple_biedge<KmerType>) <<
	        			  " kmer size " << sizeof(KmerType) <<
						  " node size " << sizeof(std::pair<KmerType, ::bliss::debruijn::simple_biedge<KmerType> >) << std::endl;
	          }


			if (comm.rank() == 0) printf("MAKE CHAINMAP\n");
			{
		//		BL_BENCH_START(chain);
		//		// find chain nodes
		//		auto chain_nodes = idx.find_if(::bliss::debruijn::filter::graph::IsChainNode());  // not isolated, not branching, not palindrome
		//		BL_BENCH_COLLECTIVE_END(chain, "get_chains", chain_nodes.size(), comm);

				// insert into local container inside chainmap.
				BL_BENCH_START(chain);
				::bliss::debruijn::filter::graph::IsChainNode is_chain;
				auto iter = idx.cbegin();
				auto end = idx.cend();

				::bliss::debruijn::to_simple_biedge<KmerType> to_biedge;
				// TODO: convert to using filter_transform_iterator.

				// initialize the chain map.  note that key kmers are canonical, same as in DBG.
				// note also that edge k-mers have same orientation as key kmers, and may not be canonical.
		//		for (auto t : chain_nodes) {
				for (; iter != end; ++iter) {
					if (!is_chain(*iter)) continue;

					// use to_simple_biedge to convert from compact_multi_biedge to simple_biedge.

					// chainmap uses the same distribution hash function and transform, so can insert locally.
					map.get_local_container().insert(to_biedge(*iter));

				}
				BL_BENCH_COLLECTIVE_END(chain, "insert in chainmap", map.local_size(), comm);

				//========= report.
				//      auto result = chainmap.find(::bliss::debruijn::filter::chain::IsTerminus());
				//      auto result2 = chainmap.find(::bliss::debruijn::filter::chain::IsIsolated());
				//      printf("chain map contains %lu chained termini and  %lu isolated\n", result.size(), result2.size());

			}


			// communicating in this phase.
			if (comm.rank() == 0) printf("MARK TERMINI NEXT TO BRANCHES\n");

			BL_BENCH_START(chain);
			setup_chain_termini(idx);
			BL_BENCH_COLLECTIVE_END(chain, "update_termini", map.local_size(), comm);

			BL_BENCH_REPORT_MPI_NAMED(chain, "chain_map", comm);

			modified = true;
		}


		size_t separate_cycles(size_t iterations) {
			::bliss::debruijn::filter::chain::IsCycleNode is_cycle(iterations);

			if (comm.rank() == 0) printf("REMOVE CYCLES\n");

			BL_BENCH_INIT(cycle);

			BL_BENCH_START(cycle);
			auto res = map.find(is_cycle);
			BL_BENCH_COLLECTIVE_END(cycle, "find cycle nodes", res.size(), comm);

			BL_BENCH_START(cycle);
			cycle_nodes.insert(res);
			BL_BENCH_COLLECTIVE_END(cycle, "copy to cycle nodes map", cycle_nodes.size(), comm);


			// remove cycle nodes from map.
			BL_BENCH_START(cycle);
			size_t cycle_node_count = map.erase(is_cycle);
			assert(cycle_node_count == res.size());
			BL_BENCH_COLLECTIVE_END(cycle, "erase cycle", cycle_node_count, comm);

			//printf("REMOVED %ld cycle nodes\n", cycle_node_count);

			cycle_node_count = mxx::allreduce(cycle_node_count, comm);
			if (comm.rank() == 0) printf("REMOVED %ld cycle nodes\n", cycle_node_count);

			BL_BENCH_REPORT_MPI_NAMED(cycle, "rem_cycle", comm);

			return cycle_node_count;
		}


		size_t separate_cycles() {
			::bliss::debruijn::filter::chain::IsInternal is_internal;

			if (comm.rank() == 0) printf("REMOVE CYCLES\n");

			BL_BENCH_INIT(cycle);

			BL_BENCH_START(cycle);
			auto res = map.find(is_internal);
			BL_BENCH_COLLECTIVE_END(cycle, "find cycle nodes", res.size(), comm);

			BL_BENCH_START(cycle);
			cycle_nodes.insert(res);
			BL_BENCH_COLLECTIVE_END(cycle, "copy to cycle nodes map", cycle_nodes.size(), comm);


			// remove cycle nodes from map.
			BL_BENCH_START(cycle);
			size_t cycle_node_count = map.erase(is_internal);
			assert(cycle_node_count == res.size());
			BL_BENCH_COLLECTIVE_END(cycle, "erase cycle", cycle_node_count, comm);

			//printf("REMOVED %ld cycle nodes\n", cycle_node_count);

			cycle_node_count = mxx::allreduce(cycle_node_count, comm);
			if (comm.rank() == 0) printf("REMOVED %ld cycle nodes\n", cycle_node_count);

			BL_BENCH_REPORT_MPI_NAMED(cycle, "rem_cycle", comm);

			return cycle_node_count;
		}

		size_t separate_isolated() {
			::bliss::debruijn::filter::chain::IsIsolated is_isolated;

			if (comm.rank() == 0) printf("REMOVE ISOLATED\n");

			BL_BENCH_INIT(isolated);

			BL_BENCH_START(isolated);
			auto res = map.find(is_isolated);
			BL_BENCH_COLLECTIVE_END(isolated, "find isolated nodes", res.size(), comm);

			BL_BENCH_START(isolated);
			isolated.insert(res);
			BL_BENCH_COLLECTIVE_END(isolated, "copy to isolated nodes map", isolated.size(), comm);


			// remove cycle nodes from map.
			BL_BENCH_START(isolated);
			size_t isolated_node_count = map.erase(is_isolated);
			assert(isolated_node_count == res.size());
			BL_BENCH_COLLECTIVE_END(isolated, "erase isolated", isolated_node_count, comm);

			//printf("REMOVED %ld isolated nodes\n", isolated_node_count);

			isolated_node_count = mxx::allreduce(isolated_node_count, comm);
			if (comm.rank() == 0) printf("REMOVED %ld isolated nodes\n", isolated_node_count);

			BL_BENCH_REPORT_MPI_NAMED(isolated, "rem_isolated", comm);

			return isolated_node_count;
		}


	   /**
	    * @brief 	processes the chains so that each points to the 2 ends of the chain, along with distances.
	    * @details	uses parallel list ranking algorithm internally.
	    * 			any cycles found (by checking that all remaining active nodes have left and right distance
	    * 			that are 2^iter.
	    *
	    * @note		this is the version used for initial list_ranking.
	    * 			this function would be inefficient for recompacting after some change.
		*           TODO: finding unfinished can be merged with the actual loop, and cycle count can be accumulated without using a separate loop.
	    */
	   size_t list_rank() {

			size_t iterations = 0;
			size_t cycle_node_count = 0;

			// optimizations:
			//		1. only send to non-end edges.  reduces update comm volume for termini.  this is probably normal
			//		2. only send if we are active (i.e. chain ranked).   reduces number of active chains.
			//		3. only send to the unfinished side(s), or to the finished side if the total distance is expected to be max (2^(iter + 1), or entire length of chain.)
			//		4. use vector of iterators to local chain map.  we are not changing the local chain map, so nothing is invalidated.
			//			then we can operate on the unfinished with partition, etc, also should save memory.
			// to prove:
			//		0. distance between left and right doubles each iteration. (middle nodes)
			//		1. all nodes between 2^iter and 0 have one end pointing to remote
			//			by induction:  any number can be expressed as sum of powers of 2.
			//			subproof:  number of nodes having one end pointing to remote doubles each iteration
			//		2. end node label is propagated to other end via 2^iter nodes in each iteration
			//		3. time: logarithmic in length of longest chain
			//		4. work: more complicated as chains drop out - upperbound is longest chain.

			// TODO: modify so that if left and right are not the same distance (shorter distance implies pointing to terminus
			//       we only send to both ends if left and right are the same distance.
			//		 we send to the longer end only if the distances are not the same.
			// 		 complexity?
			{
				if (comm.rank() == 0) printf("LIST RANKING\n");

				BL_BENCH_INIT(p_rank);
				// NOW: do the list ranking

				// search unfinished
				BL_BENCH_START(p_rank);
				::bliss::debruijn::filter::chain::PointsToInternalNode is_unfinished;
				auto unfinished = map.find(is_unfinished);
				BL_BENCH_COLLECTIVE_END(p_rank, "unfinished", unfinished.size(), comm);

				// check if all chains are compacted
				BL_BENCH_START(p_rank);
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				// set up chain node update metadata vector
				using update_md = bliss::debruijn::operation::chain::chain_update_md<KmerType>;
				std::vector<std::pair<KmerType, update_md > > updates;
				updates.reserve(unfinished.size());  // initially reserve same size.  will encounter a resize at least.

				uint dist = 0;
				KmerType ll, rr;
				bliss::debruijn::simple_biedge<KmerType> md;
				::bliss::debruijn::operation::chain::chain_update<KmerType> chain_updater;
				BL_BENCH_COLLECTIVE_END(p_rank, "init ranking", updates.capacity(), comm);

				// check if all chains are compacted
				BL_BENCH_START(p_rank);
				uint dist_in, dist_out;
				bool term_in, term_out;

				// loop until everything is compacted (except for cycles)
				while (!all_compacted) {

					updates.clear();

					// get left and right edges, generate updates
					for (auto t : unfinished) {
						md = t.second;

						// each is a pair with kmer, <in kmer, out kmer, in dist, out dist>
						// constructing 2 edges <in, out> and <out, in>  distance is sum of the 2.
						// indication of whether edge destination is a terminus depend only on the sign of distance to that node.
						dist_in = ::bliss::debruijn::get_chain_dist(std::get<2>(md));
						dist_out = ::bliss::debruijn::get_chain_dist(std::get<3>(md));
						term_in = ::bliss::debruijn::is_chain_terminal(std::get<2>(md));
						term_out = ::bliss::debruijn::is_chain_terminal(std::get<3>(md));
						dist = dist_in + dist_out;   // this double the distance...

						// and below pointer jumps.

						// construct forward edge, from ll to rr, only if current node is not a terminus for the "in" side
						if (!term_in)  {
							// send rr to ll.  also let ll know if rr is a terminus.  orientation is OUT
							updates.emplace_back(std::get<0>(md),
												update_md(term_out ? t.first : std::get<1>(md),
															::bliss::debruijn::points_to_chain_node(std::get<3>(md)) ? dist : ::bliss::debruijn::mark_as_point_to_terminal(dist),     // if md.3 <= 0, then finished, so use negative dist.
															bliss::debruijn::operation::OUT
														)
												);		// update target (md.0)'s out edge
							// if out dist is 0, then this node is a terminal node.  sent self as target.  else use right kmer.
							// if out dist is negative, then out kmer (rr) points to a terminus, including self (dist = 0), set update distance to negative to indicate so.

						}  // else case is same as below

						// construct backward edge, from out to in, only if current node is not a terminus for the "out" side
						if (!term_out) {
							// send ll to rr.  also let rr know if ll is a terminus.  orientation is IN
							updates.emplace_back(std::get<1>(md),
									update_md(term_in ? t.first : std::get<0>(md),
											::bliss::debruijn::points_to_chain_node(std::get<2>(md)) ? dist : ::bliss::debruijn::mark_as_point_to_terminal(dist),  // if md.3 <= 0, then finished, so use negative dist.
											bliss::debruijn::operation::IN));  // udpate target (md.1)'s in edge
							// if target is a terminus, then set self as target.  else use left kmer
							// if target points to a terminus, including self (dist = 0), set update distance to negative to indicate so.
						}  // else case is same as above


						// if ((std::get<2>(md) == 0) && (std::get<3>(md) == 0)) continue;  // singleton.   next.
					}

					comm.barrier();  // wait for all updates to be constructed.

					// now perform update
					size_t count = map.update( updates, false, chain_updater );

					// search unfinished.
					map.find(is_unfinished).swap(unfinished);

					// at this point, the new distances in lists are 2^(iterations + 1)
					++iterations;

					//std::cout << "rank " << comm.rank() << " iterations " << iterations << std::endl;
					// find cycles
					cycle_node_count = std::count_if(unfinished.begin(), unfinished.end(),
							::bliss::debruijn::filter::chain::IsCycleNode(iterations));

					// going over 30 makes the max_dist in IsCycleNode go to -1, then it is no longer valid as distances are int.  stop at 30
					// FORCE STOP AT iteration >= 30
					if (iterations >= 30) {

						// print the remaining non-cycle nodes locally.
						::bliss::debruijn::filter::chain::IsCycleNode check_cycle(iterations);
						for (auto t : unfinished) {
							if (check_cycle(t)) continue;

							auto md = t.second;

							std::cout << "rank " << comm.rank() << " max iter " << iterations <<
									"\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) <<
									"\tkmer: " << bliss::utils::KmerUtils::toASCIIString(t.first) << " rc: " << bliss::utils::KmerUtils::toASCIIString(t.first.reverse_complement()) <<
									"\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
						}

						//        	printf("rank %d max iter %lu updated %lu, unfinished %lu cycle nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_nodes);
						all_compacted = true;
					} else {

						all_compacted = (count == 0) || (cycle_node_count == unfinished.size());
						if (!all_compacted) printf("rank %d iter %lu updated %lu, unfinished %lu internal chain nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_node_count);
					}
					all_compacted = ::mxx::all_of(all_compacted, comm);

				}
				BL_BENCH_COLLECTIVE_END(p_rank, "compact", cycle_node_count, comm);

				BL_BENCH_REPORT_MPI_NAMED(p_rank, "list_rank", comm);
			}


//			// =========  this part should not be needed during list ranking, because all nodes should be covered due to proof 1.
//			// previous iteration, the unique source node is distance 2^(iter-1) toward the terminus
//			//
//			// ==== finally update all nodes that are not yet set to pointing to termini.
//			// should just be the ones that are pointing to termini but values are not yet marked as negative.
//			// previously has been sending out remote updates based on local information.
//			// due to pointer doubling, and local updating to the farthest remote jump,  some remote
//			// sources of updates may not be updated again, near the ends, even when they are already pointing to terminal
//			// these are resolved using a single query.  note that previous code terminates when nothing is updated or when only cycles remain.
//			{
//				if (comm.rank() == 0) printf("UPDATE ANY SKIPPED NODES (ones whose end points are handled by another kmer\n");
//
//				BL_BENCH_INIT(list_update);
//
//				// search unfinished
//				BL_BENCH_START(list_update);
//				auto unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
//				BL_BENCH_COLLECTIVE_END(list_update, "unfinished2", unfinished.size(), comm);
//
//				BL_BENCH_START(list_update);
//
//				// get global unfinished count
//				bool all_compacted = (unfinished.size() == 0);
//				all_compacted = ::mxx::all_of(all_compacted, comm);
//
//				// while have unfinished,  run.  qq contains kmers not necessarily canonical
//
//				std::vector<KmerType> qq;
//				qq.reserve(unfinished.size() * 2);
//
//				if (!all_compacted) {
//
//					qq.clear();
//					// get the end points.
//					bliss::debruijn::simple_biedge<KmerType> md;
//
//					for (auto t : unfinished) {
//						md = t.second;
//
//						// add the in edge target if it needs to be updated.
//						if (std::get<2>(md) > 0) {
//							qq.emplace_back(std::get<0>(md));
//						}
//						if (std::get<3>(md) > 0) {
//							qq.emplace_back(std::get<1>(md));
//						}
//					}
//
//					comm.barrier();
//
//					// now query.
//					auto results = chainmap.find(qq);
//
//					// put query results in a map.  key is CANONICAL
//					typename ChainMapType::local_container_type res_map(results.begin(), results.end());
//
//					// now update.  for each unfinished
//					KmerType kk;
//					::bliss::kmer::transform::lex_less<KmerType> lexless;
//					for (auto t : unfinished) {
//						// get left (in) kmer.
//						kk = std::get<0>(t.second);
//
//						// lookup left by canonical
//						auto it = res_map.find(lexless(kk));
//
//						// if left is at terminus, and we haven't updated it,
//						if (((std::get<2>(it->second) == 0) || (std::get<3>(it->second) == 0)) &&
//								(std::get<2>(t.second) > 0)) {
//							// update the current left
//							std::get<2>((*(chainmap.get_local_container().find(t.first))).second) = -(std::get<2>(t.second));
//						}
//
//						// get right (out) kmer
//						kk = std::get<1>(t.second);
//
//						// lookup left by canonical
//						it = res_map.find(lexless(kk));
//
//						// if right is at terminus, and we haven't updated it,
//						if (((std::get<2>(it->second) == 0) || (std::get<3>(it->second) == 0)) &&
//								(std::get<3>(t.second) > 0)) {
//							// update the current left
//							std::get<3>((*(chainmap.get_local_container().find(t.first))).second) = -(std::get<3>(t.second));
//						}
//
//					}
//
//				}
//				BL_BENCH_COLLECTIVE_END(list_update, "cleanup", unfinished.size(), comm);
//
//				size_t unfin = mxx::allreduce(unfinished.size(), comm);
//				if (comm.rank() == 0) printf("FINAL UPDATE for %ld nodes\n", unfin);
//
//				BL_BENCH_REPORT_MPI_NAMED(list_update, "finalize", comm);
//			}

			{

				auto unfinished = map.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				if (!all_compacted)
					separate_cycles(iterations);

			}
			{
				// VERIFY ALL DONE.
				auto unfinished = map.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				assert(all_compacted);
			}

			// separate out the isolated.
			separate_isolated();

			listranked = true;
			modified = false;

			return iterations;
		}



	   /**
	    * @brief 	processes the chains so that each points to the 2 ends of the chain, along with distances.
	    * @details	uses parallel list ranking algorithm internally.
	    * 			any cycles found (by checking that all remaining active nodes have left and right distance
	    * 			that are 2^iter.
	    *
	    * @note		this is the version used for initial list_ranking.
	    * 			this function would be inefficient for recompacting after some change.
	    *
	    * 			this is the version that implements optimization 3 below.
	    * 				for this there are 5 types of edges:
	    * 					1. neither edges are finished. so update both sides.  this has work 2 * (N - 2 * 2^iter), excluding half finished nodes and 1 for each edge.
	    * 						both > 0
	    * 					2. sending local or non-local terminal to non-terminal.  should always do this.  - this has work 2 * 2^iter and double the number of nodes with 1 side terminal.
	    * 						1 side <= 0, other > 0.  send to >0 side.  terminal node goes to 2^iter
	    * 					3. sending non-terminal to terminal node (to get terminal to point to other terminal).  in list_rank above,
	    * 						we always send, so require 2 * 2^iter work.
	    * 						but since we know we have to do a reduction for largest dist, which can be at most 2^iter,
	    * 							we can choose dist of 2*2^iter node to send.
	    *
	    * 						1 side < 0, other > 0.  if dist == 2 * 2^iter, send to < 0 side.  note that dist = 2^(iter+1) implies (l/r)dist != 0, so we can check for <= like in 2.
	    * 					4. at the last iteration, where 2*2^iter > length, criteria 3 may not be satisfied, but we will have nodes with both edges pointing to termini.
	    * 						there should be 1 where the left and right distances are equal, or 2 with left and right distances differ by 1.
	    * 						so in this iteration we will have at most 2 messages.
	    * 						both sides < 0.  send to both sides (left dist == right dist, or left dist == right dist + 1)  (middle of chain, so not going to be == 0)
	    *
	    * 						THIS IS NEVER SATISFIED BECAUSE WE AE WORKING ON UNFINISHED ONLY.
	    * 					5. both sides 0.  isolated node. do nothing since already done.
	    *
	    * 				first 2 types are normal.
	    * 			comm savings should be i = 0..log(N), sum(2 * 2^i) = 2N over the entire run, and compute savings is up to N during last stages

		*           TODO: finding unfinished can be merged with the actual loop, and cycle count can be accumulated without using a separate loop.
	    */
	   size_t list_rank_min_update() {

			size_t iterations = 0;
			size_t cycle_node_count = 0;

			// optimizations:
			//		1. only send to non-end edges.  reduces update comm volume for termini.  this is probably normal
			//		2. only send if a node is active (i.e. chain ranked).   reduces number of active chains.
			//		3. only send to the unfinished side(s), or to the finished side if the total distance is expected to be max (2^(iter + 1), or entire length of chain.)
			// to prove:
			//		0. distance between left and right doubles each iteration. (middle nodes)
			//		1. all nodes between 2^iter and 0 have one end pointing to remote
			//			by induction:  any number can be expressed as sum of powers of 2.
			//			subproof:  number of nodes having one end pointing to remote doubles each iteration
			//		2. end node label is propagated to other end via 2^iter nodes in each iteration
			//		3. time: logarithmic in length of longest chain
			//		4. work: more complicated as chains drop out - upperbound is longest chain.
			// NOTE:    it takes an extra iteration for an edge to get negative distance after pointer to terminal
			//      target of terminal edge is + 2^iter away.
			//      in next iteration, terminal is responsible to negate that edge.
			//      between [0, 2^iter), distances are 0, updated by nodes [0, 2^(iter-1) ).
			//      so checking for one side being negative and the total being 2^(iter+1) is not correct.
			//        the node at 2^(iter) will update terminal, either to 2^(iter+1) if it does not point to other terminal, else a half completed node would send out only, and that'd be okay too.
			// TODO: modify so that if left and right are not the same distance (shorter distance implies pointing to terminus
			//       we only send to both ends if left and right are the same distance.
			//		 we send to the longer end only if the distances are not the same.
			// 		 complexity?
			{
				if (comm.rank() == 0) printf("LIST RANKING\n");

				BL_BENCH_INIT(p_rank);
				// NOW: do the list ranking

				// search unfinished
				BL_BENCH_START(p_rank);
				::bliss::debruijn::filter::chain::PointsToInternalNode is_unfinished;
				auto unfinished = map.find_iterators(is_unfinished);
				BL_BENCH_COLLECTIVE_END(p_rank, "unfinished", unfinished.size(), comm);

				// check if all chains are compacted
				BL_BENCH_START(p_rank);
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				// set up chain node update metadata vector
				using update_md = bliss::debruijn::operation::chain::chain_update_md<KmerType>;
				std::vector<std::pair<KmerType, update_md > > updates;
				updates.reserve(unfinished.size());  // initially reserve a large block - at least same number as unfinisehd
						// will grow

				if (comm.rank() == 0)
				std::cout << "SIZES chain_update_md size is " << sizeof(update_md) << std::endl;

				uint dist = 0, ldist = 0, rdist = 0;
				KmerType ll, rr;
				bliss::debruijn::simple_biedge<KmerType> md;
				KmerType km;
				::bliss::debruijn::operation::chain::chain_update<KmerType> chain_updater;
				BL_BENCH_COLLECTIVE_END(p_rank, "init ranking", updates.capacity(), comm);

				// loop until everything is compacted (except for cycles)
				while (!all_compacted) {

					updates.clear();

					// get left and right edges, generate updates
					for (auto t : unfinished) {
						md = (*t).second;
						km = (*t).first;

						// each is a pair with kmer, <in kmer, out kmer, in dist, out dist>
						// constructing 2 edges <in, out> and <out, in>  distance is sum of the 2.
						// indication of whether edge destination is a terminus depend only on the sign of distance to that node.

						ldist = std::get<2>(md);
						rdist = std::get<3>(md);
						dist = ::bliss::debruijn::get_chain_dist(ldist) + ::bliss::debruijn::get_chain_dist(rdist);   // this double the distance...

						// and below pointer jumps.

						// switch based on node type
						if (::bliss::debruijn::points_to_chain_node(ldist) && ::bliss::debruijn::points_to_chain_node(rdist)) {
							// type 1.
							// update left
							// send rr to ll.  also let ll know if rr is a terminus.  orientation is OUT
							updates.emplace_back(std::get<0>(md),
									update_md(std::get<1>(md),
												dist,     // if md.3 <= 0, then finished, so use negative dist.
											bliss::debruijn::operation::OUT));		// update target (md.0)'s out edge
							// if out dist is 0, then this node is a terminal node.  sent self as target.  else use right kmer.
							// if out dist is negative, then out kmer (rr) points to a terminus, including self (dist = 0), set update distance to negative to indicate so.

							// udpate right.
							updates.emplace_back(std::get<1>(md),
									update_md(std::get<0>(md),
											 	 dist,  // if md.3 <= 0, then finished, so use negative dist.
											bliss::debruijn::operation::IN));  // udpate target (md.1)'s in edge

						} else if (::bliss::debruijn::points_to_terminal(ldist) && ::bliss::debruijn::points_to_terminal(rdist)) {
//							// type 4.  because of filter, this would never be called.
//							if ((ldist == rdist) || (ldist == (rdist - 1))) {
//								// update left
//								updates.emplace_back(std::get<0>(md),
//										update_md(std::get<1>(md),
//													-dist,     // if md.3 <= 0, then finished, so use negative dist.
//												bliss::debruijn::operation::OUT));		// update target (md.0)'s out edge
//								// update right.
//								updates.emplace_back(std::get<1>(md),PointsToInternalNode
//										update_md(std::get<0>(md),
//													-dist,  // if md.3 <= 0, then finished, so use negative dist.
//												bliss::debruijn::operation::IN));  // udpate target (md.1)'s in edge
//
//							}
			              throw std::logic_error("both distances are negative. should be filtered already and not get here");

						} else if (dist == 0){  // types 5. do nothing. both distances are 0.  never be called either.
						} else {  // types 2 and 3.  one of 2 edges points to terminus but not the other.
							// ((ldist <= 0) && (rdist > 0)) || ((ldist < 0) && (rdist >= 0)) ||
							// ((ldist > 0) && (rdist <= 0)) || ((ldist >= 0) && (rdist < 0))
							// handle type 2 here
						  // should only happen for nodes that are less than 2^iter away from terminal.
							if (::bliss::debruijn::points_to_chain_node(rdist)) {  // rdist > 0 and ldist <= 0, in side is done.
								// left is terminus.  send to right then
								// update right.
								updates.emplace_back(std::get<1>(md),
										update_md((::bliss::debruijn::get_chain_dist(ldist) == 0) ? km : std::get<0>(md),
												::bliss::debruijn::mark_as_point_to_terminal(dist),  // if md.3 <= 0, then finished, so use negative dist.
												bliss::debruijn::operation::IN));  // udpate target (md.1)'s in edge

								// handle type 3 here - these are never satistfied.
								if (dist == (static_cast<uint>(2) << iterations)) { // iter i updates a node at max of 2^(i+1) distance away.
//									// update left
//									updates.emplace_back(std::get<0>(md),
//											update_md(std::get<1>(md),  // this cannot be a terminus
//													dist,     // if md.3 <= 0, then finished, so use negative dist.
//													bliss::debruijn::operation::OUT));		// update target (md.0)'s out edge
								  throw std::logic_error("total distance at 2^t and trying to update terminal. should not get here");
								}

							} else {  // ldist > 0 and rdist <= 0
								// right is terminus.  send to left then
								// update left
								updates.emplace_back(std::get<0>(md),
										update_md((::bliss::debruijn::get_chain_dist(rdist) == 0) ? km : std::get<1>(md),
												::bliss::debruijn::mark_as_point_to_terminal(dist),     // if md.3 <= 0, then finished, so use negative dist.
												bliss::debruijn::operation::OUT));		// update target (md.0)'s out edge
								// handle type 3 here
								if (dist == (static_cast<uint>(2) << iterations)) { // iter i updates a node at max of 2^(i+1) distance away.
//									// update right.
//									updates.emplace_back(std::get<1>(md),
//											update_md(std::get<0>(md),
//													dist,  // if md.3 <= 0, then finished, so use negative dist.
//													bliss::debruijn::operation::IN));  // udpate target (md.1)'s in edge
				                  throw std::logic_error("total distance at 2^t and trying to update terminal. should not get here");
								}
							}
						}


						// if ((std::get<2>(md) == 0) && (std::get<3>(md) == 0)) continue;  // singleton.   next.
					}

//					comm.barrier();  // wait for all updates to be constructed.

					// now perform update
					size_t count = map.update( updates, false, chain_updater );

					// search unfinished.
					//map.find(is_unfinished).swap(unfinished);  // this is scanning the whole list every time...
					auto new_end = ::std::partition(unfinished.begin(), unfinished.end(), is_unfinished);
//							[&is_unfinished](typename std::iterator_traits<decltype(unfinished)>::value_type t){
//						return is_unfinished(*t);
//					});
					unfinished.erase(new_end, unfinished.end());   // remove d the finished part.


					// at this point, the new distances in lists are 2^(iterations + 1)
					++iterations;

					//std::cout << "rank " << comm.rank() << " iterations " << iterations << std::endl;
					// find cycles
					cycle_node_count = std::count_if(unfinished.begin(), unfinished.end(),
							::bliss::debruijn::filter::chain::IsCycleNode(iterations));

					// going over 30 makes the max_dist in IsCycleNode go to -1, then it is no longer valid as distances are int.  stop at 30
					// FORCE STOP AT iteration >= 30
					if (iterations >= 30) {

						// print the remaining non-cycle nodes locally.
						::bliss::debruijn::filter::chain::IsCycleNode check_cycle(iterations);
						for (auto t : unfinished) {
							if (check_cycle(*t)) continue;

							auto md = (*t).second;

							std::cout << "rank " << comm.rank() << " max iter " << iterations <<
									"\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) <<
									"\tkmer: " << bliss::utils::KmerUtils::toASCIIString((*t).first) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString((*t).first.reverse_complement()) <<
									"\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
						}

						//        	printf("rank %d max iter %lu updated %lu, unfinished %lu cycle nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_nodes);
						all_compacted = true;
					} else {

						all_compacted = (count == 0) || (cycle_node_count == unfinished.size());
						if (!all_compacted) printf("rank %d iter %lu updated %lu, unfinished %lu chain internal nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_node_count);
					}
					all_compacted = ::mxx::all_of(all_compacted, comm);

				} // while
				BL_BENCH_COLLECTIVE_END(p_rank, "compact", cycle_node_count, comm);

				BL_BENCH_REPORT_MPI_NAMED(p_rank, "list_rank", comm);
			}

			{

				auto unfinished = map.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				if (!all_compacted)
					separate_cycles(iterations);

			}
			{
				// VERIFY ALL DONE.
				auto unfinished = map.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				assert(all_compacted);
			}

			// separate out the isolated.
			separate_isolated();


			listranked = true;
			modified = false;

			return iterations;
		}


		// alternative list ranking algorithm that defines cycle nodes as ones with left and right edges NOT pointing to termini.
		// when all unfinished are cycle nodes, then we terminate.
	   size_t list_rank2() {

			size_t iterations = 0;
			size_t terminal_count = 0;

			// optimizations:
			//		1. only send to non-end edges.  reduces update comm volume for termini.  this is probably normal
			//		2. only send if we are active (i.e. chain ranked).   reduces number of active chains.
			//		3. only send to the unfinished side(s), or to the finished side if the total distance is expected to be max (2^(iter + 1), or entire length of chain.)
			//		4. use vector of iterators to local chain map.  we are not changing the local chain map, so nothing is invalidated.
			//			then we can operate on the unfinished with partition, etc, also should save memory.
			// to prove:
			//		0. distance between left and right doubles each iteration. (middle nodes)
			//		1. all nodes between 2^iter and 0 have one end pointing to remote
			//			by induction:  any number can be expressed as sum of powers of 2.
			//			subproof:  number of nodes having one end pointing to remote doubles each iteration
			//		2. end node label is propagated to other end via 2^iter nodes in each iteration
			//		3. time: logarithmic in length of longest chain
			//		4. work: more complicated as chains drop out - upperbound is longest chain.

			// TODO: modify so that if left and right are not the same distance (shorter distance implies pointing to terminus
			//       we only send to both ends if left and right are the same distance.
			//		 we send to the longer end only if the distances are not the same.
			// 		 complexity?
			{
				if (comm.rank() == 0) printf("LIST RANKING\n");

				BL_BENCH_INIT(p_rank);
				// NOW: do the list ranking

				// search unfinished
				BL_BENCH_START(p_rank);
				::bliss::debruijn::filter::chain::PointsToInternalNode is_unfinished;  // left or right not pointing to termus.
				auto unfinished = map.find_iterators(is_unfinished);
				BL_BENCH_COLLECTIVE_END(p_rank, "unfinished", unfinished.size(), comm);

				// check if all chains are compacted
				BL_BENCH_START(p_rank);
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				// set up chain node update metadata vector
				using update_md = bliss::debruijn::operation::chain::chain_update_md<KmerType>;
				std::vector<std::pair<KmerType, update_md > > updates;
				updates.reserve(unfinished.size());  // initially reserve same size.  will encounter a resize at least.

				uint dist = 0;
				KmerType ll, rr;
				bliss::debruijn::simple_biedge<KmerType> md;
				KmerType km;
				::bliss::debruijn::operation::chain::chain_update<KmerType> chain_updater;
				BL_BENCH_COLLECTIVE_END(p_rank, "init ranking", updates.capacity(), comm);

				// check if all chains are compacted
				BL_BENCH_START(p_rank);
				uint dist_in, dist_out;
				bool term_in, term_out;

				
				::bliss::debruijn::filter::chain::IsTerminusOrIsolated is_semifinished;
				// loop until everything is compacted (except for cycles)
				while (!all_compacted) {

					updates.clear();

					// get left and right edges, generate updates
					for (auto t : unfinished) {
						md = (*t).second;
						km = (*t).first;

						// each is a pair with kmer, <in kmer, out kmer, in dist, out dist>
						// constructing 2 edges <in, out> and <out, in>  distance is sum of the 2.
						// indication of whether edge destination is a terminus depend only on the sign of distance to that node.
						dist_in = ::bliss::debruijn::get_chain_dist(std::get<2>(md));
						dist_out = ::bliss::debruijn::get_chain_dist(std::get<3>(md));
						term_in = ::bliss::debruijn::is_chain_terminal(std::get<2>(md));
						term_out = ::bliss::debruijn::is_chain_terminal(std::get<3>(md));
						dist = dist_in + dist_out;   // this double the distance...

						// and below pointer jumps.

						// construct forward edge, from ll to rr, only if current node is not a terminus for the "in" side
						if (!term_in)  {
							// send rr to ll.  also let ll know if rr is a terminus.  orientation is OUT
							updates.emplace_back(std::get<0>(md),
												update_md(term_out ? km : std::get<1>(md),
															::bliss::debruijn::points_to_chain_node(std::get<3>(md)) ? dist : ::bliss::debruijn::mark_as_point_to_terminal(dist),     // if md.3 <= 0, then finished, so use negative dist.
															bliss::debruijn::operation::OUT
														)
												);		// update target (md.0)'s out edge
							// if out dist is 0, then this node is a terminal node.  sent self as target.  else use right kmer.
							// if out dist is negative, then out kmer (rr) points to a terminus, including self (dist = 0), set update distance to negative to indicate so.

						}  // else case is same as below

						// construct backward edge, from out to in, only if current node is not a terminus for the "out" side
						if (!term_out) {
							// send ll to rr.  also let rr know if ll is a terminus.  orientation is IN
							updates.emplace_back(std::get<1>(md),
									update_md(term_in ? km : std::get<0>(md),
											::bliss::debruijn::points_to_chain_node(std::get<2>(md)) ? dist : ::bliss::debruijn::mark_as_point_to_terminal(dist),  // if md.3 <= 0, then finished, so use negative dist.
											bliss::debruijn::operation::IN));  // udpate target (md.1)'s in edge
							// if target is a terminus, then set self as target.  else use left kmer
							// if target points to a terminus, including self (dist = 0), set update distance to negative to indicate so.
						}  // else case is same as above


						// if ((std::get<2>(md) == 0) && (std::get<3>(md) == 0)) continue;  // singleton.   next.
					}

					//comm.barrier();  // wait for all updates to be constructed.

					// now perform update
					size_t count = map.update( updates, false, chain_updater );

					// search unfinished.
					//map.find(is_unfinished).swap(unfinished);  // this is scanning the whole list every time...
					auto new_end = ::std::partition(unfinished.begin(), unfinished.end(), is_unfinished);
//							[&is_unfinished](typename std::iterator_traits<decltype(unfinished)>::value_type t){
//						return is_unfinished(*t);
//					});
					unfinished.erase(new_end, unfinished.end());   // remove d the finished part.

					// at this point, the new distances in lists are 2^(iterations + 1)
					++iterations;

					//std::cout << "rank " << comm.rank() << " iterations " << iterations << std::endl;
					// find cycles
					terminal_count = std::count_if(unfinished.begin(), unfinished.end(),
							is_semifinished);

					// going over 30 makes the max_dist in IsCycleNode go to -1, then it is no longer valid as distances are int.  stop at 30
					// FORCE STOP AT iteration >= 30
					if (iterations >= 30) {
						::bliss::debruijn::filter::chain::IsCycleNode check_cycle(iterations);
						// print the remaining non-cycle nodes locally.
						for (auto t : unfinished) {
							if (check_cycle(*t)) continue;

							auto md = (*t).second;

							std::cout << "rank " << comm.rank() << " max iter " << iterations <<
									"\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) <<
									"\tkmer: " << bliss::utils::KmerUtils::toASCIIString((*t).first) << 
									" rc: " << bliss::utils::KmerUtils::toASCIIString((*t).first.reverse_complement()) <<
									"\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
						}

						//        	printf("rank %d max iter %lu updated %lu, unfinished %lu cycle nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_nodes);
						all_compacted = true;
					} else {

						all_compacted = (count == 0) || (terminal_count == 0);
						if (!all_compacted) printf("rank %d iter %lu updated %lu, unfinished %lu terminals of unfinished chains %lu\n", comm.rank(), iterations, count, unfinished.size(), terminal_count);
					}
					all_compacted = ::mxx::all_of(all_compacted, comm);

				}
				BL_BENCH_COLLECTIVE_END(p_rank, "compact", unfinished.size(), comm);

				BL_BENCH_REPORT_MPI_NAMED(p_rank, "list_rank", comm);
			}


			listranked = true;
			modified = false;

			return iterations;
		}

		void separate_isolated_and_cycles() {

			// separate out the isolated.
			separate_isolated();

			{
				auto unfinished = map.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				if (!all_compacted)
					separate_cycles();

			}
			{
				// VERIFY ALL DONE.
				auto unfinished = map.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
				bool all_compacted = (unfinished.size() == 0);
				all_compacted = ::mxx::all_of(all_compacted, comm);

				assert(all_compacted);
			}

		}


		/**
		 *  merge one chain into another.
		 */
		void merge(debruijn_chain_graph const & other) {  // only local operations since dist hash is identical.
			auto it = other.get_map().get_local_container().cbegin();
			auto end = other.get_map().get_local_container().cend();
			auto lend = map.get_local_container().end();
			for (; it != end; ++it) {
				auto found = map.get_local_container().find((*it).first);
				if (found != lend) {  // if existing, replace value
					(*found).second = (*it).second;
				} else {  // insert if previous not there.
					map.get_local_container().insert((*it));
				}
			}

			// merge the cycle nodes
			it = other.cycle_nodes.get_local_container().cbegin();
			end = other.cycle_nodes.get_local_container().cend();
			for (; it != end; ++it) {
				cycle_nodes.get_local_container().insert(*it);
			}

			// merge the isolated nodes
			it = other.isolated.get_local_container().cbegin();
			end = other.isolated.get_local_container().cend();
			for (; it != end; ++it) {
				isolated.get_local_container().insert(*it);
			}
			
		}

	   //====== OUTPUT....


/**
 * 	convert t oranked chain
 */
	   std::vector<::bliss::debruijn::chain::listranked_chain_node<KmerType> >
	   to_ranked_chain_nodes() {

		   using chain_vec = std::vector<::bliss::debruijn::chain::listranked_chain_node<KmerType> >;
		   chain_vec compacted_chain;
		   compacted_chain.reserve(map.local_size());
		   ::fsc::back_emplace_iterator<chain_vec > back_emplacer(compacted_chain);

			//== first transform nodes so that we are pointing to canonical terminus k-mers.
			std::transform(map.get_local_container().cbegin(), map.get_local_container().cend(), back_emplacer,
					::bliss::debruijn::operation::chain::to_listranked_chain_node<KmerType>());

			return compacted_chain;
	   }

/**
 * 	compress chain
 */
	   std::vector<::std::string> to_compressed_chains() {
	     // ========== construct new graph with compacted chains and junction nodes.
	     BL_BENCH_INIT(compress_chain);

	     if (comm.rank() == 0) printf("Compress Chain\n");

	     auto compacted_chains = to_ranked_chain_nodes();

	     int has_data = (compacted_chains.size() == 0) ? 0 : 1;
	     int all_has_data = mxx::allreduce(has_data, comm);
	     if (all_has_data == 0) {
		     BL_BENCH_REPORT_MPI_NAMED(compress_chain, "compress_chain", comm);

	    	 return std::vector<::std::string>();
	     }


	     if (comm.size() > 1) {
			// reshuffle data
			BL_BENCH_START(compress_chain);
			::bliss::debruijn::operation::chain::chain_node_to_proc<map_params_template, KmerType> mapper(comm.size());
			// distribute the data

			std::vector<size_t> recv_counts;
     		 std::vector<::bliss::debruijn::chain::listranked_chain_node<KmerType> > distributed_chains;
//					::dsc::assign_and_bucket(compacted_chains, mapper, this->comm.size());
//			BL_BENCH_END(compress_chain, "assign and bucket", compacted_chains.size());
//
//			// and distribute the data.
//			BL_BENCH_START(compress_chain);
//			::dsc::distribute_bucketed(compacted_chains, recv_counts, this->comm).swap(compacted_chains);
			::imxx::distribute(compacted_chains, mapper, recv_counts, distributed_chains, comm);
      		distributed_chains.swap(compacted_chains);
			BL_BENCH_COLLECTIVE_END(compress_chain, "distribute nodes", compacted_chains.size(), comm);  // this is for output ordering.
	     }

	     // next local sort the data by terminus kmer and position
	     BL_BENCH_START(compress_chain);
	     ::std::sort(compacted_chains.begin(), compacted_chains.end(), ::bliss::debruijn::operation::chain::chain_rep_less<KmerType>());
	     BL_BENCH_COLLECTIVE_END(compress_chain, "sort lmer nodes", compacted_chains.size(), comm);   // this is for constructing the chains


	     BL_BENCH_START(compress_chain);
	     std::vector<::std::string> compressed_chains;

	     // compressing.
	     std::stringstream ss;
	     //=== approach is to scan for the end of a chain
	     auto curr = compacted_chains.begin();
	     auto next = curr;
	     compressed_chains.clear();
	     while (curr != compacted_chains.end()) {
	       // scan forward for start of next chain, using adjacent find.
	       next = ::std::adjacent_find(curr, compacted_chains.end(), [](
	         ::bliss::debruijn::chain::listranked_chain_node<KmerType> const & x,
	          ::bliss::debruijn::chain::listranked_chain_node<KmerType> const & y
	       ){
	         // adjacent_find returns first occurrence of 2 consecutive elements that the predicate evaluates to true.
	         return std::get<1>(x) != std::get<1>(y);
	       });

	       if (next != compacted_chains.end()) ++next;  // get the second of the pair, == start of next.

	       // now do something with the <curr, next> range.
	       ss.clear();
	       ss.str(std::string());
	       std::for_each(curr, next, ::bliss::debruijn::operation::chain::print_chain_as_fasta<KmerType>(ss, true));  // chain_only
	       compressed_chains.push_back(ss.str());

	       // get ready for next segment.
	       curr = next;

	     }
	     BL_BENCH_COLLECTIVE_END(compress_chain, "compress", compressed_chains.size(), comm);

	     BL_BENCH_REPORT_MPI_NAMED(compress_chain, "compress_chain", comm);

	     return compressed_chains;
	   }



		// make a separate chain_graph that only has terminal chain nodes.  exclude completely isolated nodes, but leave in unit chain nodes.(chain length of 1.)
	   void make_terminal_chain_graph(debruijn_chain_graph<KmerType, DistHash> & termini) {
		   // assume cycles and isolated are excluded.
		   auto temp = this->get_terminal_nodes();
			termini.get_map().get_local_container().insert(temp);  // local insert only.
	   }
		void make_representative_chain_graph(debruijn_chain_graph<KmerType, DistHash> & reps) {
		   // assume cycles and isolated are excluded.
		   auto temp = this->get_chain_representatives();
			reps.get_map().get_local_container().insert(temp); //, ::bliss::transform::identity<mutable_value_type>(), true);  // local insert only.
	   }

	/**
	 * @brief return a vector of chain terminals from deadends.  isolated chain nodes are ignored as they are NOT deadends.
	 * 		
	 */
	std::vector<::bliss::debruijn::chain::summarized_chain<KmerType> > to_summarized_chains() const {
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

		std::vector<mutable_value_type> termini = this->get_terminal_nodes();
		std::vector<::bliss::debruijn::chain::summarized_chain<KmerType> > results;

		// evenly redistribute.
		mxx::distribute_inplace(termini, comm);
		
		// participate only if there is at least one entry.
		bool has_data = (termini.size() > 0);
		bool all_has_data = ::mxx::all_of(has_data, comm);
		mxx::comm subcomm = all_has_data ? comm.copy() : comm.split(has_data);

		if (has_data) {
			// now sort by chain representatives
			mxx::sort(termini.begin(), termini.end(),
				[](mutable_value_type const & lhs, mutable_value_type const & rhs){
					// get chain representatives, then compare.
					return ::bliss::debruijn::chain::get_chain_rep(lhs) < ::bliss::debruijn::chain::get_chain_rep(rhs);
			}, subcomm);
			
			// termini comes in pairs, except for isolated kmers.  ship first element to prev
			mutable_value_type first = termini.front();
			auto first_chain_rep = ::bliss::debruijn::chain::get_chain_rep(first);

			mutable_value_type second = mxx::left_shift(first, subcomm);
			auto curr_chain_rep = ::bliss::debruijn::chain::get_chain_rep(second);
			// if same as last, then insert to the back.  would not introduce a new chain, or duplicate an isolated.
			if (::bliss::debruijn::chain::get_chain_rep(termini.back()) == curr_chain_rep) {
				termini.emplace_back(second);
			}  // this means that the last entries will be isolated, or paired.

			// insert the first one in the chain if it is length 1 (guarantee next is going to have different chain rep)
			if (::bliss::debruijn::is_chain_terminal(std::get<2>(first.second)) &&  // if both 2 and 3 are termini, then this is an isolated one.
				::bliss::debruijn::is_chain_terminal(std::get<3>(first.second)) ) { // Forward, Reverse.
					results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
						first.first, first.second,
						first.first, first.second
						)); 
				}  // else the first one is actually the second terminus, or the first terminus.

			// iterate over all termini to summarize deadend chains
			size_t i = 1;
			// then scan to locate deadends.  we must have a start and an end.  compact in place.  isolated are ignored.  process current 
			for (; i < termini.size(); ++i) {

				curr_chain_rep = ::bliss::debruijn::chain::get_chain_rep(termini[i]);
				// decide if we are at first, or second terminus
				if (curr_chain_rep != first_chain_rep ) {  // 
					first = termini[i];   // new chain, reset first.
					first_chain_rep = curr_chain_rep;

					// insert if length 1.  next one is guaranteed to have different chain rep.
					if (::bliss::debruijn::is_chain_terminal(std::get<2>(first.second)) &&  // if both 2 and 3 are termini, then this is an isolated one.
						::bliss::debruijn::is_chain_terminal(std::get<3>(first.second)) ) { // Forward, Reverse.
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								first.first, first.second,
								first.first, first.second
								)); 
						}  // mutually exclusive with successive "chain rep" being equal...
				}  // end case with switching to new chain terminus (potentially pairs). 
				else {
					second = termini[i];   // found the termini of 1 chain.  

					// begin conversion.  check ordering.

					// <kmer1, <L1, R1, uint, uint> >   and <kmer2, <L2, R2, uint, uint> >.  the chain reps are same, so termini k-mers are same.
					// just need to determine the order. within terminus we are on the same strand.  between termini, terminal dist plus chain rep can be used to infer orientation
					// one of these has the chain representative as the kmer.
					if (::bliss::debruijn::is_chain_terminal(std::get<2>(first.second)) &&  // 2 and 3 are mutually exclusive within a terminus
						::bliss::debruijn::is_chain_terminal(std::get<3>(second.second)) ) { // Forward, Reverse.
						if (first.first == curr_chain_rep) {  // and first is rep.  order is 1, 2
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								first.first, first.second,
								second.first, second.second
								)); 
						} else if (second.first.reverse_complement() == curr_chain_rep) { // and second is rep.  order is rc(2), rc(1)
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								second.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(second.second),
								first.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(first.second)
								)); 
						} else {
							throw std::logic_error("ERROR! chain rep does not match the kmers on the expected strand. case 1");
						}
					} else if (::bliss::debruijn::is_chain_terminal(std::get<2>(first.second)) &&  // 2 and 3 are mutually exclusive within a terminus
						::bliss::debruijn::is_chain_terminal(std::get<2>(second.second)) ) { // Forward, Forward.
						if (first.first == curr_chain_rep) {  // and first is rep order is 1, rc(2)
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								first.first, first.second,
								second.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(second.second)
								));
						} else if (second.first == curr_chain_rep) { // and second is rep, then order is 2, rc(1)
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								second.first, second.second,
								first.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(first.second)
								));
						} else {
							throw std::logic_error("ERROR! chain rep does not match the kmers on the expected strand. case 2");
						}
					} else if (::bliss::debruijn::is_chain_terminal(std::get<3>(first.second)) &&  // 2 and 3 are mutually exclusive within a terminus
						::bliss::debruijn::is_chain_terminal(std::get<3>(second.second)) ) { // Reverse, Reverse.
						if (first.first.reverse_complement() == curr_chain_rep) {  // and first is rep.  order is rc(1), 2
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								first.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(first.second),
								second.first, second.second
								)); 
						
						} else if (second.first.reverse_complement() == curr_chain_rep) { // and second is rep.  order is rc(2), 1
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								second.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(second.second),
								first.first, first.second
								));
						} else {
							throw std::logic_error("ERROR! chain rep does not match the kmers on the expected strand. case 3");

						}
					} else if (::bliss::debruijn::is_chain_terminal(std::get<3>(first.second)) &&  // 2 and 3 are mutually exclusive within a terminus
						::bliss::debruijn::is_chain_terminal(std::get<2>(second.second)) ) { // reverse, Forward.
						if (first.first.reverse_complement() == curr_chain_rep) {  // and first is rep order is rc(1), rc(2)
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								first.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(first.second),
								second.first.reverse_complement(), ::bliss::debruijn::transform::reverse_complement(second.second)
								)); 

						} else if (second.first == curr_chain_rep) { // and second is rep, then order is 2, 1
							results.emplace_back(::bliss::debruijn::chain::make_summarized_chain(
								second.first, second.second,
								first.first, first.second
								)); 

						} else {
							throw std::logic_error("ERROR! chain rep does not match the kmers on the expected strand. case 4");
						}
					} else {
						throw std::logic_error("ERROR! orientation of the two termini is confused.");
					}

					
				} // end case for chain with 2 termini

			}

		}
		return results;
	}


	/**
	 * @brief return a vector of chain terminals from deadends.  isolated chain nodes are ignored as they are NOT deadends.
	 * 		
	 */
	template <typename DBG>
	std::vector<::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> > to_summarized_chains(DBG const & dbg) const {
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

		// need edge frequency.  Instead of sorting termini then do pair conversion,
		//   convert, sort, then merge.  both termini should share the same chain rep, so sort by chain rep.  mergeing should be a max operation on edge freq.



		std::vector<::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> > results;

		{
			std::vector<mutable_value_type> termini = this->get_terminal_nodes();
			// convert termini to summarized chain
			::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> summary;
			
			for (mutable_value_type terminus : termini) {
				// get the node from dbg.  note that the processing is simplified because 
					// 1. we are targetting chains, so frequency extraction is easy.
					// 2. the chain map key is distributed the same way as dbg's keys, so these operations are local
					// 3. canonicalization are identical, so we can directly query.
					// 4. relies on freq of edge being recorded the same way on both kmers.
				auto it = dbg.get_map().get_local_container().find(terminus.first);
				assert(((*it).second.get_in_edge_count() < 2) && ((*it).second.get_out_edge_count() < 2));  // should be an edge

				// if this is an isolated or unit chain,
				if (::bliss::debruijn::is_chain_terminal(std::get<2>(terminus.second)) && 
					::bliss::debruijn::is_chain_terminal(std::get<3>(terminus.second)) ) {  // 5' terminus

					summary = ::bliss::debruijn::chain::make_summarized_chain(
							terminus.first, terminus.second, terminus.first, terminus.second);
					// set the in-edge frequency.
					std::get<5>(summary) = (*it).second.get_in_edge_total_frequency();  // looking for in edge freq, set as in freq.
					std::get<6>(summary) = (*it).second.get_out_edge_total_frequency();  // looking for in edge freq, set as in freq.
					results.emplace_back(summary);
				} else if (::bliss::debruijn::is_chain_terminal(std::get<2>(terminus.second))) {  // 5' terminus
					if (terminus.first < std::get<1>(terminus.second).reverse_complement() ) {   // 5' is chain rep
						summary = ::bliss::debruijn::chain::make_summarized_chain(
								terminus.first, terminus.second, static_cast<int>(0));
						// set the in-edge frequency.
						std::get<5>(summary) = (*it).second.get_in_edge_total_frequency();  // looking for in edge freq, set as in freq.
						results.emplace_back(summary);
					} else if (terminus.first > std::get<1>(terminus.second).reverse_complement() ) {  // rc(3') is chain rep
						summary = ::bliss::debruijn::chain::make_summarized_chain(
								static_cast<int>(0), terminus.first.reverse_complement(),
								::bliss::debruijn::transform::reverse_complement(terminus.second));
						std::get<6>(summary) = (*it).second.get_in_edge_total_frequency();  // looking for in edge freq, set as out frequency
						results.emplace_back(summary);			
					} else {  // 5' and rc(3') are the same.  this is a mobius strip?
								//  note that 5' and 3' kmers are on same strand, and we do not allow palindromes as chains, so dist can't be 0.
								// mobius strips in a chain here means there one of the node is also a branch.  not possible.
						throw std::logic_error("ERROR: we may have a MOBIUS STRIP, or a palindrome");
					}
				} else if (::bliss::debruijn::is_chain_terminal(std::get<3>(terminus.second))) {  // 3' terminus
					if (std::get<0>(terminus.second) < terminus.first.reverse_complement() ) {   // 5' is chain rep
						summary = ::bliss::debruijn::chain::make_summarized_chain(
								static_cast<int>(0), terminus.first, terminus.second);
						std::get<6>(summary) = (*it).second.get_out_edge_total_frequency();  // looking for out edge freq, set as out frequency
						results.emplace_back(summary);
					} else if (std::get<0>(terminus.second) > terminus.first.reverse_complement() ) {  // rc(3') is chain rep
						summary = ::bliss::debruijn::chain::make_summarized_chain(
								terminus.first.reverse_complement(),
								::bliss::debruijn::transform::reverse_complement(terminus.second),
								static_cast<int>(0));
						std::get<5>(summary) = (*it).second.get_out_edge_total_frequency();  // looking for out edge freq, set as in freq.
						results.emplace_back(summary);
					} else {  // 5' and rc(3') are the same.  this is a mobius strip?
								//  note that 5' and 3' kmers are on same strand, and we do not allow palindromes as chains, so dist can't be 0.
								// mobius strips in a chain here means there one of the node is also a branch.  not possible.
						throw std::logic_error("ERROR: we may have a MOBIUS STRIP, or a palindrome");
					}
				} else {
					throw std::logic_error("ERROR:  NOT A TERMINUS.");
				}
			}
		}
		// evenly redistribute.
		mxx::distribute_inplace(results, comm);

		// participate only if there is at least one entry.
		bool has_data = (results.size() > 0);
		bool all_has_data = ::mxx::all_of(has_data, comm);
		mxx::comm subcomm = all_has_data ? comm.copy() : comm.split(has_data);

		// distribute, sort by chain rep, and merge (max should do it provided that no edge is marked as 0.)
		if (has_data) {
			// now sort by chain representatives (position 1.)
			mxx::sort(results.begin(), results.end(),
				[](::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> const & lhs,
				  ::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> const & rhs){
					// compare position 1 kmer, which are chain reps.   then use position 3 kmer to further sort.
					return (std::get<1>(lhs) < std::get<1>(rhs)) || 
							((std::get<1>(lhs) == std::get<1>(rhs)) && (std::get<3>(lhs) < std::get<3>(rhs)));
			}, subcomm);
			
			// termini comes in pairs, except for isolated kmers.  need to shift left to if split so we can merge.
			::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> first = results.front();
			::bliss::debruijn::chain::summarized_chain<KmerType, typename DBG::count_type> last = results.back();
			first = mxx::left_shift(first, subcomm);
			last = mxx::right_shift(last, subcomm);

			// if the last in results is same as first in next proc, then remove it and process on the next proc
			size_t max_pos = results.size();  // in here because we have at least 1.
			if (std::get<1>(results.back()) == std::get<1>(first)) {
				--max_pos;
			}  
			size_t insert_at = 0;
			// now begin scanning through. and merge. elements could be single node, or paired.
			// if the first in results is same as last in the prev proc, then merge with the first element and increment position.
			for (size_t pos = 0; pos < max_pos; ++pos) {
				// check to see if we have a unit or isolated chain node.
				if (std::get<4>(results[pos]) == 0) {
					results[insert_at] = results[pos];
					++insert_at;
				} else if (std::get<1>(results[pos]) == std::get<1>(last)) {
					// two termini of the same chain, merge.
					results[insert_at] = ::bliss::debruijn::chain::merge_summarized_chains(results[pos], last);
					++insert_at;
				}  // else not same and not unit /isolated, must be a new pair
				last = results[pos];
			}
			results.erase(results.begin() + insert_at, results.end());

		}
		return results;
	}

	};




} // ns: graph

} // ns: debruijn

} // ns: bliss

#endif // DEBRUIJN_CHAIN_GRAPH_HPP_
