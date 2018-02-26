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
 * debruijn_graph.hpp
 *
 *  Created on: Oct 2, 2016
 *      Author: tony pan
 */

#ifndef DEBRUIJN_GRAPH_HPP_
#define DEBRUIJN_GRAPH_HPP_

// #include "index/kmer_index.hpp"
#include "utils/function_traits.hpp"

#include "debruijn/debruijn_graph_map.hpp"
#include "debruijn/debruijn_graph_filters.hpp"


namespace bliss
{

namespace debruijn
{

namespace graph
{


	/// a class representing a debruijn graph
	template <typename KmerType, typename FreqValueType>
	class debruijn_graph {

	public:
		template <typename KMER>
		using map_params_template = ::bliss::debruijn::CanonicalDeBruijnHashMapParams<KMER>;
		using map_params_type =
				map_params_template<KmerType>;

		// full debruijn graph map - with or without frequency
		using map_type =
				::bliss::debruijn::graph::debruijn_graph_map<KmerType,
				 ::bliss::debruijn::graph::compact_multi_biedge<typename KmerType::KmerAlphabet, FreqValueType>,
				  map_params_template>;

		using key_type =  typename map_type::key_type;
		using kmer_type =  typename map_type::key_type;
		using edge_type =  typename map_type::mapped_type;
		using value_type = typename map_type::value_type;
		using mutable_value_type = ::std::pair<kmer_type, edge_type>;
		using const_iterator = typename map_type::const_iterator;

	protected:
		/// distributed map fo storage.
		map_type map;

		mxx::comm comm;

	public:

		debruijn_graph(mxx::comm const & _comm) : map(_comm), comm(_comm.copy()) {};

		virtual ~debruijn_graph() {};

		/// get internal map
		map_type & get_map() { return map; }
		map_type const & get_map() const { return map; }

		/// insert kmers with edges.
		template <typename T>
		size_t insert(std::vector<T> & nodes,
				bool is_local = false) {
			// insert into underlying map
			static_assert(std::is_convertible<T, KmerType>::value ||
					std::is_convertible<T, ::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> >::value,
					"only support input type of Kmers, or Kmer+biedge pairs");

			// do not local reserve since repeats may be high and the map is a reduction map.
			//map.local_reserve(node.size());
			// take the realloc hit.

			bool all_local = mxx::all_of(is_local, comm);
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
		template <typename T, typename Predicate>
		size_t insert(std::vector<T> & nodes,
				Predicate const & pred,
				bool is_local = false) {

			static_assert(std::is_convertible<T, KmerType>::value ||
					std::is_convertible<T, ::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> >::value,
					"only support input type of Kmers, or Kmer+biedge pairs");

			// doing it this way to ensure that we get the intersection and not a complement
			bool all_local = mxx::all_of(is_local, comm);
			if (all_local) {
				return map.local_insert(nodes.begin(), nodes.end(), pred);
			} else {
				return map.insert(nodes, false, pred);  // includes distribution.
			}
		}


		/// insert kmers with edges.
		template <typename IT>
		size_t insert_incremental(IT start, IT endd, 
				size_t block_size = std::numeric_limits<size_t>::max(),
				bool is_local = false) {
			// insert into underlying map
			static_assert(std::is_convertible<typename ::std::iterator_traits<IT>::value_type, KmerType>::value ||
					std::is_convertible<typename ::std::iterator_traits<IT>::value_type, ::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> >::value,
					"only support input type of Kmers, or Kmer+biedge pairs");

			// do not local reserve since repeats may be high and the map is a reduction map.
			//map.local_reserve(node.size());
			// take the realloc hit.

			bool all_local = mxx::all_of(is_local, comm);
			if (all_local) {
				return map.local_insert(start, endd);
			} else {
				return map.insert_incremental(start, endd, block_size);  // includes distribution.
			}
		}

		/**
		 * @brief	insert the selected node into index.  predicate returns true means the element should be inserted.
		 * @details	example of predicate:  using another index along with a subpredicate, evaluate 1. if the element is in that index, and 2. is the subpredicate satisfied.
		 * 					we can insert e.g. just branch nodes that are in the graph, into the node.
		 */
		template <typename IT, typename Predicate>
		size_t insert_incremental(IT start, IT endd,
				Predicate const & pred,
				size_t block_size = std::numeric_limits<size_t>::max(),
				bool is_local = false) {

			static_assert(std::is_convertible<typename ::std::iterator_traits<IT>::value_type, KmerType>::value ||
					std::is_convertible<typename ::std::iterator_traits<IT>::value_type, ::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> >::value,
					"only support input type of Kmers, or Kmer+biedge pairs");

			// doing it this way to ensure that we get the intersection and not a complement
			bool all_local = mxx::all_of(is_local, comm);
			if (all_local) {
				return map.local_insert(start, endd, pred);
			} else {
				return map.insert_incremental(start, endd, block_size, false, pred);  // includes distribution.
			}
		}

		//==========
		//  new threshold version:  use the underlying map directly please.  2 calls:  compute_biedge_frequencies followed by local_insert_by_frequencies
		//==========

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
		std::vector<const_iterator> find_if_iterators(Predicate const &pred) const {
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
			return map.update(val, sorted, updater);
		}

		// no compact function


	   typename map_type::const_iterator cbegin() const
	   {
		 return map.get_local_container().cbegin();
	   }

	   typename map_type::const_iterator cend() const {
		 return map.get_local_container().cend();
	   }

	   size_t size() const {
		 return map.size();
	   }

	   size_t local_size() const {
		 return map.local_size();
	   }


	   //======  convenience functions

	   /// return local branch nodes
	   std::vector<mutable_value_type> get_branch_nodes() const {
		   return this->find_if(::bliss::debruijn::filter::graph::IsBranchPoint());
	   }
	   /// return local branch node's kmer
	   std::vector<KmerType> get_branch_node_kmers() const {
		   return this->find_if_and_transform(::bliss::debruijn::filter::graph::IsBranchPoint(),
				   [](mutable_value_type const & x){
			   return x.first;
		   });
	   }


	};





} // ns: graph

} // ns: debruijn

} // ns: bliss

#endif // DEBRUIJN_GRAPH_HPP_
