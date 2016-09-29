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
 * debruijn_graph_map.hpp
 *
 *  Created on: Aug 6, 2015
 *      Author: yongchao
 *      Author: tony pan
 */

#ifndef DEBRUIJN_GRAPH_MAP_HPP_
#define DEBRUIJN_GRAPH_MAP_HPP_

#include "bliss-config.hpp"

#include <utility> 			  // for std::pair

#include "containers/distributed_densehash_map.hpp"

#include "debruijn/debruijn_graph_node.hpp"	//node trait data structure storing the linkage information to the node

#include "utils/benchmark_utils.hpp"  // for timing.

#include "utils/logging.h"


namespace bliss {
namespace debruijn {

namespace graph {
/**
 * de bruijn map.  essentially a reduction map, but with slight differences.
 */
template<typename Kmer, typename Edge,
template <typename> class MapParams,
class Alloc = ::std::allocator< ::std::pair<const Kmer, Edge> >
>
class debruijn_graph_map :
public ::dsc::densehash_map<Kmer, Edge, MapParams,
::bliss::kmer::hash::sparsehash::special_keys<Kmer, true>,
 Alloc> {
	using Base = ::dsc::densehash_map<Kmer, Edge, MapParams,
			::bliss::kmer::hash::sparsehash::special_keys<Kmer, true>,
			 Alloc>;

 public:
	using local_container_type = typename Base::local_container_type;

	// std::unordered_multimap public members.
	using key_type              = typename local_container_type::key_type;
	using mapped_type           = typename local_container_type::mapped_type;
	using value_type            = typename local_container_type::value_type;
	using hasher                = typename local_container_type::hasher;
	using key_equal             = typename local_container_type::key_equal;
	using allocator_type        = typename local_container_type::allocator_type;
	using reference             = typename local_container_type::reference;
	using const_reference       = typename local_container_type::const_reference;
	using pointer               = typename local_container_type::pointer;
	using const_pointer         = typename local_container_type::const_pointer;
	using iterator              = typename local_container_type::iterator;
	using const_iterator        = typename local_container_type::const_iterator;
	using size_type             = typename local_container_type::size_type;
	using difference_type       = typename local_container_type::difference_type;

 protected:
	Edge dummy;

 public:

	/**
	 * @brief insert new elements in the distributed unordered_multimap.
	 * @param first
	 * @param last
	 */
	template <class InputIterator>
	size_t local_insert(InputIterator first, InputIterator last) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			auto result = this->c.insert(::std::make_pair((*it).first, Edge()));   // TODO: reduce number of allocations.
			// failed insertion - means an entry is already there, so reduce
			result.first->second.update((*it).second);
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}


	/**
	 * @brief insert new elements in the distributed unordered_multimap.
	 * @param first
	 * @param last
	 */
	template <class InputIterator, class Predicate>
	size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if (pred(*it)) {

				auto result = this->c.insert(::std::make_pair(it->first, Edge()));   // TODO: reduce number of allocations.
				// failed insertion - means an entry is already there, so reduce
				result.first->second.update(it->second);
			}
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;

	}

	debruijn_graph_map(const mxx::comm& _comm) : Base(_comm) {/*do nothing*/}

	virtual ~debruijn_graph_map() {/*do nothing*/};

	/*transform function*/

	/**
	 * @brief insert new elements in the distributed unordered_multimap.
	 * @param first
	 * @param last
	 */
	template <typename InputEdgeType, typename Predicate = ::fsc::TruePredicate>
	size_t insert(std::vector<::std::pair<Kmer, InputEdgeType> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
		// even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
		BL_BENCH_INIT(insert);


		if (::dsc::empty(input, this->comm)) {
			BL_BENCH_REPORT_MPI_NAMED(insert, "debruijnmap:insert", this->comm);
			return 0;
		}

		BL_BENCH_START(insert);
		this->transform_input(input);
		BL_BENCH_END(insert, "transform_input", input.size());

		// communication part
		if (this->comm.size() > 1) {
			BL_BENCH_START(insert);
			::std::vector<size_t> recv_counts =
					::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm);
			BLISS_UNUSED(recv_counts);
			BL_BENCH_END(insert, "dist_data", input.size());
		}

		BL_BENCH_START(insert);
		// local compute part.  called by the communicator.
		size_t count = 0;
		if (!::std::is_same<Predicate, ::fsc::TruePredicate>::value)
			count = this->local_insert(input.begin(), input.end(), pred);
		else
			count = this->local_insert(input.begin(), input.end());
		BL_BENCH_END(insert, "local_insert", this->c.size());

		BL_BENCH_REPORT_MPI_NAMED(insert, "debruijnmap:insert", this->comm);


		return count;
	}
};

template<typename Kmer >
using simple_hash_debruijn_graph_map = ::bliss::debruijn::graph::debruijn_graph_map<Kmer,
		::bliss::debruijn::graph::compact_multi_biedge<typename Kmer::KmerAlphabet, bool>,
		 ::bliss::debruijn::CanonicalDeBruijnHashMapParams>;

template<typename Kmer, typename Count = uint32_t >
using count_hash_debruijn_graph_map = ::bliss::debruijn::graph::debruijn_graph_map<Kmer,
		::bliss::debruijn::graph::compact_multi_biedge<typename Kmer::KmerAlphabet, Count>,
		 ::bliss::debruijn::CanonicalDeBruijnHashMapParams>;

} // namespace graph
}/*debruijn*/
}/*bliss*/

#endif /* DEBRUIJN_GRAPH_MAP_HPP_ */
