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
#include <iterator>

#include "containers/distributed_densehash_map.hpp"

// #include "debruijn/debruijn_graph_node.hpp"	//node trait data structure storing the linkage information to the node

// #include "utils/benchmark_utils.hpp"  // for timing.
// #include "utils/function_traits.hpp"
// #include "utils/filter_utils.hpp"

// #include "utils/logging.h"


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

	// wrapper for special keys for k2 counter.
	template <typename Key>
	struct k2_special_keys {
		static_assert(::std::is_same<Key, ::std::pair<key_type, typename Edge::EdgeInputType> >::value, "support only <kmer, edge_input_type> pairs");
		
		using k_special_keys = ::bliss::kmer::hash::sparsehash::special_keys<Kmer, true>;

		inline Key generate(uint8_t id = 0) {
			return std::make_pair(k_special_keys().generate(id), typename Edge::EdgeInputType());
		}

		inline Key invert(Key const &x) {
			return std::make_pair(k_special_keys().invert(x.first), x.second);;
		}

		inline Key get_splitter() {
			return std::make_pair(k_special_keys().get_splitter(), typename Edge::EdgeInputType());
		}

		static constexpr bool need_to_split = k_special_keys::need_to_split;
	};
	// wrapper for hash function for ::std::pair<key_type, typename Edge::EdgeInputType>


	/**
	 * @brief Kmer specialization for MurmurHash.  generated hash is 128 bit.
	 *
	 * TODO: move KMER type template param to operator.
	 * TODO: change h to member variable.
	 */
	template <typename Key>
	class k2_murmur {

	protected:
		static constexpr unsigned int nBytes = sizeof(Key);
		uint32_t seed;

	public:
		static constexpr uint8_t batch_size = 1;

		static const unsigned int default_init_value = 24U;  // allow 16M processors.  but it's ignored here.

		k2_murmur(const unsigned int prefix_bits = default_init_value, uint32_t const & _seed = ((1 << 17) - 1) ) : seed(_seed) {};

		inline uint64_t operator()(const Key & k) const	{
			// produces 128 bit hash.
			uint64_t h[2];
			// let compiler optimize out all except one of these.
			if (sizeof(void*) == 8)
				MurmurHash3_x64_128(reinterpret_cast<unsigned char const *>(&k), nBytes, seed, h);
			else if (sizeof(void*) == 4)
				MurmurHash3_x86_128(reinterpret_cast<unsigned char const *>(&k), nBytes, seed, h);
			else
				throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

			// use the upper 64 bits.
			return h[0];
		}

      };

 public:
	// for thresholded insert.  default frequency type is uint16_t, which should be sufficient for most use cases but does use more memory.
	using FreqType = 
		typename std::conditional<std::is_same<typename Edge::CountType, bool>::value,
			uint16_t, typename Edge::CountType>::type;
	// step 1:  build k2mer counter. only need it for aggregation. use murmurhash.
	using LocalK2CountMapType = ::fsc::densehash_map<
		::std::pair<key_type, typename Edge::EdgeInputType>, FreqType,
		k2_special_keys<::std::pair<key_type, typename Edge::EdgeInputType> >,
		::bliss::transform::identity,   // should not need to change relative to the transform used for the graph.
		k2_murmur<::std::pair<key_type, typename Edge::EdgeInputType> >,    	// need specialization
		::fsc::sparsehash::compare<::std::pair<key_type, typename Edge::EdgeInputType>, std::equal_to, ::bliss::transform::identity>,  // should not need to change
		::std::allocator<::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> >,
		k2_special_keys<::std::pair<key_type, typename Edge::EdgeInputType> >::need_to_split >;


	/**
	 * @brief insert new elements in the distributed unordered_multimap.
	 * @param first
	 * @param last
	 */
	template <class InputIterator,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			::std::pair<key_type, typename Edge::EdgeInputType>
		  >::value, int>::type = 1
		  >
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
	template <class InputIterator,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			key_type
		  >::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			this->c.insert(::std::make_pair(*it, Edge()));   // TODO: reduce number of allocations.
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}

	template <typename CT = typename Edge::CountType, class InputIterator,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType>
		  >::value && 
		  ::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			// try inserting a stub first
			auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
			// then reduce.
			result.first->second.update((*it).first.second);
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}

	// frequency version
	template <typename CT = typename Edge::CountType, class InputIterator,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType>
		  >::value && 
		  !::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
			// failed insertion - means an entry is already there, so reduce
			result.first->second.update((*it).first.second, (*it).second);
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}



	/**
	 * @brief insert new elements in the distributed unordered_multimap.
	 * @param first
	 * @param last
	 */
	template <class InputIterator, class Predicate,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			::std::pair<key_type, typename Edge::EdgeInputType>
		  >::value, int>::type = 1
	>
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

	/**
	 * @brief insert new elements in the distributed unordered_multimap, with default constructed mapped data. useful as pre-filter for updating.
	 * @param first
	 * @param last
	 */
	template <class InputIterator, class Predicate,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			key_type
		  >::value, int>::type = 1
	>
	size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if (pred(*it)) {

				this->c.insert(::std::make_pair(*it, Edge()));   // TODO: reduce number of allocations.
			}
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;

	}


	template <typename CT = typename Edge::CountType, class InputIterator, class Predicate,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType>
		  >::value && 
		  ::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if (pred(*it)) {
				// try inserting a stub first
				auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
				// then reduce.
				result.first->second.update((*it).first.second);
			}
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}

	// frequency version
	template <typename CT = typename Edge::CountType, class InputIterator, class Predicate,
		typename ::std::enable_if<
		  ::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType>
		  >::value && 
		  !::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if (pred(*it)) {
				auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
				// failed insertion - means an entry is already there, so reduce
				result.first->second.update((*it).first.second, (*it).second);
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
	 * @details InputElemType is ::std::pair<Kmer, InputEdgeType>
	 * @param first
	 * @param last
	 */
	template <typename InputElemType, typename Predicate = ::bliss::filter::TruePredicate>
	size_t insert(std::vector<InputElemType >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
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

			::std::vector<size_t> recv_counts;
	    	std::vector<InputElemType > output;

			::imxx::distribute(input, this->key_to_rank, recv_counts, output, this->comm);
      		output.swap(input);
			BL_BENCH_END(insert, "dist_data", output.size());
		}

		BL_BENCH_START(insert);
		// this->c.resize(input.size() / 2);
		if (this->comm.rank() == 0)
		std::cout << "rank " << this->comm.rank() <<
		" BEFFORE input=" << input.size() << " size=" << this->local_size() << " buckets=" << this->c.bucket_count() << std::endl;

		// local compute part.  called by the communicator.
		size_t count = 0;
		if (!::std::is_same<Predicate, ::bliss::filter::TruePredicate>::value)
			count = this->local_insert(input.begin(), input.end(), pred);
		else
			count = this->local_insert(input.begin(), input.end());

		if (this->comm.rank() == 0)
		  std::cout << "rank " << this->comm.rank() <<
        " AFTER input=" << input.size() << " size=" << this->local_size() << " reported=" << count << " buckets=" << this->c.bucket_count() << std::endl;

		BL_BENCH_END(insert, "local_insert", this->local_size());

		BL_BENCH_REPORT_MPI_NAMED(insert, "debruijnmap:insert", this->comm);


		return count;
	}


	/**
	 * @brief insert new elements in the distributed unordered_multimap.
	 * @details InputElemType is ::std::pair<Kmer, InputEdgeType>
	 * @param first
	 * @param last
	 */
	template <typename IT, typename Predicate = ::bliss::filter::TruePredicate>
	size_t insert_incremental(IT start, IT endd, 
		size_t block_size,
		bool sorted_input = false, Predicate const & pred = Predicate()) {
		// even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
		BL_BENCH_INIT(insert);

		bool emp = (start == endd);
		if (::mxx::all_of(emp, this->comm)) {  // empty
			BL_BENCH_REPORT_MPI_NAMED(insert, "debruijnmap:insert_incr", this->comm);
			return 0;
		}

		BL_BENCH_START(insert);
		using InputElemType = typename ::std::iterator_traits<IT>::value_type;
		std::vector<InputElemType> input;
		input.reserve(block_size);
		std::vector<InputElemType> output;
		output.reserve(block_size);
		::fsc::back_emplace_iterator<std::vector<InputElemType> > input_emplace_iter(input);
		::fsc::back_emplace_iterator<std::vector<InputElemType> > output_emplace_iter(output);

		size_t src_count = 0, dest_count = 0;

		size_t i = 0;
		IT bstart = start;
		IT bend = start;
		for (i = 0; (i < block_size) && (bend != endd); ++i, ++bend) {};
		bool all_done = (i == 0);
		all_done = mxx::all_of(all_done, this->comm);

		::std::vector<size_t> recv_counts(this->comm.size());

		BL_BENCH_COLLECTIVE_END(insert, "init", block_size, this->comm);

		if (this->comm.rank() == 0)
			std::cout << "rank " << this->comm.rank() <<
			" BEFFORE size=" << this->local_size() << " buckets=" << this->c.bucket_count() << std::endl;

		BL_BENCH_LOOP_START(insert, 0);  // for parse and transform
		BL_BENCH_LOOP_START(insert, 1);  // for distribute
		BL_BENCH_LOOP_START(insert, 2);  // for insert
		BL_BENCH_LOOP_START(insert, 3);  // for cleanup

		size_t count = 0;
		while (! all_done) {
			
			if (this->comm.size() > 1) {
				BL_BENCH_LOOP_RESUME(insert, 0);
				this->transform_input(bstart, bend, input_emplace_iter);
				src_count += input.size();
				BL_BENCH_LOOP_PAUSE(insert, 0);

			// communication part
				BL_BENCH_LOOP_RESUME(insert, 1);
				::imxx::distribute(input, this->key_to_rank, recv_counts, output, this->comm);
				BL_BENCH_LOOP_PAUSE(insert, 1);
			} else {
				BL_BENCH_LOOP_RESUME(insert, 0);
				this->transform_input(bstart, bend, output_emplace_iter);
				src_count += output.size();
				BL_BENCH_LOOP_PAUSE(insert, 0);
			}
			dest_count += output.size();

			BL_BENCH_LOOP_RESUME(insert, 2);
			if (!::std::is_same<Predicate, ::bliss::filter::TruePredicate>::value)
				count += this->local_insert(output.begin(), output.end(), pred);
			else
				count += this->local_insert(output.begin(), output.end());
			BL_BENCH_LOOP_PAUSE(insert, 2);

			BL_BENCH_LOOP_RESUME(insert, 3);
			bstart = bend;
			for (i = 0; (i < block_size) && (bend != endd); ++i, ++bend) {};
			all_done = (i == 0);
			all_done = mxx::all_of(all_done, this->comm);
			input.clear();
			output.clear();
			BL_BENCH_LOOP_PAUSE(insert, 3);

		}

		BL_BENCH_LOOP_END(insert, 0, "transform_input", src_count);  // for init
		BL_BENCH_LOOP_END(insert, 1, "dist_data", dest_count);  // for parse
		BL_BENCH_LOOP_END(insert, 2, "local_insert", this->local_size());  // for insert
		BL_BENCH_LOOP_END(insert, 3, "cleanup", block_size);  // for insert

		if (this->comm.rank() == 0)
			std::cout << "rank " << this->comm.rank() <<
			" AFTER src=" << src_count << " dest=" << dest_count << " size=" << this->local_size() << " reported=" << count << " buckets=" << this->c.bucket_count() << std::endl;

		BL_BENCH_REPORT_MPI_NAMED(insert, "debruijnmap:insert_incr", this->comm);


		return count;
	}



	/**
	 * @brief insert new elements in the debruijn graph map, with frequency filtering
	 * 
	// canonicalizing k2mer and k-mer maintains symmetry.  k1mer does not.  the lack of symmetry can result in double counting, etc, and requires communication
	// features of this "build".  
	//			1. reduces computation load by first count k2-mers.  reduce from N (nonunique K2mers) down to M (unique canonical k2mers)
	//			2. avoid communication during k1mer filtering by making k1mer count construction and query local. 
	//				a. by partitioning to processors using central canonical k-mers
	//				b. by local sorting using central canonical k-mers, then counting k1mer only within the range of k2mers with same central k-mer
	//						(ensures no double counting or missing counts from canonicalization, because each k-mer's neighbors are counted separately from the other k-mers.)
	//						(also, separating in and out edges as well as maintaining uncanonicalized k2mer are no longer necessary.)
	//			3. by computing k2, k1, and k-mer counts, can filter by all three.  k1-mer filtering correspond to accurate k-mer count
	//				a. k2-mer counting uses <k, [edges]> as KEY. 
	//				b. k1-mer counting can sort or just scan using the edge array WITHIN a central k-mer range.
	//			4. maintains original k2mer count during accumulation of k1 and k-mer counting.  filtering is recorded but applied at the end.
	//			5. maintains original input and output tuple formats, thus supporting double vs single strand graphs, and exist/freq notation of edges.
	//							

	// given 2 k2mers v1 and v2 that share a k1mer edge e, the edge frequencies as calculated using v1 and v2 frequencies should be identical  
	//			
	//	user choose one of k2mer, k1mer, or kmer threshold for filter.  otherwise a heuristic about which to apply first is needed.
	//          NOTE THAT reduction thresholding, and compaction steps are local, without communication.
	//  		NOTE that user specified frequency is for canonicalized (should be this way else could get mismatched forward and backward k and k2).
	//				k2, k1, and k counts can be computed either canonicalized or uncanonicalized.
	//  broken up into 2 functions so the k2mer counting can work with multiple files.
	// 
	// ==== step 1.  build k2mer counter from received InputElemType (<kmer, edge_tuples>), key is <kmer, edge_tuples>
	// ==== step 2.  scan the k2mer counter (to_vector)  <<k, edge_tuples>, c>
	// ==== step 3.  sort by k. (MlogM)  - simpler, less memory, and multiple summaries.
	// ==== step 4.  do per-block sum
	// ==== step 4a.  accumulate k-mer count
	// ==== step 4b.  scan and accumulate k1-mer count
	// ==== step 5.  do per-block threshold
	// ==== step 5a.  mark k2mer in delete array with 4
	// ==== step 5b.  mark kmer in delete array with 8
	// ==== step 5c.  mark k1mer in delete array with 1 (out) and 2 (in).  in/out follow convention of k-mer construction.
	// ==== step 6.  use delete array to compact the k2mer counters.  heuristic here assumes thresholds independently applied.
					//  if > 2, remove.  if 1, zero in.  if 2 zero out.
	// ==== step 7.  build dbg:  direct insertion into local hash map.
	 * @param first
	 * @param last
	 */
	size_t compute_biedge_freqencies(std::vector<::std::pair<key_type, typename Edge::EdgeInputType> >& input, 
		LocalK2CountMapType & k2_counter) {

		static_assert(std::is_same<typename Edge::EdgeInputType, bliss::debruijn::biedge::compact_simple_biedge>::value, "currently only supporitng compact_simple_biedge as input type.");

		// even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
		BL_BENCH_INIT(k2freq);

		if (::dsc::empty(input, this->comm)) {
			BL_BENCH_REPORT_MPI_NAMED(k2freq, "debruijnmap:compute_k2freq", this->comm);
			return 0;
		}

		BL_BENCH_START(k2freq);
		this->transform_input(input);
		BL_BENCH_END(k2freq, "transform_input", input.size());

		// communication part
		if (this->comm.size() > 1) {
			BL_BENCH_START(k2freq);

			::std::vector<size_t> recv_counts;
	    	std::vector<::std::pair<key_type, typename Edge::EdgeInputType> > output;

			::imxx::distribute(input, this->key_to_rank, recv_counts, output, this->comm);
      		output.swap(input);
			BL_BENCH_END(k2freq, "dist_data", input.size());
		}

		BL_BENCH_START(k2freq);

		{  // scope to clean up k2_counter later

			// step 1:  build k2mer counter. only need it for aggregation. use murmurhash.
			dsc::sat_plus<FreqType> k2_sat_add;

			// insert into counter.  need to make this a saturating counter.
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> kv = 
				std::make_pair(std::make_pair(key_type(), typename Edge::EdgeInputType()), FreqType(1));
			for (auto it = input.begin(); it != input.end(); ++it) {
				kv.first = *it;
				auto result = k2_counter.insert(kv);
				if (!(result.second)) {
					result.first->second = k2_sat_add(result.first->second, 1);
				}
			}
		} // automatically cleans up.
		BL_BENCH_END(k2freq, "count_k2mers", k2_counter.size());

		BL_BENCH_REPORT_MPI_NAMED(k2freq, "debruijnmap:compute_k2freq", this->comm);


		return k2_counter.size();

	}

	template <typename IT,
		typename ::std::enable_if<::std::is_same<
			typename ::std::iterator_traits<IT>::value_type, 
			::std::pair<key_type, typename Edge::EdgeInputType> >::value,
			int>::type = 1>
	size_t compute_biedge_freqencies_incremental( IT start, IT endd, 
		LocalK2CountMapType & k2_counter, size_t block_size) {

		static_assert(std::is_same<typename Edge::EdgeInputType, bliss::debruijn::biedge::compact_simple_biedge>::value, "currently only supporitng compact_simple_biedge as input type.");

		// even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
		BL_BENCH_INIT(k2freq);

		bool emp = (start == endd);
		if (::mxx::all_of(emp, this->comm)) {  // empty
			BL_BENCH_REPORT_MPI_NAMED(k2freq, "debruijnmap:compute_k2freq_incr", this->comm);
			return 0;
		}

		{  // scope to clean up k2_counter later
			BL_BENCH_START(k2freq);
			using InputElemType = typename ::std::iterator_traits<IT>::value_type;
			std::vector<InputElemType> input;
			input.reserve(block_size);
			std::vector<InputElemType> output;
			output.reserve(block_size);
			::fsc::back_emplace_iterator<std::vector<InputElemType> > input_emplace_iter(input);
			::fsc::back_emplace_iterator<std::vector<InputElemType> > output_emplace_iter(output);

			size_t src_count = 0, dest_count = 0;

			size_t i = 0;
			IT bstart = start;
			IT bend = start;
			for (i = 0; (i < block_size) && (bend != endd); ++i, ++bend) {};
			bool all_done = (i == 0);
			all_done = mxx::all_of(all_done, this->comm);

			::std::vector<size_t> recv_counts(this->comm.size());

			// step 1:  build k2mer counter. only need it for aggregation. use murmurhash.
			::dsc::sat_plus<FreqType> k2_sat_add;

			BL_BENCH_COLLECTIVE_END(k2freq, "init", block_size, this->comm);

			if (this->comm.rank() == 0)
				std::cout << "rank " << this->comm.rank() <<
				" BEFFORE size=" << this->local_size() << " buckets=" << this->c.bucket_count() << std::endl;

			BL_BENCH_LOOP_START(k2freq, 0);  // for parse and transform
			BL_BENCH_LOOP_START(k2freq, 1);  // for distribute
			BL_BENCH_LOOP_START(k2freq, 2);  // for insert into k2mer counter
			BL_BENCH_LOOP_START(k2freq, 3);  // for cleanup

			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> kv = std::make_pair(InputElemType(), FreqType(1));
			while (! all_done) {
				
				if (this->comm.size() > 1) {
					BL_BENCH_LOOP_RESUME(k2freq, 0);
					this->transform_input(bstart, bend, input_emplace_iter);
					src_count += input.size();
					BL_BENCH_LOOP_PAUSE(k2freq, 0);

				// communication part
					BL_BENCH_LOOP_RESUME(k2freq, 1);
					::imxx::distribute(input, this->key_to_rank, recv_counts, output, this->comm);
					BL_BENCH_LOOP_PAUSE(k2freq, 1);
				} else {
					BL_BENCH_LOOP_RESUME(k2freq, 0);
					this->transform_input(bstart, bend, output_emplace_iter);
					src_count += output.size();
					BL_BENCH_LOOP_PAUSE(k2freq, 0);
				}
				dest_count += output.size();

				BL_BENCH_LOOP_RESUME(k2freq, 2);
				// insert into counter.  need to make this a saturating counter.
				for (auto it = output.begin(); it != output.end(); ++it) {
					kv.first = *it;
					auto result = k2_counter.insert(kv);
					if (!(result.second)) {
						result.first->second = k2_sat_add(result.first->second, 1);
					}
				}
				BL_BENCH_LOOP_PAUSE(k2freq, 2);

				BL_BENCH_LOOP_RESUME(k2freq, 3);
				bstart = bend;
				for (i = 0; (i < block_size) && (bend != endd); ++i, ++bend) {};
				all_done = (i == 0);
				all_done = mxx::all_of(all_done, this->comm);
				input.clear();
				output.clear();
				BL_BENCH_LOOP_PAUSE(k2freq, 3);

			}

			BL_BENCH_LOOP_END(k2freq, 0, "transform_input", src_count);  // for init
			BL_BENCH_LOOP_END(k2freq, 1, "dist_data", dest_count);  // for parse
			BL_BENCH_LOOP_END(k2freq, 2, "count_k2mer", k2_counter.size());  // for insert
			BL_BENCH_LOOP_END(k2freq, 3, "cleanup", block_size);  // for insert
		} // automatically cleans up.


		BL_BENCH_REPORT_MPI_NAMED(k2freq, "debruijnmap:compute_k2freq_incr", this->comm);


		return k2_counter.size();

	}

	// actual insertion using k2frequencies.
	template <typename Predicate = ::bliss::filter::TruePredicate>
	size_t local_insert_by_freqencies(
		LocalK2CountMapType & k2_counter,
		std::vector<FreqType> const & threshes,  // order is k_low, k_hi, k1_low, k1_hi, k2_low, k2_hi
		Predicate const & pred = Predicate()) {

		BL_BENCH_INIT(local_insert);

		// step 2: get the content of the k2 count map.
		BL_BENCH_START(local_insert);
		::std::vector<::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> > k2counts;
		k2_counter.to_vector(k2counts);
		k2_counter.reset(); // get rid of allocated space.
		BL_BENCH_END(local_insert, "to_vec", k2counts.size());

		// step 3: sort by kmer.  no transform of the kmer.
		BL_BENCH_START(local_insert);
		std::sort(k2counts.begin(), k2counts.end(),
			[](::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> const & x,
			   ::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> const & y) {
				   return x.first.first < y.first.first;
		});
		BL_BENCH_END(local_insert, "k2mer_sort", k2counts.size());

		{ // scope to clear filter_flag later.
			// step 4: do per block sum
			// step 5: do per block threshold (interleaved.)
			BL_BENCH_START(local_insert);
			std::vector<unsigned char> filter_flag(k2counts.size(), 0);
			auto flag_it = filter_flag.begin();
			{  // scope to clear local arrays later
				size_t k_f = 0;
				size_t k1in_f[4];  memset(k1in_f, 0, sizeof(size_t) << 2);   // ignore empty.
				size_t k1out_f[4];  memset(k1out_f, 0, sizeof(size_t) << 2);  // ignore empty.
				size_t k1cnt;
				key_type k;
				FreqType c;
				uint8_t edges;
				auto block_it = k2counts.begin();
				auto it = k2counts.begin();
				FreqType t0_lo = threshes[0],
					t0_hi = threshes[1],
					t1_lo = threshes[2],
					t1_hi = threshes[3],
					t2_lo = threshes[4],
					t2_hi = threshes[5];
				if (it != k2counts.end()) k = it->first.first;  // init
			
				for (; it != k2counts.end(); ++it) {
					if (it->first.first != k) {  // start of new block
						// accumulation for block complete. filter and record results now.
						for (; block_it != it; ++block_it, ++flag_it) {
							// step 5: do per block threshold
							c = block_it->second;
							
							// step 5a: filter k2mer 
							if ((c < t2_lo) || (c >= t2_hi)) *flag_it |= 4;
							// step 5b: filter kmer
							if ((k_f < t0_lo) || (k_f >= t0_hi)) *flag_it |= 8;
							// step 5c: filter k1mer
							edges = block_it->first.second.getData()[0];

							switch (edges & 0xF0) {
								case 0x10: k1cnt = k1in_f[0]; break;
								case 0x20: k1cnt = k1in_f[1]; break;
								case 0x40: k1cnt = k1in_f[2]; break;
								case 0x80: k1cnt = k1in_f[3]; break;
								default: k1cnt = t1_hi; break;
							}
							if ((k1cnt < t1_lo) || (k1cnt >= t1_hi)) *flag_it |= 2;  // in

							switch (edges & 0x0F) {
								case 0x01: k1cnt = k1out_f[0]; break;
								case 0x02: k1cnt = k1out_f[1]; break;
								case 0x04: k1cnt = k1out_f[2]; break;
								case 0x08: k1cnt = k1out_f[3]; break;
								default: k1cnt =  t1_hi; break;
							}
							if ((k1cnt < t1_lo) || (k1cnt >= t1_hi)) *flag_it |= 1;  // out

						}
						
						// reinitialize for next block
						k = block_it->first.first;
						k_f = 0;
						memset(k1in_f, 0, sizeof(size_t) << 2);
						memset(k1out_f, 0, sizeof(size_t) << 2);
					}

					// step 4: do per block sum
					c = it->second;
					// step 4a: accumulate the kmer counts. 
					k_f += c;   // accumulate the kmer count
					// step 4b: accumulate the k1mer in and out edge counts. 
					edges = it->first.second.getData()[0];
					switch (edges & 0xF0) {
						case 0x10: k1in_f[0] += c; break;
						case 0x20: k1in_f[1] += c; break;
						case 0x40: k1in_f[2] += c; break;
						case 0x80: k1in_f[3] += c; break;
						default:
							assert("ERROR in edge has a combination of characters");
							break;
					}
					switch (edges & 0x0F) {
						case 0x01: k1out_f[0] += c; break;
						case 0x02: k1out_f[1] += c; break;
						case 0x04: k1out_f[2] += c; break;
						case 0x08: k1out_f[3] += c; break;
						default:
							assert("ERROR out edge has a combination of characters");
							break;
					}
				}
				// do the final block
				for (; block_it != it; ++block_it, ++flag_it) {
					// step 5: do per block threshold
					c = block_it->second;
					
					// step 5a: filter k2mer 
					if ((c < t2_lo) || (c >= t2_hi)) *flag_it |= 4;
					// step 5b: filter kmer
					if ((k_f < t0_lo) || (k_f >= t0_hi)) *flag_it |= 8;
					// step 5c: filter k1mer
					edges = block_it->first.second.getData()[0];

					switch (edges & 0xF0) {
						case 0x10: k1cnt = k1in_f[0]; break;
						case 0x20: k1cnt = k1in_f[1]; break;
						case 0x40: k1cnt = k1in_f[2]; break;
						case 0x80: k1cnt = k1in_f[3]; break;
						default: k1cnt = t1_hi; break;
					}
					if ((k1cnt < t1_lo) || (k1cnt >= t1_hi)) *flag_it |= 2;  // in

					switch (edges & 0x0F) {
						case 0x01: k1cnt = k1out_f[0]; break;
						case 0x02: k1cnt = k1out_f[1]; break;
						case 0x04: k1cnt = k1out_f[2]; break;
						case 0x08: k1cnt = k1out_f[3]; break;
						default: k1cnt =  t1_hi; break;
					}
					if ((k1cnt < t1_lo) || (k1cnt >= t1_hi)) *flag_it |= 1;  // out
				}

			}  // end block for accumulation and thresholding.
			BL_BENCH_END(local_insert, "count_thresh", k2counts.size());


			BL_BENCH_START(local_insert);
			{ // scope to clear stats array
				// step 6: modify and compact.
				auto it = k2counts.data();
				auto write_it = k2counts.begin();
				flag_it = filter_flag.begin();
				uint8_t flag = 0;
				size_t i = 0, imax = k2counts.size();
				size_t k_del[16];  memset(k_del, 0, sizeof(size_t) * 16);
				for (; i < imax; ++it, ++flag_it, ++i) {
					flag = *flag_it;
					++k_del[flag];

					if (flag < 3) {
						memmove(&(*write_it), it, sizeof(::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, FreqType> ));
						switch (flag) {
							case 1:
								write_it->first.second.getDataRef()[0] &= 0xF0;
								break;  // remove out and keep
							case 2: 
								write_it->first.second.getDataRef()[0] &= 0x0F;
								break;  // remove in and keep
							default: break; 
						}
						++write_it;
					} // else discard
				}
				if (this->comm.rank() == 0) {
					std::cout << "rank " << this->comm.rank() <<
					" BEFFORE filtering unique k2mers=" << k2counts.size() << " filter count:" << std::endl;
					for (int j = 0; j < 16; ++j) {
						std::cout << "\t[" << j << " : " << k_del[j] << "]" << std::endl;
					}
				}
				k2counts.erase(write_it, k2counts.end());
			} // end scope to delete stats array
			BL_BENCH_END(local_insert, "filter", k2counts.size());
		} // end scope to delete filter_flag		

		// step 7: insert k2counts now
		BL_BENCH_COLLECTIVE_START(local_insert, "local_insert", this->comm);
		if (this->comm.rank() == 0)
			std::cout << "rank " << this->comm.rank() <<
			" BEFFORE filtered unique k2mers =" << k2counts.size() << 
			" size=" << this->local_size() << " buckets=" << this->c.bucket_count() << std::endl;

		// local compute part.  called by the communicator.
		size_t count = 0;
		if (!::std::is_same<Predicate, ::bliss::filter::TruePredicate>::value)
			count = this->local_insert(k2counts.begin(), k2counts.end(), pred);
		else
			count = this->local_insert(k2counts.begin(), k2counts.end());

		if (this->comm.rank() == 0)
		  std::cout << "rank " << this->comm.rank() <<
        	" AFTER size=" << this->local_size() << " buckets=" << this->c.bucket_count() << std::endl;

		BL_BENCH_END(local_insert, "local_insert", this->c.size());

		BL_BENCH_REPORT_MPI_NAMED(local_insert, "debruijnmap:local_insert_freq", this->comm);


		return count;
	}

};

template<typename Kmer>
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
