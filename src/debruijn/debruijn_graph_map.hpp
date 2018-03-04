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
#include "kmer_traits.hpp"

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
	using mutable_value_type    = std::pair<key_type, mapped_type>;
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
	 * @brief K2mer specialization for MurmurHash.  generated hash is 128 bit.
	 *  does hash combining, using boost version.  from 
	 * https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
	 * TODO: move KMER type template param to operator.
	 * TODO: change h to member variable.
	 */
	template <typename KT, typename ET>
	class k2_murmur {

	protected:
		static constexpr unsigned int KTBytes = sizeof(KT);
		static constexpr unsigned int ETBytes = sizeof(ET);
		uint32_t seed;

	public:
		static constexpr uint8_t batch_size = 1;

		static const unsigned int default_init_value = 24U;  // allow 16M processors.  but it's ignored here.

		k2_murmur(const unsigned int prefix_bits = default_init_value, uint32_t const & _seed = ((1 << 17) - 1) ) : seed(_seed) {};

		inline uint64_t operator()(std::pair<KT, ET> const & k) const	{
			// produces 128 bit hash.
			uint64_t h[2], f[2];
			// let compiler optimize out all except one of these.
			if (sizeof(void*) == 8) {
				MurmurHash3_x64_128(reinterpret_cast<unsigned char const *>(&(k.first)), KTBytes, seed, h);
				MurmurHash3_x64_128(reinterpret_cast<unsigned char const *>(&(k.second)), ETBytes, seed, f);
			} else if (sizeof(void*) == 4) {
				MurmurHash3_x86_128(reinterpret_cast<unsigned char const *>(&(k.first)), KTBytes, seed, h);
				MurmurHash3_x86_128(reinterpret_cast<unsigned char const *>(&(k.second)), ETBytes, seed, f);
			} else
				throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

			// use the upper 64 bits.
			// mixing hash functions.  uses Phttps://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
			return h[0] ^ (f[0] + 0x517cc1b727220a95 +(h[0] << 6) + (h[0] >> 2));
		}

      };

 public:
	// for thresholded insert.  default frequency type is uint8_t, which should be sufficient for most use cases but does use more memory.
	using FreqType = 
		typename std::conditional<std::is_same<typename Edge::CountType, bool>::value,
			uint8_t, typename Edge::CountType>::type;
	// step 1:  build k2mer counter. only need it for aggregation. use murmurhash.

	using K2FreqType = uint32_t;
	using LocalK2CountMapType = ::fsc::densehash_map<
		::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType,
		k2_special_keys<::std::pair<key_type, typename Edge::EdgeInputType> >,
		::bliss::transform::identity,   // should not need to change relative to the transform used for the graph.
		k2_murmur<key_type, typename Edge::EdgeInputType>,    	// need specialization
		::fsc::sparsehash::compare<::std::pair<key_type, typename Edge::EdgeInputType>, std::equal_to, ::bliss::transform::identity>,  // should not need to change
		::std::allocator<::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> >,
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
		  (::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType>
		  >::value) && 
		  ::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			// try inserting a stub first
			if ((*it).second > 0) {
				auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
				// then reduce.
				// cap to max of FeqType
				result.first->second.update((*it).first.second);
			}
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}

	// frequency version.  cap to FreqType max.
	template <typename CT = typename Edge::CountType, class InputIterator,
		typename ::std::enable_if<
		  (::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType>
		  >::value)  && 
		  !::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if ((*it).second > 0) {
				auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
				// failed insertion - means an entry is already there, so reduce
				result.first->second.update((*it).first.second, 
					(std::numeric_limits<FreqType>::max() >= std::numeric_limits<K2FreqType>::max()) ? 
						static_cast<FreqType>((*it).second) :
						static_cast<FreqType>(std::min(static_cast<K2FreqType>(std::numeric_limits<FreqType>::max()), (*it).second))
					);
			}
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
				result.first->second.update((*it).second);
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
		  (::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType>
		  >::value)  && 
		  ::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if (pred(*it) && ((*it).second > 0)) {
				// try inserting a stub first
				auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
				// then reduce.
				result.first->second.update((*it).first.second);
			}
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}

	// frequency version.  cap to FreqType max.
	template <typename CT = typename Edge::CountType, class InputIterator, class Predicate,
		typename ::std::enable_if<
		  (::std::is_same<
			typename ::std::iterator_traits<InputIterator>::value_type,
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType>
		  >::value) && 
		  !::std::is_same<CT, bool>::value, int>::type = 1
		  >
	size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
		size_t before = this->c.size();

		// first/last could be a lot, especially for highly repeated reads.
		//this->local_reserve(before + ::std::distance(first, last));

		for (auto it = first; it != last; ++it) {
			if (pred(*it) && ((*it).second > 0)) {
				auto result = this->c.insert(::std::make_pair((*it).first.first, Edge()));   // TODO: reduce number of allocations.
				// failed insertion - means an entry is already there, so reduce
				result.first->second.update((*it).first.second, 
					(std::numeric_limits<FreqType>::max() >= std::numeric_limits<K2FreqType>::max()) ? 
					static_cast<FreqType>((*it).second) :
					static_cast<FreqType>(std::min(static_cast<K2FreqType>(std::numeric_limits<FreqType>::max()), (*it).second))
				);
			}
		}

		if (this->c.size() != before) this->local_changed = true;

		return this->c.size() - before;
	}



	debruijn_graph_map(const mxx::comm& _comm) : Base(_comm) {/*do nothing*/}

	virtual ~debruijn_graph_map() {/*do nothing*/};


	/// ====== edge and node erase and find operations.

	protected:
		// given multi-biedges, compute the reverse edges.
		std::vector<std::pair<key_type, typename Edge::EdgeInputType> >
		get_remote_edges(std::vector<mutable_value_type> & targets) {

			// ===  then disconnect the edges
			// ===== extract edges (in reverse form, still 5' to 3')
			std::vector<std::pair<key_type, typename Edge::EdgeInputType> > edges_to_remove;
			std::vector<std::pair<key_type, typename Edge::CountType> > neighbors;
			typename Edge::EdgeInputType out;
			typename Edge::EdgeInputType in;
			
			for (auto x : targets) {
				// process in edges
				out.setCharsAtPos(
					Edge::EdgeInputType::KmerAlphabet::FROM_ASCII[
						key_type::KmerAlphabet::TO_ASCII[x.first.getCharsAtPos(0, 1)]],
					0, 1);
				neighbors.clear();
				x.second.get_in_neighbors(x.first, neighbors);
				for (auto n : neighbors) {
					edges_to_remove.emplace_back(n.first, out);  // reverse version of 
				}
				// process out edges
				in.setCharsAtPos(
					Edge::EdgeInputType::KmerAlphabet::FROM_ASCII[
						key_type::KmerAlphabet::TO_ASCII[x.first.getCharsAtPos(key_type::size - 1, 1)]],
					1, 1);
				neighbors.clear();
				x.second.get_out_neighbors(x.first, neighbors);
				for (auto n : neighbors) {
					edges_to_remove.emplace_back(n.first, in);  // reverse version of 
				}
			}

			return edges_to_remove;
		}

		// given simple-biedges, compute the reverse edges.
		std::vector<std::pair<key_type, typename Edge::EdgeInputType> >
		get_remote_edges(std::vector<std::pair<key_type, typename Edge::EdgeInputType> > & targets) {

			// ===  then disconnect the edges
			// ===== extract edges (in reverse form, still 5' to 3')
			std::vector<std::pair<key_type, typename Edge::EdgeInputType> > edges_to_remove;
			edges_to_remove.reserve(targets.size());
			typename Edge::EdgeInputType out, in;
			uint8_t ch;
			
			for (auto x : targets) {
				// process in edges
				ch = x.second.getCharsAtPos(1,1);
				if ((ch > 0) && ((ch & (ch-1)) == 0)) {  // has 1 and only 1 character at the in edge position (DNA16)
					out.setCharsAtPos(
						Edge::EdgeInputType::KmerAlphabet::FROM_ASCII[
							key_type::KmerAlphabet::TO_ASCII[x.first.getCharsAtPos(0, 1)]],
						0, 1);
					key_type y = x.first;
					y.nextReverseFromChar(
						key_type::KmerAlphabet::FROM_ASCII[ Edge::EdgeInputType::KmerAlphabet::TO_ASCII[ch]]);
					edges_to_remove.emplace_back(y, out);
				}
				ch = x.second.getCharsAtPos(0,1);
				if ((ch > 0) && ((ch & (ch-1)) == 0)) {  // has 1 and only 1 character at the in edge position (DNA16)
					in.setCharsAtPos(
						Edge::EdgeInputType::KmerAlphabet::FROM_ASCII[
							key_type::KmerAlphabet::TO_ASCII[x.first.getCharsAtPos(key_type::size - 1, 1)]],
						1, 1);
					key_type y = x.first;
					y.nextFromChar(
						key_type::KmerAlphabet::FROM_ASCII[ Edge::EdgeInputType::KmerAlphabet::TO_ASCII[ch]]);
					edges_to_remove.emplace_back(y, in);
				}
			}

			return edges_to_remove;
		}

		// erase simple biedges from map.  distributed.
		inline void erase_simple_biedges(std::vector<std::pair<key_type, typename Edge::EdgeInputType> > & query) {
			this->update(query, false, 
				[](mapped_type & curr, typename Edge::EdgeInputType const & val){
					size_t before = curr.get_in_edge_count() + curr.get_out_edge_count();
					curr.erase(val);
					return curr.get_in_edge_count() + curr.get_out_edge_count() - before;
			});
		}
		// erase multi biedge from map.  distributed.
		inline void erase_multi_biedges(std::vector<mutable_value_type > & query) {
			this->update(query, false, 
				[](mapped_type & curr, mapped_type const & val){
					size_t before = curr.get_in_edge_count() + curr.get_out_edge_count();
					curr.erase(val);
					return curr.get_in_edge_count() + curr.get_out_edge_count() - before;

			});
		}

	public:
		//================== find edges given predicate.  predicate is also a filter operation that modifies the results..
		/// find edges, only ones that exist are returned.  this is probably not the most efficient
		/// however, the alternative would be the dig into distributed hash map to allow per element filter.
		std::vector<std::pair<key_type, typename Edge::EdgeInputType> > 
		find_edges(std::vector<std::pair<key_type, typename Edge::EdgeInputType>> &query) const {
			local_container_type localc;

			{
				// first get kmer vector.
				std::vector<key_type> kmers;   
				kmers.reserve(query.size());
				for (size_t i = 0; i < query.size(); ++i) {
					kmers.emplace_back(query[i].first);
				}

				// then do a standard query and insert results into a local hash table.  TODO: if we have one-to-one matching, then this is not needed.
				auto res = this->find(kmers);
				localc.insert(res);
			}

			// canonicalize so it matches with the local input.
			this->transform_input(query);

			// then match them up 
			typename local_container_type::iterator it;
			std::vector<std::pair<key_type, typename Edge::EdgeInputType> > results;
			uint8_t ch;
			for (size_t i = 0; i < query.size(); ++i) {
				it = localc.find(query[i].first);

				if (it != localc.end()) {
					// compare left and right.
					ch = query[i].second.getCharsAtPos(1, 1);
					if ((ch > 0) && ((ch & (ch-1)) == 0)) {  // has 1 and only 1 character at the in edge position (DNA16)
						if ((*it).second.get_in_edge_frequency(key_type::KmerAlphabet::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[ch]]) == 0)
							query[i].second.setCharsAtPos(static_cast<uint8_t>(0), 1, 1); 
					}
					ch = query[i].second.getCharsAtPos(0, 1);
					if ((ch > 0) && ((ch & (ch-1)) == 0)) {  // has 1 and only 1 character at the in edge position (DNA16)
						if ((*it).second.get_out_edge_frequency(key_type::KmerAlphabet::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[ch]]) == 0)
							query[i].second.setCharsAtPos(static_cast<uint8_t>(0), 0, 1); 
					}
					if (query[i].second.getData()[0] > 0) {
						results.emplace_back(query[i]);
					}
				}
			}
			return results;
		}

		/// matching edges are kept.  others are removed.
		std::vector<mutable_value_type> 
		find_edges(std::vector<mutable_value_type> &query) const {
			local_container_type localc;

			{
				// first get kmer vector.
				std::vector<key_type> kmers;   
				kmers.reserve(query.size());
				for (size_t i = 0; i < query.size(); ++i) {
					kmers.emplace_back(query[i].first);
				}

				// then do a standard query and insert results into a local hash table.  TODO: if we have one-to-one matching, then this is not needed.
				auto res = this->find(kmers);
				localc.insert(res);
			}

			// canonicalize so it matches with the local input.
			this->transform_input(query);

			// then match them up 
			typename local_container_type::iterator it;
			std::vector<std::pair<key_type, typename Edge::EdgeInputType> > results;
			for (size_t i = 0; i < query.size(); ++i) {
				it = localc.find(query[i].first);

				if (it != localc.end()) {
					query[i].second.intersect((*it).second);

					// compare left and right.
					if ((query[i].second.get_in_edge_count() + query[i].second.get_out_edge_count()) > 0) {
						results.emplace_back(query[i]);
					}
				}
			}
			return results;
		}

		/// find edges matching the specified edges and satisfy the predicates.
		template <typename Predicate>
		std::vector<std::pair<key_type, typename Edge::EdgeInputType> > 
		find_edges(std::vector<std::pair<key_type, typename Edge::EdgeInputType> > &query, 
		Predicate const &pred) const {
			local_container_type localc;

			{
				// first get kmer vector.
				std::vector<key_type> kmers;   
				kmers.reserve(query.size());
				for (size_t i = 0; i < query.size(); ++i) {
					kmers.emplace_back(query[i].first);
				}

				// then do a standard query and insert results into a local hash table.  TODO: if we have one-to-one matching, then this is not needed.
				auto res = this->find(kmers);
				localc.insert(res);
			}

			// canonicalize so it matches with the local input.
			this->transform_input(query);

			// then match them up 
			typename local_container_type::iterator it;
			std::vector<std::pair<key_type, typename Edge::EdgeInputType> > results;
			uint8_t ch;
			mapped_type val;
			for (size_t i = 0; i < query.size(); ++i) {
				it = localc.find(query[i].first);
				
				if (it != localc.end()) {
					val = (*it).second.select(pred);
					// compare left and right.
					ch = query[i].second.getCharsAtPos(1, 1);
					if ((ch > 0) && ((ch & (ch-1)) == 0)) {  // has 1 and only 1 character at the in edge position (DNA16)
						if (val.get_in_edge_frequency(key_type::KmerAlphabet::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[ch]]) == 0)
							query[i].second.setCharsAtPos(static_cast<uint8_t>(0), 1, 1); 
					}
					ch = query[i].second.getCharsAtPos(0, 1);
					if ((ch > 0) && ((ch & (ch-1)) == 0)) {  // has 1 and only 1 character at the in edge position (DNA16)
						if (val.get_out_edge_frequency(key_type::KmerAlphabet::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[ch]]) == 0)
							query[i].second.setCharsAtPos(static_cast<uint8_t>(0), 0, 1); 
					}
					if (query[i].second.getData()[0] > 0) {
						results.emplace_back(query[i]);
					}
				}
			}
			return results;
		}

		/// matching edges are kept.  others are removed.
		template <typename Predicate>
		std::vector<mutable_value_type> 
		find_edges(std::vector<mutable_value_type> &query, 
		Predicate const &pred) const {
			local_container_type localc;

			{
				// first get kmer vector.
				std::vector<key_type> kmers;   
				kmers.reserve(query.size());
				for (size_t i = 0; i < query.size(); ++i) {
					kmers.emplace_back(query[i].first);
				}

				// then do a standard query and insert results into a local hash table.  TODO: if we have one-to-one matching, then this is not needed.
				auto res = this->find(kmers);
				localc.insert(res);
			}

			// canonicalize so it matches with the local input.
			this->transform_input(query);

			// then match them up 
			typename local_container_type::iterator it;
			std::vector<std::pair<key_type, typename Edge::EdgeInputType> > results;
			for (size_t i = 0; i < query.size(); ++i) {
				it = localc.find(query[i].first);

				if (it != localc.end()) {
					query[i].second.intersect((*it).second.select(pred));

					// compare left and right.
					if ((query[i].second.get_in_edge_count() + query[i].second.get_out_edge_count()) > 0) {
						results.emplace_back(query[i]);
					}
				}
			}
			return results;
		}

		// find edges matching criteria.  LOCAL.
		template <typename Predicate>
		std::vector<mutable_value_type> find_edges(Predicate const &pred) const {

          ::std::vector<mutable_value_type > results;

          if (this->local_empty()) {
            //printf("rank %d local is empty\n", this->comm.rank());
            return results;
          }
          results.reserve(this->c.size() / 2);
		  mutable_value_type val;
          for (auto it = this->c.begin(); it != this->c.end(); ++it) {
			// make a copy.
			val.first = (*it).first;
			val.second = (*it).second.select(pred);  // check the value field.
			// if some edges are left, and push into results.
			if ((val.second.get_out_edge_count() + val.second.get_in_edge_count()) > 0) {
				results.emplace_back(val);
			} 
          }

          return results;
		}

		// ALL erase functions need to erase edges in both directions, so query first then erase in both directions.

	protected:
		//================== erase edges.  WE ASSUME THAT REVERSE EDGES ALSO NEED TO BE REMOVED.
		/// unconditional erase.  no predicate.  NOTE THAT DEST NODE's FREQUENCY IS NOT MAINTAINED.
		void erase_edges_internal(std::vector<std::pair<key_type, typename Edge::EdgeInputType> > & query) {
			// get the reverse query  - these should be simple_biedges.
			{
				std::vector<std::pair<key_type, typename Edge::EdgeInputType> > rev_edges =
					this->get_remote_edges(query);  // do before distributed query - to avoid query being modified.

				this->erase_simple_biedges(rev_edges);  // distributed.
			}

			this->erase_simple_biedges(query);  // distributed.

		}
		void erase_edges_internal(std::vector<mutable_value_type > & query) {
			// get the reverse query  - these should be simple_biedges.
			{
				std::vector<std::pair<key_type, typename Edge::EdgeInputType> > rev_edges =
					this->get_remote_edges(query);  // do before distributed query - to avoid query being modified.

				this->erase_simple_biedges(rev_edges);
			}

			this->erase_multi_biedges(query);
		}
	public:
		void erase_edges(std::vector<std::pair<key_type, typename Edge::EdgeInputType> > & query) {
			erase_edges_internal(query);
		}
		void erase_edges(std::vector<mutable_value_type > & query) {
			erase_edges_internal(query);
		}
		/// conditional erase. predicate operates on the edges.
		/// use for deleting specified edges that also is low frequency, for example.
		/// predicate is BINARY, , using val to select data in compact_multi_biedge for evaluation.
		template <typename Predicate>
		void erase_edges(std::vector<std::pair<key_type, typename Edge::EdgeInputType> > & query, Predicate const &pred) {
			// apply predicate to find edges and not erase_edges_internal, so that we don't have half edges left.
			auto res = this->find_edges(query, pred);
			this->erase_edges_internal(res);
		}

		template <typename Predicate>
		void erase_edges(std::vector<mutable_value_type > & query, Predicate const &pred) {
			// apply predicate to find edges and not erase_edges, so that we don't have half edges left.
			auto res = this->find_edges(query, pred);
			this->erase_edges_internal(res);
		}

		/// predicate evaluates for each edge.  only edges matching predicate are erased.
		/// used for deleting all nodes with low frequency, for example.
		/// predicate needs to provide a UNARY operator that select the nodes (param is compact_multi_biedge)
		///    and a UNARY operator that selects the 
		template <typename Predicate>
		void erase_edges(Predicate const &pred) {
			// apply predicate to find edges and not erase_edges, so that we don't have half edges left.
			std::vector<mutable_value_type > query = this->find_edges(pred);   // local.

			{
				std::vector<std::pair<key_type, typename Edge::EdgeInputType> > rev_edges =
					this->get_remote_edges(query);  // do before distributed query - to avoid query being modified.
				this->erase_simple_biedges(rev_edges);     // distributed.
			}
			// this part is local, because we know the query refer to local content.
			this->c.update(query,
				[](mapped_type & curr, mapped_type const & val){
					size_t before = curr.get_in_edge_count() + curr.get_out_edge_count();
					curr.erase(val);
					return curr.get_in_edge_count() + curr.get_out_edge_count() - before;
			});
		}
		//=================== erase nodes.  maintains graph consistency. =================================

	protected:

		/// unconditional erase.  no predicate.  NOTE THAT DEST NODE's FREQUENCY IS NOT MAINTAINED.
		void erase_nodes_and_edges(std::vector<mutable_value_type > & q) {
			// now erase the nodes.
			{
				// ===  first find the nodes and get the edges.
				std::vector<key_type> q2;
				q2.reserve(q.size());
				for (size_t i = 0; i < q.size(); ++i) {
					q2.push_back(q[i].first);
				}
				// ===  finally erase the nodes
				this->erase(q2);    // distributed
			}

			// erase the reverse edges.
			this->erase_simple_biedges(this->get_remote_edges(q));  // distributed
		}
		
	public:


		/// unconditional erase.  no predicate.  NOTE THAT DEST NODE's FREQUENCY IS NOT MAINTAINED.
		void erase_nodes(std::vector<key_type> &query) {
			this->erase_nodes_and_edges(this->find(query));    // distributed
		}

		/// conditional erase. predicate is UNARY and operates on the nodes, including all edges of the nodes.
		/// use for deleting specified nodes that also is low frequency, for example.
		template <typename Predicate>
		void erase_nodes(std::vector<key_type> &query, Predicate const &pred) {
			this->erase_nodes_and_edges(this->find(query, false, pred));  // distributed
		}

		/// predicate is UNARY evaluates for the entire node.  all edges of the nodes are erased.
		/// used for deleting all nodes with low frequency, for example.
		template <typename Predicate>
		void erase_nodes(Predicate const &pred) {
			// ===  erase remote edges of the matching nodes.
			this->erase_simple_biedges(  // distributed
				this->get_remote_edges(
					this->find(pred)));  // local

			// remove the nodes
			this->erase(pred);   // local.
		}




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
	//  		NOTE: that user specified frequency is for canonicalized (should be this way else could get mismatched forward and backward k and k2).
	//				k2, k1, and k counts can be computed either canonicalized or uncanonicalized.
	//      NOTE: k, k1, and k2 mers may be palindromic.  for k2 mer to be palindromic, k mer also needs to be.  
				palindromic k-mer is always canonical - and so is k2-mer.  no canonicalization would be applied and the edge frequencies can exist as is.
				palindromic k1-mer means that 2 successive k2-mer have central k-mers that are on opposite strands, and when canonicalized, the palindromic k1-mer would be double counted.
				freqency of k1-mer palindrome in a k2-mer should be added to both sides and then halved.
				if both sides are k1-palindromes, each k1-mer is doubled in frequency (from prev and next k2-mers) so each must be halved and added to both sides.
				unless the central k-mer is also a palindrome (DNA5 and DNA16), then no canonicalization happens, and we can skip the halving and adding to both sides, in order to preserve freq in real data
				another special case occurs when some k1-palindrome and some non-palindrome share the same central k-mer.  if they are counted and thresholded separately,
					the the total count would be off for common edges, as threshold is applied to partial counts.
					the solution to this special case is to count k2-mers together without saturation (uint32_t) and then do threshold and saturation only during graph freq aggregation.
				assume count is less than 4B so uint32_t works.  (if palindromes and k2-mers are counted separately, a sort and merge would need to be done
					only to save 3 bytes per unique k2-mer (subject to load factor of count map.))

	//  broken up into 2 functions so the k2mer counting can work with limited memory.
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
			// insert into counter.  need to make this a saturating counter.
			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> kv = 
				std::make_pair(std::make_pair(key_type(), typename Edge::EdgeInputType()), K2FreqType(1));
			for (auto it = input.begin(); it != input.end(); ++it) {
 
				// no saturation add.  K2Freq is uint32_t.  
				// NO check for rc_palindrome.  could be k, k+1, or k+2 (requires k)

				kv.first = *it;
				auto result = k2_counter.insert(kv);
				if (!(result.second)) {
					++(result.first->second);
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
			BL_BENCH_COLLECTIVE_END(k2freq, "init", block_size, this->comm);

			if (this->comm.rank() == 0)
				std::cout << "rank " << this->comm.rank() <<
				" BEFFORE size=" << this->local_size() << " buckets=" << this->c.bucket_count() << std::endl;

			BL_BENCH_LOOP_START(k2freq, 0);  // for parse and transform
			BL_BENCH_LOOP_START(k2freq, 1);  // for distribute
			BL_BENCH_LOOP_START(k2freq, 2);  // for insert into k2mer counter
			BL_BENCH_LOOP_START(k2freq, 3);  // for cleanup

			std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> kv = 
				std::make_pair(InputElemType(), K2FreqType(1));
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

					// no saturation add.  K2Freq is uint32_t.  
					// NO check for rc_palindrome.  could be k, k+1, or k+2 (requires k)
					kv.first = *it;
					auto result = k2_counter.insert(kv);
					if (!(result.second)) {
						++(result.first->second);
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

	// actual insertion using k2frequencies for palindromes.   actually, cant do separately from normal - filter on WHOLE GROUP of same kmers
	// seaparate would require storing partial sums.
	template <typename Predicate = ::bliss::filter::TruePredicate>
	size_t local_insert_by_freqencies(
		LocalK2CountMapType & k2_counter,
		std::vector<size_t> const & threshes,  // order is k_low, k_hi, k1_low, k1_hi, k2_low, k2_hi
		Predicate const & pred = Predicate()) {

		BL_BENCH_INIT(local_insert);

		// step 2: get the content of the k2 count map.
		BL_BENCH_START(local_insert);
		::std::vector<::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> > k2counts;
		k2_counter.to_vector(k2counts);
		k2_counter.reset(); // get rid of allocated space.
		BL_BENCH_END(local_insert, "to_vec", k2counts.size());

		// step 3: sort by kmer.  no transform of the kmer.
		BL_BENCH_START(local_insert);
		std::sort(k2counts.begin(), k2counts.end(),
			[](::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> const & x,
			   ::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> const & y) {
				   return x.first.first < y.first.first;
		});
		BL_BENCH_END(local_insert, "k2mer_sort", k2counts.size());

		// key_type trouble;
		// CTATATACTATATACTAGTATATAGTATATA
		// trouble.nextFromChar(1);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(1);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(1);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(2);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(2);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
		// trouble.nextFromChar(3);
		// trouble.nextFromChar(0);
/*
		// TACTATATACTAGTATATAGTATATTAAATA
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(1);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(1);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(2);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(2);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(0);
		 trouble.nextFromChar(3);
		 trouble.nextFromChar(0);
*/


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
				K2FreqType cnt, cnt_unsplit;
				uint8_t edges;
				auto block_it = k2counts.begin();
				auto it = k2counts.begin();
				size_t t0_lo = threshes[0],
					t0_hi = threshes[1],
					t1_lo = threshes[2],
					t1_hi = threshes[3],
					t2_lo = threshes[4],
					t2_hi = threshes[5];
				size_t cnt_max = ::std::numeric_limits<FreqType>::max();
		    	size_t thresh_max = cnt_max + 1;
				bool k_is_palindrome = false, km1_low_is_palindrome = false, km1_high_is_palindrome = false;
				if (it != k2counts.end()) {
					k = it->first.first;  // init
					k_is_palindrome = ::bliss::common::kmer::kmer_traits<key_type>::is_rc_palindrome(k);
					km1_high_is_palindrome = ::bliss::common::kmer::kmer_traits<key_type>::is_kminus1_rc_palindrome_high(k);
					km1_low_is_palindrome = ::bliss::common::kmer::kmer_traits<key_type>::is_kminus1_rc_palindrome_low(k);
				}
				uint8_t ch;
				bool split_counts = false;

				for (; it != k2counts.end(); ++it) {
					if (it->first.first != k) {  // start of new block

/*					 if ((k == trouble) || (k == trouble.reverse_complement())) {
					 	std::cout << "summary " <<
							k << " f " << k_f << " in [" <<
					 		static_cast<size_t>(k1in_f[0]) << "," << 
					 		static_cast<size_t>(k1in_f[1]) << "," <<
					 		static_cast<size_t>(k1in_f[2]) << "," <<
					 		static_cast<size_t>(k1in_f[3]) << "], out [" <<
					 		static_cast<size_t>(k1out_f[0]) << "," << 
					 		static_cast<size_t>(k1out_f[1]) << "," <<
					 		static_cast<size_t>(k1out_f[2]) << "," <<
					 		static_cast<size_t>(k1out_f[3]) << "]" << std::endl;
					 }

*/					  // cap the counts, could be data type max, and larger than thresh_max.
					  k_f = std::min(cnt_max, k_f);
					  // average by 2 since this is from palindromic k+1 mer.
					  	for (ch = 0; ch < 4; ++ch) {
							k1in_f[ch] = std::min(cnt_max, k1in_f[ch]);
							k1out_f[ch] = std::min(cnt_max, k1out_f[ch]);
						  }

						// accumulation for block complete. filter and record results now.
						for (; block_it != it; ++block_it, ++flag_it) {
//					 if ((k == trouble) || (k == trouble.reverse_complement())) {
//							 	std::cout << (*block_it) <<  std::endl;
//							 }

							// step 5: do per block threshold
							cnt = std::min(cnt_max, static_cast<size_t>(block_it->second));	

							// step 5a: filter k2mer 
							if ((cnt < t2_lo) || (cnt >= t2_hi)) *flag_it |= 4;
							if ((cnt >= t2_hi) && (this->comm.rank()==0)) std::cout << block_it->first << " k2 count " << cnt << std::endl;

							// step 5b: filter kmer
							if ((k_f < t0_lo) || (k_f >= t0_hi)) *flag_it |= 8;
							if ((k_f >= t0_hi) && (this->comm.rank()==0)) std::cout << block_it->first << " k count " << k_f << std::endl;
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
							if ((k1cnt >= t1_hi) && (k1cnt != thresh_max) && (this->comm.rank()==0)) std::cout << block_it->first << " in k1 count " << k1cnt << " t1_hi= " << t1_hi << std::endl;

							switch (edges & 0x0F) {
								case 0x01: k1cnt = k1out_f[0]; break;
								case 0x02: k1cnt = k1out_f[1]; break;
								case 0x04: k1cnt = k1out_f[2]; break;
								case 0x08: k1cnt = k1out_f[3]; break;
								default: k1cnt =  t1_hi; break;
							}
							if ((k1cnt < t1_lo) || (k1cnt >= t1_hi)) *flag_it |= 1;  // out
							if ((k1cnt >= t1_hi) && (k1cnt != thresh_max) && (this->comm.rank()==0)) std::cout << block_it->first << " out k1 count " << k1cnt << " t1_hi= " << t1_hi << std::endl;

						}
						
						// reinitialize for next block
						k = block_it->first.first;
						k_is_palindrome = ::bliss::common::kmer::kmer_traits<key_type>::is_rc_palindrome(k);
						km1_high_is_palindrome = ::bliss::common::kmer::kmer_traits<key_type>::is_kminus1_rc_palindrome_high(k);
						km1_low_is_palindrome = ::bliss::common::kmer::kmer_traits<key_type>::is_kminus1_rc_palindrome_low(k);
						k_f = 0;
						memset(k1in_f, 0, sizeof(size_t) << 2);
						memset(k1out_f, 0, sizeof(size_t) << 2);
					}

					// step 4: do per block sum
					cnt_unsplit = it->second;  // show frequency in real data, so leave as is, even k1palindromic ones.
					// step 4a: accumulate the kmer counts. 
					k_f += cnt_unsplit;   // accumulate the kmer count  - each k+2 mer occurs once, even palindromic ones.
					// step 4b: accumulate the k1mer in and out edge counts.
					// if k is palindromic, then revcomp of edges may be present.  AVERAGE. 
					edges = it->first.second.getData()[0];

					ch = edges >> 4;
					if (ch > 0) {
					  if ((ch & (ch - 1)) == 0) {  // power of 2 
					    ch = key_type::KmerAlphabet::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[ch]]; // convert to target alphabet
					    split_counts = ::bliss::common::kmer::kmer_traits<key_type>::is_k1_rc_palindrome(ch, k, km1_high_is_palindrome) && (!k_is_palindrome);
//					 if ((k == trouble) || (k == trouble.reverse_complement())) {
//					 	std::cout << "km1 high palindrome ? " << (km1_high_is_palindrome ? "Y" : "N") << 
//					 	" k palindrome ? " << (k_is_palindrome ? "Y" : "N")  << " ch " << 
//					 	static_cast<size_t>(ch) << ", low val " << k.getCharsAtPos(0, 1) << " split ? " << (split_counts ? "Y" : "N") << std::endl;
//					    }
					    cnt = split_counts ? ((cnt_unsplit+1) >> 1) : cnt_unsplit; 
					    k1in_f[ch] += cnt;
		//			    if (split_counts) {  // TODO: only handling k+1 palindrome right now.
		//				k1out_f[(3-ch)] += (cnt_unsplit - cnt);
		//			    }
					  } else {
					    throw std::logic_error("combination of edges, so not a k2mer.");
					  }
					}
					

					ch = edges & 0xF;
					if (ch > 0) {
					  if ((ch & (ch - 1)) == 0) { // power of 2
					    ch = key_type::KmerAlphabet::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[ch]];
					    split_counts = ::bliss::common::kmer::kmer_traits<key_type>::is_k1_rc_palindrome(k, ch, km1_low_is_palindrome) && (!k_is_palindrome);
//					 if ((k == trouble) || (k == trouble.reverse_complement())) {
//					 	std::cout << "km1 low palindrome ? " << (km1_low_is_palindrome ? "Y" : "N") << 
//					 	" k palindrome ? " << (k_is_palindrome ? "Y" : "N")  << " ch " << 
//					 	static_cast<size_t>(ch) << ", high val " << k.getCharsAtPos(key_type::size - 1, 1) << " split ? " << (split_counts ? "Y" : "N") << std::endl;
//					    }
					    cnt = split_counts ? ((cnt_unsplit+1) >> 1) : cnt_unsplit; 
					    k1out_f[ch] += cnt;
		//			    if (split_counts) {  // TODO: only handling k+1 palindrome right now.
		//				k1in_f[3-ch] += (cnt_unsplit - cnt);
		//			    }
					  } else {
					    throw std::logic_error("combination of edges, so not a k2mer." );
					  }
					}
				}	
				// do the final block
				// cap the counts
				k_f = std::min(cnt_max, k_f);
				// since palindromic, average the two.
				for (ch = 0; ch < 4; ++ch) {
					k1in_f[ch] = std::min(cnt_max, k1in_f[ch]);
					k1out_f[ch] = std::min(cnt_max, k1out_f[ch]);
				}

				for (; block_it != it; ++block_it, ++flag_it) {
//					 if ((k == trouble) || (k == trouble.reverse_complement())) {
//					 	std::cout << (*block_it) <<  std::endl;
//					 }

					// step 5: do per block threshold
					cnt = std::min(cnt_max, static_cast<size_t>(block_it->second));
					
					// step 5a: filter k2mer 
					if ((cnt < t2_lo) || (cnt >= t2_hi)) *flag_it |= 4;
							if (cnt >= t2_hi) std::cout << block_it->first << " k2 count " << cnt << std::endl;
					// step 5b: filter kmer
					if ((k_f < t0_lo) || (k_f >= t0_hi)) *flag_it |= 8;
							if (k_f >= t0_hi) std::cout << block_it->first << " k count " << k_f << std::endl;
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
          if ((k1cnt >= t1_hi) && (k1cnt != thresh_max) && (this->comm.rank()==0)) std::cout << block_it->first << " in k1 count " << k1cnt << " t1_hi= " << t1_hi << std::endl;

					switch (edges & 0x0F) {
						case 0x01: k1cnt = k1out_f[0]; break;
						case 0x02: k1cnt = k1out_f[1]; break;
						case 0x04: k1cnt = k1out_f[2]; break;
						case 0x08: k1cnt = k1out_f[3]; break;
						default: k1cnt =  t1_hi; break;
					}
					if ((k1cnt < t1_lo) || (k1cnt >= t1_hi)) *flag_it |= 1;  // out
          if ((k1cnt >= t1_hi) && (k1cnt != thresh_max) && (this->comm.rank()==0)) std::cout << block_it->first << " out k1 count " << k1cnt << " t1_hi= " << t1_hi << std::endl;
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
						memmove(&(*write_it), it, sizeof(::std::pair<::std::pair<key_type, typename Edge::EdgeInputType>, K2FreqType> ));
						switch (flag) {
							case 1:
								write_it->first.second.getDataRef()[0] &= 0xF0;
								break;  // remove out and keep in
							case 2: 
								write_it->first.second.getDataRef()[0] &= 0x0F;
								break;  // remove in and keep out
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

template<typename Kmer, template <typename> MapParams = ::bliss::debruijn::CanonicalDeBruijnHashMapParamsMurmur>
using simple_hash_debruijn_graph_map = ::bliss::debruijn::graph::debruijn_graph_map<Kmer,
		::bliss::debruijn::graph::compact_multi_biedge<typename Kmer::KmerAlphabet, bool>,
		MapParams>;

template<typename Kmer, typename Count = uint32_t, template <typename> MapParams = ::bliss::debruijn::CanonicalDeBruijnHashMapParamsMurmur >
using count_hash_debruijn_graph_map = ::bliss::debruijn::graph::debruijn_graph_map<Kmer,
		::bliss::debruijn::graph::compact_multi_biedge<typename Kmer::KmerAlphabet, Count>,
		MapParams>;

} // namespace graph
}/*debruijn*/
}/*bliss*/

#endif /* DEBRUIJN_GRAPH_MAP_HPP_ */
