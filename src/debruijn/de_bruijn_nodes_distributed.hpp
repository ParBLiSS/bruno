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
 * de_bruijn_nodes_distributed.hpp
 *
 *  Created on: Aug 6, 2015
 *      Author: yongchao
 */

#ifndef DE_BRUIJN_NODES_DISTRIBUTED_HPP_
#define DE_BRUIJN_NODES_DISTRIBUTED_HPP_

#include "bliss-config.hpp"

#include <unordered_map>  // local storage hash table  // for multimap
#include <unordered_set>  // local storage hash table  // for multimap
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // not a multimap, where we need it most.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.

#include <type_traits>
#include "debruijn/de_bruijn_node_trait.hpp"	//node trait data structure storing the linkage information to the node
#include "containers/distributed_map_base.hpp"
#include "containers/distributed_unordered_map.hpp"

#include "utils/benchmark_utils.hpp"  // for timing.
#include "utils/logging.h"
#include "utils/filter_utils.hpp"
#include "utils/transform_utils.hpp"


namespace bliss{
	namespace de_bruijn{


	// NOTE:  DOES NOT work with canonicalized. in other words, input transform should be Identity.

		template<typename Key, typename T,
			template <typename> class MapParams,
		class Alloc = ::std::allocator< ::std::pair<const Key, T> >
		  >
		  class de_bruijn_nodes_distributed : public ::dsc::unordered_map<Key, T, MapParams, Alloc> {
			  using Base = ::dsc::unordered_map<Key, T, MapParams, Alloc>;

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


			  /**
			   * @brief insert new elements in the distributed unordered_multimap.
			   * @param first
			   * @param last
			   */
			  template <class InputIterator>
			  size_t local_insert(InputIterator first, InputIterator last) {
				  
          if (first == last) return 0;
          
          int32_t relative_strand;
				  size_t before = this->c.size();

          

				  /*reserve space*/
				  this->local_reserve(before + ::std::distance(first, last));

				  // allocate a node
          auto node = this->c.find(first->first);

          using EdgeType = decltype(node->second);
          using Alphabet = typename EdgeType::Alphabet;

				  /*iterate each input tuple*/
				  for (auto it = first; it != last; ++it) {

				    node = this->c.find(it->first);

				    /*tranform from <key, int> to <node, node_info>*/
					  if(node == this->c.end()){
						  /*create a new node*/
						  auto ret = this->c.emplace(::std::make_pair(it->first, T()));
						  if(ret.second == false){
							  cerr << "Insertion failed at line " << __LINE__ << " in file " << __FILE__ << endl;
							  exit(-1);
						  }

#if 0
              ::std::cerr << "Not exist in the hash table" << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(ret.first->first)  << ::std::endl;
              ::std::cerr << "0x" << std::hex << it->second << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(it->first) << ::std::endl << ::std::endl;
#endif
						  node = ret.first;  // now use it.

	            // determine if the node in the graph has the same orientation as the one being inserted.
						  relative_strand = bliss::de_bruijn::node::SENSE;
					  } else {
#if 0
             /*update the node*/
              ::std::cerr << "Exist in the hash table" << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(node->first)  << ::std::endl;
              ::std::cerr << "0x" << std::hex << it->second << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(it->first) << ::std::endl << ::std::endl;
#endif

              // determine if the node in the graph has the same orientation as the one being inserted.
              relative_strand = (node->first == it->first) ? bliss::de_bruijn::node::SENSE : bliss::de_bruijn::node::ANTI_SENSE;

					  }


					  // if different, swap and reverse complement the edges.  else, use as is.
					  if (relative_strand == bliss::de_bruijn::node::ANTI_SENSE) {
              node->second.update(bliss::de_bruijn::node::input_edge_utils::reverse_complement_edges<Alphabet>(it->second));

					  } else {
              node->second.update(it->second);
					  }

				  }
				  return this->c.size() - before;
			  }

			  /**
			   * @brief insert new elements in the distributed unordered_multimap.
			   * @param first
			   * @param last
			   */
			  template <class InputIterator, class Predicate>
			  size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {


          if (first == last) return 0;
				  int32_t relative_strand;
				  size_t before = this->c.size();

				  this->local_reserve(before + ::std::distance(first, last));
				  auto node = this->c.find(first->first);

          using EdgeType = decltype(node->second);
          using Alphabet = typename EdgeType::Alphabet;

				  for (auto it = first; it != last; ++it) {
					if (pred(*it)) {

					  node = this->c.find(it->first);

					  /*tranform from <key, int> to <node, node_info*/
					  if(node == this->c.end()){
						  /*create a new node*/
						  auto ret = this->c.emplace(::std::make_pair(it->first, T()));
						  if(ret.second == false){
							  cerr << "Insertion failed at line " << __LINE__ << " in file " << __FILE__ << endl;
							  exit(-1);
						  }

						  node = ret.first;  // now use it.

						  /*update the node*/
						  relative_strand = bliss::de_bruijn::node::SENSE;
					  }else{
						 /*update the node*/
						  relative_strand = (node->first == it->first) ? bliss::de_bruijn::node::SENSE : bliss::de_bruijn::node::ANTI_SENSE;
					  }

            // if different, swap and reverse complement the edges.  else, use as is.
            if (relative_strand == bliss::de_bruijn::node::ANTI_SENSE) {
              node->second.update(bliss::de_bruijn::node::input_edge_utils::reverse_complement_edges<Alphabet>(it->second));

            } else {
              node->second.update(it->second);
            }


					}
				  }
				  return this->c.size() - before;

			  }

			public:
			  de_bruijn_nodes_distributed(const mxx::comm& _comm) : Base(_comm) {/*do nothing*/}

			  virtual ~de_bruijn_nodes_distributed() {/*do nothing*/};

			  /*transform function*/

			  /**
				* @brief insert new elements in the distributed unordered_multimap.
				* @param first
				* @param last
				*/
			   template <typename InputEdgeType, typename Predicate = ::bliss::filter::TruePredicate>
			   size_t insert(std::vector<::std::pair<Key, InputEdgeType> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
				 // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
				 BL_BENCH_INIT(insert);

//				 BL_BENCH_START(insert);
//				 this->transform_input(input);
//				 BL_BENCH_END(insert, "start", input.size());
				 static_assert(std::is_same<typename Base::Base::Base::InputTransform, ::bliss::transform::identity<Key>>::value,
						 "de bruijn graph does not support transform of input Kmers. (e.g. canonicalizing).  Hash can use transformed values, though.");


				 // communication part
				 if (this->comm.size() > 1) {
					 BL_BENCH_COLLECTIVE_START(insert, "distribute", this->comm);
				   ::std::vector<size_t> recv_counts =
						   ::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm);
				   BLISS_UNUSED(recv_counts);
				   BL_BENCH_END(insert, "distribute", input.size());

				 }

				 BL_BENCH_START(insert);
				 // local compute part.  called by the communicator.
				 size_t count = 0;
				 if (!::std::is_same<Predicate, ::bliss::filter::TruePredicate>::value)
				   count = this->local_insert(input.begin(), input.end(), pred);
				 else
				   count = this->local_insert(input.begin(), input.end());
				 BL_BENCH_END(insert, "insert", this->c.size());

				 BL_BENCH_REPORT_MPI(insert, this->comm.rank(), this->comm);

				 return count;
			   }
		};
	}/*de_bruijn*/
}/*bliss*/

#endif /* DE_BRUIJN_NODES_DISTRIBUTED_HPP_ */
