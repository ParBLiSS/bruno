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

/**
 * @file    distributed_hashed_vec.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the distributed_multimap, distributed map, and distributed_reduction_map
 *          data structures.
 *
 *          implementation is hash-base (O(1) lookup). later will support sort-based (load balanced).
 *
 *          for now, input and output via local vectors.
 *          (later, allow remote vectors,
 *            which can have remote ranges  (all to all to "sort" remote ranges to src proc,
 *            then src proc compute target procs for each elements in the remote ranges,
 *            communicate remote ranges to target procs.  target proc can then materialize the data.
 *            may not be efficient though if we don't have local spatial coherence..
 *          )
 *
 *          most create-find-delete operations support remote filtering via predicates.
 *          most create-find-delete oeprations support remote transformation.
 *
 *          signature of predicate is bool pred(T&).  if predicate needs to access the local map, it should be done via its constructor.
 *

 */

#ifndef BLISS_DISTRIBUTED_HASHED_VEC_HPP
#define BLISS_DISTRIBUTED_HASHED_VEC_HPP


#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // not a multimap, where we need it most.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.

#include <type_traits>

#include <mxx/collective.hpp>
#include <mxx/reduction.hpp>
#include <mxx/algos.hpp> // for bucketing

#include "containers/distributed_map_base.hpp"
#include "containers/unordered_vecmap.hpp"
#include "containers/hashed_vecmap.hpp"

#include "utils/benchmark_utils.hpp"  // for timing.
#include "utils/logging.h"

#include "common/kmer_transform.hpp"

#include "containers/dsc_container_utils.hpp"

namespace dsc  // distributed std container
{

	// =================
	// NOTE: when using this, need to further alias so that only Key param remains.
	// =================


  /**
   * @brief  distributed unordered map following std unordered map's interface.
   * @details   This class is modeled after the std::unordered_map.
   *         it has as much of the same methods of std::unordered_map as possible.  however, all methods consider the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.  Also since we
   *         are working with 'distributed' data, batched operations are preferred.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   *  this class and its subclasses rely on 2 hash function for data distribution and 1 equal and 1 less comparators.  these are specified
   *  as template parameters.  the less comparator is for preprocessing queries and inserts.  keys (e.g.) kmers are transformed before they
   *  are hashed/compared.
   *    an alternative approach is to hold only canonical keys in the map.
   *    yet another alternative approach is to perform 2 queries for every key.- 2x computation but communication is spread out.
   *
   *  note: KeyTransform is applied before Hash, and Equal operators.  These operators should have NO KNOWLEDGE of any transform applied, including kmolecule to kmer mapping.
   *
   *  any operation that uses sort is not going to scale well.  this includes "hash_unique_key, hash_unique_tuple, local_reduction"...
   *  to reduce this, we can try by using a hash set instead.  http://www.vldb.org/pvldb/2/vldb09-257.pdf, http://www.vldb.org/pvldb/vol7/p85-balkesen.pdf
   *
   *  conditional version of insert/erase/find/count supports predicate that operate on INTERMEDIATE RESULTS.  input (key,key-value pair) can be pre-filtered.
   *    output (query result, e.g.) can be post filtered (and optionally reduce comm volume).
   *    intermediate results (such as counting in multimap only if T has certain value) can only be filtered at during local_operation.
   *
   *
   * key to proc assignment can be done as hash or splitters in sorted range.
   * tuples can be sotred in hash table or as sorted array.
   *   hash-hash combination works
   *  sort-sort combination works as well
   *  hash-sort combination can work.  advantage is in range query.
   *  sort-hash combination would be expensive for updating splitters
   *
   * This version is the hash-hash.
   *
   * @tparam Key
   * @tparam T
   * @tparam Container  default to unordered_map and unordered multimap, requiring 5 template params.
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
    template <typename, typename, typename, typename, typename...> class Container,
    template <typename> class MapParams,
    class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class hashed_vec_base :
		  public ::dsc::map_base<Key, T, MapParams, Alloc> {

    protected:
      using Base = ::dsc::map_base<Key, T, MapParams, Alloc>;

//      using TransformedHash = ::fsc::TransformedHash<Key, Hash<Key, false>, KeyTransform>;
//      TransformedHash hash;

      struct KeyToRank {
          typename Base::DistTransformedFunc proc_trans_hash;
          const int p;

          // 2x comm size to allow more even distribution?
          KeyToRank(int comm_size) :
        	  proc_trans_hash(typename Base::DistFunc(ceilLog2(comm_size)),
        			  	  	  typename Base::DistTrans()),
        			  p(comm_size) {};

          inline int operator()(Key const & x) const {
            //            printf("KeyToRank operator. commsize %d  key.  hashed to %d, mapped to proc %d \n", p, proc_hash(Base::trans(x)), proc_hash(Base::trans(x)) % p);
            return proc_trans_hash(x) % p;
          }
          template<typename V>
          inline int operator()(::std::pair<Key, V> const & x) const {
            return this->operator()(x.first);
          }
          template<typename V>
          inline int operator()(::std::pair<const Key, V> const & x) const {
            return this->operator()(x.first);
          }
      } key_to_rank;


      /**
       * @brief count elements with the specified keys in the distributed sorted_multimap.
       * @note  input cannot have duplicate elements.
       *
       * @param first
       * @param last
       */
      struct QueryProcessor {  // assume unique, always.

          // assumes that container is sorted. and exact overlap region is provided.  do not filter output here since it's an output iterator.
          template <class DB, class QueryIter, class OutputIter, class Operator, class Predicate = ::fsc::TruePredicate>
          static size_t process(DB &db,
                                QueryIter query_begin, QueryIter query_end,
                                OutputIter &output, Operator & op,
                                bool sorted_query = false, Predicate const &pred = Predicate()) {

              if (query_begin == query_end) return 0;

              size_t count = 0;  // before size.
              if (!::std::is_same<Predicate, ::fsc::TruePredicate>::value)
                for (auto it = query_begin; it != query_end; ++it) {
                  count += op(db, *it, output, pred);
                }
              else
                for (auto it = query_begin; it != query_end; ++it) {
                  count += op(db, *it, output);
                }
              return count;
          }

      };




    public:
      using local_container_type = Container<Key, T,
    		  typename Base::StoreTransformedFunc,
    		  typename Base::StoreTransformedLess,
    		  typename Base::StoreTransformedEqual,
    		  Alloc>;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      local_container_type c;

      mutable bool local_changed;

      struct LocalCount {
          // unfiltered.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) const {
              *output = ::std::move(::std::make_pair(v, db.count(v)));
              ++output;
              return 1;
          }
          // filtered element-wise.
          template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &output,
                            Predicate const& pred) const {
              auto range = db.equal_range(v);

              // add the output entry.
              size_t count = 0;
              if (pred(range.first, range.second))  // operator to decide if range matches.
                count = ::std::count_if(range.first, range.second, pred);  // operator for each element in range.

              *output = ::std::move(::std::make_pair(v, count));
              ++output;
              return 1;
          }
          // no filter by range AND elemenet for now.
      } count_element;

      struct LocalErase {
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &) {
              size_t before = db.size();
              db.erase(v);
              return before - db.size();
          }
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &,
                            Predicate const & pred) {
              auto range = (const_cast<DB const &>(db)).equal_range(v);

              // check range first.  then erase.  only removed iterators are invalidated.
              // order of remaining elements preserved.  true for map/multimap, unordered or not.
              size_t count = 0;
              auto tmp = range.first;
              if (pred(range.first, range.second)) { // operator to decide if range matches.
                for (auto it = range.first; it != range.second;) {
                  if (pred(*it)) {
                    // advance, then remove.
                    tmp = it;  ++it;
                    // remove.
                    db.erase(tmp);  // erase entry at 1 iterator position. (not by value).
                    ++count;
                  } else {
                    // keep.  so just advance
                    ++it;
                  }
                }
              }

              return count;
          }
          // no filter by range AND elemenet for now.
      } erase_element;

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_insert(InputIterator first, InputIterator last) {
    	  BL_BENCH_INIT(local_insert);

    	  BL_BENCH_START(local_insert);
          this->local_reserve(c.size() + ::std::distance(first, last));  // before branching, because reserve calls collective "empty()"
          BL_BENCH_END(local_insert, "reserve", this->c.size());

          if (first == last) return 0;

          size_t before = c.size();

          BL_BENCH_START(local_insert);
          for (auto it = first; it != last; ++it) {
            c.emplace(*it);
          }
          BL_BENCH_END(local_insert, "emplace", this->c.size());

          if (c.size() != before) local_changed = true;


          BL_BENCH_REPORT_MPI_NAMED(local_insert, "base_hashmap:local_insert", this->comm);

          //          c.insert(first, last);  // mem usage?
          return c.size() - before;
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.  example use: stop inserting if more than x entries.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      size_t local_insert(InputIterator first, InputIterator last, Predicate const &pred) {

          auto new_end = std::partition(first, last, pred);
//
//          this->local_reserve(c.size() + ::std::distance(first, last));   // before branching, because reserve calls collective "empty()"
//
          if (first == last) return 0;

          size_t before = c.size();

          this->local_insert(first, new_end);
//
//          for (auto it = first; it != last; ++it) {
//            if (pred(*it)) c.emplace(*it);
//          }
//
          if (c.size() != before) local_changed = true;

          return c.size() - before;
      }





      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param keys  content will be changed and reordered
       * @param last
       */
      template <class LocalFind, typename Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_a2a(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
          BL_BENCH_INIT(find);

          ::std::vector<::std::pair<Key, T> > results;

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_a2a", this->comm);
            return results;
          }


          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_a2a", this->comm);
            return results;
          }


          BL_BENCH_START(find);
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          this->transform_input(keys);
          BL_BENCH_END(find, "input_transform", keys.size());

          if (this->comm.size() > 1) {

            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
            // distribute (communication part)
            std::vector<size_t> recv_counts(
            		::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
            				typename Base::StoreTransformedFunc(),
            				typename Base::StoreTransformedEqual()));
            BL_BENCH_END(find, "dist_query", keys.size());


            // local find. memory utilization a potential problem.
            // do for each src proc one at a time.

            BL_BENCH_START(find);
            // 1.5s for 32 nodes 6gb file type of slow....
            //            float multi = this->get_multiplicity();
            //            if (this->comm.rank() == 0) printf("rank %d multiplicity %f\n", this->comm.rank(), multi);
            results.reserve(keys.size());                   // TODO:  should estimate coverage.
            BL_BENCH_END(find, "reserve", results.capacity());

            BL_BENCH_START(find);
            std::vector<size_t> send_counts(this->comm.size(), 0);
            auto start = keys.begin();
            auto end = start;
            for (int i = 0; i < this->comm.size(); ++i) {
              ::std::advance(end, recv_counts[i]);

              // work on query from process i.
              send_counts[i] = QueryProcessor::process(c, start, end, emplace_iter, find_element, sorted_input, pred);
              // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);

              start = end;
            }
            BL_BENCH_END(find, "local_find", results.size());
            if (this->comm.rank() == 0) printf("rank %d result size %lu capacity %lu\n", this->comm.rank(), results.size(), results.capacity());


            BL_BENCH_COLLECTIVE_START(find, "a2a2", this->comm);
            // send back using the constructed recv count
            mxx::all2allv(results, send_counts, this->comm).swap(results);
            BL_BENCH_END(find, "a2a2", results.size());

          } else {

            BL_BENCH_START(find);
            // keep unique keys
            ::fsc::unique(keys, sorted_input,
            		typename Base::StoreTransformedFunc(),
            		typename Base::StoreTransformedEqual());
            BL_BENCH_END(find, "uniq1", keys.size());

            BL_BENCH_START(find);
//            float multi = this->get_multiplicity();
//            if (this->comm.rank() == 0) printf("rank %d multiplicity %f\n", this->comm.rank(), multi);
            results.reserve(keys.size());                   // TODO:  should estimate coverage.
            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
            BL_BENCH_END(find, "reserve", keys.capacity() );

            BL_BENCH_START(find);
            QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, sorted_input, pred);
            BL_BENCH_END(find, "local_find", results.size());

            if (this->comm.rank() == 0) printf("rank %d result size %lu capacity %lu\n", this->comm.rank(), results.size(), results.capacity());

          }

          BL_BENCH_REPORT_MPI_NAMED(find, "base_hashmap:find_a2a", this->comm);

          return results;

      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       *
       * why this version that uses isend and irecv?  because all2all version requires all result data to be in memory.
       * this one can do it one source process at a time.
       *
       * @param keys    content will be changed and reordered.
       * @param last
       */
      template <class LocalFind, typename Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_overlap(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
          BL_BENCH_INIT(find);

          ::std::vector<::std::pair<Key, T> > results;

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_overlap", this->comm);
            return results;
          }


          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_overlap", this->comm);
            return results;
          }


          BL_BENCH_START(find);
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;

          ::std::vector<::std::pair<Key, T> > local_results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > local_emplace_iter(local_results);
          this->transform_input(keys);
          BL_BENCH_END(find, "transform_input", keys.size());

          if (this->comm.size() > 1) {

            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
            // distribute (communication part)
            std::vector<size_t> recv_counts(
            		::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
            				typename Base::StoreTransformedFunc(),
            				typename Base::StoreTransformedEqual()));

            BL_BENCH_END(find, "dist_query", keys.size());


            //======= local count to determine amount of memory to allocate at destination.
            BL_BENCH_START(find);
            ::std::vector<::std::pair<Key, size_t> > count_results;
            size_t max_key_count = *(::std::max_element(recv_counts.begin(), recv_counts.end()));
            count_results.reserve(max_key_count);
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

            std::vector<size_t> send_counts(this->comm.size(), 0);

            auto start = keys.begin();
            auto end = start;
            size_t total = 0;
            for (int i = 0; i < this->comm.size(); ++i) {
              ::std::advance(end, recv_counts[i]);

              // count results for process i
              count_results.clear();
              QueryProcessor::process(c, start, end, count_emplace_iter, count_element, sorted_input, pred);
              send_counts[i] =
                  ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                    [](size_t v, ::std::pair<Key, size_t> const & x) {
                return v + x.second;
              });
              //            for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
              //              send_counts[i] += it->second;
              //            }
              total += send_counts[i];
              start = end;
              //printf("Rank %d local count for src rank %d:  recv %d send %d\n", this->comm.rank(), i, recv_counts[i], send_counts[i]);
            }
            ::std::vector<::std::pair<Key, size_t> >().swap(count_results);
            BL_BENCH_END(find, "local_count", total);


            BL_BENCH_COLLECTIVE_START(find, "a2a_count", this->comm);
            std::vector<size_t> resp_counts = mxx::all2all(send_counts, this->comm);  // compute counts of response to receive
            BL_BENCH_END(find, "a2a_count", keys.size());


            //==== reserve
            BL_BENCH_START(find);
            auto resp_displs = mxx::impl::get_displacements(resp_counts);  // compute response displacements.

            auto resp_total = resp_displs[this->comm.size() - 1] + resp_counts[this->comm.size() - 1];
            auto max_send_count = *(::std::max_element(send_counts.begin(), send_counts.end()));
            results.resize(resp_total);   // allocate, not just reserve
            local_results.reserve(max_send_count);

            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
            BL_BENCH_END(find, "reserve", resp_total);

            //=== process queries and send results.  O(p) iterations
            BL_BENCH_START(find);
            auto recv_displs = mxx::impl::get_displacements(recv_counts);  // compute response displacements.
            int recv_from, send_to;
            size_t found;
            total = 0;
            std::vector<MPI_Request> reqs(2 * this->comm.size());

            mxx::datatype dt = mxx::get_datatype<::std::pair<Key, T> >();

            for (int i = 0; i < this->comm.size(); ++i) {
              recv_from = (this->comm.rank() + (this->comm.size() - i)) % this->comm.size(); // rank to recv data from

              // set up receive.
              MPI_Irecv(&results[resp_displs[recv_from]], resp_counts[recv_from], dt.type(),
                        recv_from, i, this->comm, &reqs[2 * i]);


              send_to = (this->comm.rank() + i) % this->comm.size();    // rank to send data to

              //== get data for the dest rank
              start = keys.begin();                                   // keys for the query for the dest rank
              ::std::advance(start, recv_displs[send_to]);
              end = start;
              ::std::advance(end, recv_counts[send_to]);

              local_results.clear();
              // work on query from process i.
              found = QueryProcessor::process(c, start, end, local_emplace_iter, find_element, sorted_input, pred);
              // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);
              total += found;
              //== now send the results immediately - minimizing data usage so we need to wait for both send and recv to complete right now.

              MPI_Isend(&(local_results[0]), found, dt.type(), send_to,
                        i, this->comm, &reqs[2 * i + 1]);

              // wait for both requests to complete.
              MPI_Waitall(2, &reqs[2 * i], MPI_STATUSES_IGNORE);


              // verify correct? done by comparing to previous code.


              //printf("Rank %d local find send to %d:  query %d result sent %d (%d).  recv from %d received %d\n", this->comm.rank(), send_to, recv_counts[send_to], found, send_counts[send_to], recv_from, resp_counts[recv_from]);
            }
            //printf("Rank %d total find %lu\n", this->comm.rank(), total);
            BL_BENCH_END(find, "find_send", results.size());

          } else {

            BL_BENCH_START(find);
            // keep unique keys
            ::fsc::unique(keys, sorted_input,
    				typename Base::StoreTransformedFunc(),
    				typename Base::StoreTransformedEqual());
            BL_BENCH_END(find, "uniq1", keys.size());

            // memory is constrained.  find EXACT count.
            BL_BENCH_START(find);
            ::std::vector<::std::pair<Key, size_t> > count_results;
            count_results.reserve(keys.size());
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

            // count now.
            QueryProcessor::process(c, keys.begin(), keys.end(), count_emplace_iter, count_element, sorted_input, pred);
            size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                             [](size_t v, ::std::pair<Key, size_t> const & x) {
              return v + x.second;
            });
            //          for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
            //            count += it->second;
            //          }
            BL_BENCH_END(find, "local_count", count);

            BL_BENCH_START(find);
            results.reserve(count);                   // TODO:  should estimate coverage.
            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
            BL_BENCH_END(find, "reserve", results.capacity());

            BL_BENCH_START(find);
            QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, sorted_input, pred);
            BL_BENCH_END(find, "local_find", results.size());
          }

          BL_BENCH_REPORT_MPI_NAMED(find, "base_hashmap:find_overlap", this->comm);

          return results;

      }


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param keys  content will be changed and reordered
       * @param last
       */
      template <class LocalFind, typename Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
          BL_BENCH_INIT(find);

          ::std::vector<::std::pair<Key, T> > results;

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find", this->comm);
            return results;
          }

          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find", this->comm);
            return results;
          }

          BL_BENCH_START(find);
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          this->transform_input(keys);
          BL_BENCH_END(find, "input_transform", keys.size());

          if (this->comm.size() > 1) {

            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
            // distribute (communication part)
            std::vector<size_t> recv_counts(
            		::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
            				typename Base::StoreTransformedFunc(),
            				typename Base::StoreTransformedEqual()));
            BL_BENCH_END(find, "dist_query", keys.size());


            // local find. memory utilization a potential problem.
            // do for each src proc one at a time.

            BL_BENCH_START(find);
            results.reserve(keys.size() * 10);                   // TODO:  should estimate coverage.
            BL_BENCH_END(find, "reserve", results.capacity());

            BL_BENCH_START(find);
            std::vector<size_t> send_counts(this->comm.size(), 0);
            auto start = keys.begin();
            auto end = start;
            for (int i = 0; i < this->comm.size(); ++i) {
              ::std::advance(end, recv_counts[i]);

              // work on query from process i.
              send_counts[i] = QueryProcessor::process(c, start, end, emplace_iter, find_element, sorted_input, pred);
              // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);

              // use the first one to estimate the rest.
              //if (i == std::ceil(static_cast<double>(this->comm.size()) * 0.05)) {
			  if ((results.size() + send_counts[i]) > results.capacity()) {  // guess that next batch is going to get us similar size, and that may cause results to resize too much.
            	  // count so far
            	  size_t new_est = std::ceil((static_cast<double>(results.size()) / static_cast<double>(std::distance(keys.begin(), end))) * static_cast<double>(keys.size()) * 1.1f);
            	  if (this->comm.rank() == 0) printf("rank %d nkeys %lu nresuts %lu est result size %lu original estimate %lu\n", this->comm.rank(), keys.size(), results.size(), new_est, results.capacity());
            	  results.reserve(new_est);
              }

              start = end;
            }
            BL_BENCH_END(find, "local_find", results.size());
            if (this->comm.rank() == 0) printf("rank %d result size %lu capacity %lu\n", this->comm.rank(), results.size(), results.capacity());


            BL_BENCH_COLLECTIVE_START(find, "a2a2", this->comm);
            // send back using the constructed recv count
            mxx::all2allv(results, send_counts, this->comm).swap(results);
            BL_BENCH_END(find, "a2a2", results.size());

          } else {

            BL_BENCH_START(find);
            // keep unique keys
            ::fsc::unique(keys, sorted_input,
            		typename Base::StoreTransformedFunc(),
            		typename Base::StoreTransformedEqual());
            BL_BENCH_END(find, "uniq1", keys.size());

            BL_BENCH_START(find);
            results.reserve(keys.size());                   // TODO:  should estimate coverage.
            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
            BL_BENCH_END(find, "reserve", results.capacity() );

            size_t estimating = std::ceil(static_cast<double>(keys.size()) * 0.05);

            BL_BENCH_START(find);
            QueryProcessor::process(c, keys.begin(), keys.begin() + estimating, emplace_iter, find_element, sorted_input, pred);
            BL_BENCH_END(find, "local_find_0.1", estimating);

            BL_BENCH_START(find);
            size_t est = std::ceil((static_cast<double>(results.size()) / static_cast<double>(estimating)) * static_cast<double>(keys.size()) * 1.1f);
            results.reserve(est);
            BL_BENCH_END(find, "reserve_est", results.capacity());

            BL_BENCH_START(find);
            QueryProcessor::process(c, keys.begin() + estimating, keys.end(), emplace_iter, find_element, sorted_input, pred);
            BL_BENCH_END(find, "local_find", results.size());

            if (this->comm.rank() == 0) printf("rank %d result size %lu capacity %lu\n", this->comm.rank(), results.size(), results.capacity());

          }

          BL_BENCH_REPORT_MPI_NAMED(find, "base_hashmap:find", this->comm);

          return results;

      }


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       *
       * why this version that uses isend and irecv?  because all2all version requires all result data to be in memory.
       * this one can do it one source process at a time.
       *
       * @param keys    content will be changed and reordered.
       * @param last
       */
      template <class LocalFind, typename Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_sendrecv(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
          BL_BENCH_INIT(find);

          ::std::vector<::std::pair<Key, T> > results;

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_sendrecv", this->comm);
            return results;
          }

          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_sendrecv", this->comm);
            return results;
          }


          BL_BENCH_START(find);
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;

          ::std::vector<::std::pair<Key, T> > local_results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > local_emplace_iter(local_results);
          this->transform_input(keys);
          BL_BENCH_END(find, "transform_input", keys.size());

          BL_BENCH_START(find);
          ::fsc::unique(keys, sorted_input,
            				typename Base::StoreTransformedFunc(),
            				typename Base::StoreTransformedEqual());
          size_t num_orig_keys = keys.size();
          BL_BENCH_END(find, "unique", keys.size());


          if (this->comm.size() > 1) {

            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
            // distribute (communication part)
            std::vector<size_t> recv_counts(
            		::dsc::distribute(keys, this->key_to_rank, sorted_input, this->comm));
            BL_BENCH_END(find, "dist_query", keys.size());

            //==== reserve
            BL_BENCH_START(find);
            local_results.reserve(keys.size());
            BL_BENCH_END(find, "reserve_local", keys.size());

            //=== process queries and send results.  O(p) iterations
            BL_BENCH_START(find);
            auto recv_displs = mxx::impl::get_displacements(recv_counts);  // compute response displacements.
            int recv_from, send_to;
            size_t found[2], recved[2];
            size_t found_total = 0, recv_total = 0, reqs_total = 0, est_total = 0;
            std::vector<MPI_Request> reqs(2 * this->comm.size());
            auto start = keys.begin();
            auto end = start;

            mxx::datatype dt = mxx::get_datatype<::std::pair<Key, T> >();
            mxx::datatype size_dt = mxx::get_datatype<size_t >();

            for (int i = 0; i < this->comm.size(); ++i) {
            	// wait for previous pair of requests to complete
            	if (i > 0) {
                    MPI_Waitall(2, &reqs[2 * (i - 1)], MPI_STATUSES_IGNORE);
            	}


              send_to = (this->comm.rank() + i) % this->comm.size();    // rank to send data to

              //== get data for the dest rank
              start = keys.begin();                                   // keys for the query for the dest rank
              ::std::advance(start, recv_displs[send_to]);
              end = start;
              ::std::advance(end, recv_counts[send_to]);

              local_results.clear();
              // work on query from process i.
              found[0] = recv_counts[send_to];
              found[1] = QueryProcessor::process(c, start, end, local_emplace_iter, find_element, sorted_input, pred);
              // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);
              found_total += found[1];
              //== now send the results immediately - minimizing data usage so we need to wait for both send and recv to complete right now.

              // === recv from this rank
              recv_from = (this->comm.rank() + (this->comm.size() - i)) % this->comm.size(); // rank to recv data from


              // shift the size found around.
              MPI_Sendrecv(found, 2, size_dt.type(), send_to, i, recved, 2, size_dt.type(), recv_from, i, this->comm, MPI_STATUS_IGNORE);

              // once we have the size, now we estimate the size based on the first few iterations
              reqs_total += recved[0];

              // resize if we need to
              if (results.size() < (recv_total + recved[1])) {
            	  est_total = std::ceil((static_cast<double>(recv_total + recved[1]) / static_cast<double>(reqs_total)) * static_cast<double>(num_orig_keys) * 1.1f);
            	  if (this->comm.rank() == 0)
            		  printf("rank %d resizing results: iter %d req %lu res %lu, req_total %lu, recv_total %lu. curr result size %lu, capacity %lu.  est total %lu\n",
            				  this->comm.rank(), i, recved[0], recved[1], reqs_total, recv_total + recved[1], results.size(), results.capacity(), est_total);
            	  results.resize(est_total);
              }

              // set up receive.
              MPI_Irecv(&(results[recv_total]), recved[1], dt.type(),
                        recv_from, i + this->comm.size(), this->comm, &reqs[2 * i]);
              MPI_Isend(&(local_results[0]), found[1], dt.type(),
            		  send_to, i + this->comm.size(), this->comm, &reqs[2 * i + 1]);

//              MPI_Sendrecv(&(local_results[0]), found[1], dt.type(), send_to, i + this->comm.size(),
//            		  &(results[recv_total]), recved[1], dt.type(), recv_from, i + this->comm.size(),
//            		  this->comm, MPI_STATUS_IGNORE);

              recv_total += recved[1];

            }

            // wait for last pair of requests to complete.
            MPI_Waitall(2, &reqs[2 * (this->comm.size() - 1)], MPI_STATUSES_IGNORE);


            // now the estimate may be too large, so resize it back.
            results.resize(recv_total);


            //printf("Rank %d total find %lu\n", this->comm.rank(), total);
            BL_BENCH_END(find, "find_send", results.size());

          } else {

            BL_BENCH_START(find);
            // keep unique keys
            ::fsc::unique(keys, sorted_input,
    				typename Base::StoreTransformedFunc(),
    				typename Base::StoreTransformedEqual());
            BL_BENCH_END(find, "uniq1", keys.size());

            // memory is constrained.  find EXACT count.
            BL_BENCH_START(find);
            ::std::vector<::std::pair<Key, size_t> > count_results;
            count_results.reserve(keys.size());
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

            // count now.
            QueryProcessor::process(c, keys.begin(), keys.end(), count_emplace_iter, count_element, sorted_input, pred);
            size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                             [](size_t v, ::std::pair<Key, size_t> const & x) {
              return v + x.second;
            });
            //          for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
            //            count += it->second;
            //          }
            BL_BENCH_END(find, "local_count", count);

            BL_BENCH_START(find);
            results.reserve(count);                   // TODO:  should estimate coverage.
            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
            BL_BENCH_END(find, "reserve", results.capacity());

            BL_BENCH_START(find);
            QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, sorted_input, pred);
            BL_BENCH_END(find, "local_find", results.size());
          }

          BL_BENCH_REPORT_MPI_NAMED(find, "base_hashmap:find_sendrecv", this->comm);

          return results;

      }


//      /**
//       * @brief find elements with the specified keys in the distributed unordered_multimap.
//       *
//       * why this version that uses isend and irecv?  because all2all version requires all result data to be in memory.
//       * this one can do it one source process at a time.
//       *
//       * @param keys    content will be changed and reordered.
//       * @param last
//       */
//      template <class LocalFind, typename Predicate = ::fsc::TruePredicate>
//      ::std::vector<::std::pair<Key, T> > find_irecv(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
//          BL_BENCH_INIT(find);
//
//          ::std::vector<::std::pair<Key, T> > results;
//
//      if (::dsc::empty(keys, this->comm)) {
//        BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_irecv", this->comm);
//        return results;
//      }
//          if (this->empty()) {
//            BL_BENCH_REPORT_MPI_NAMED(find, "base_hashed_vec:find_irecv", this->comm);
//            return results;
//          }
//
//
//          BL_BENCH_START(find);
//          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
//          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
//
//          std::vector<::std::vector<::std::pair<Key, T> > > local_results(2);  // allocate 2, one foreground, 1 background.
//          std::vector<::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > > local_emplace_iters {
//            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > >(local_results[0]),
//                ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > >(local_results[1])
//          };
//
//          this->transform_input(keys);
//          BL_BENCH_END(find, "transform_input", keys.size());
//
//          if (this->comm.size() > 1) {
//
//            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
//            // distribute (communication part)
//            std::vector<size_t> recv_counts(
//            		::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
//            				typename Base::StoreTransformedFunc(),
//            				typename Base::StoreTransformedEqual()));
//
//            BL_BENCH_END(find, "dist_query", keys.size());
//
//
//            //======= local count to determine amount of memory to allocate at destination.
//            BL_BENCH_START(find);
//            ::std::vector<::std::pair<Key, size_t> > count_results;
//            size_t max_key_count = *(::std::max_element(recv_counts.begin(), recv_counts.end()));
//            count_results.reserve(max_key_count);
//            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);
//
//            std::vector<size_t> send_counts(this->comm.size(), 0);
//
//            auto start = keys.begin();
//            auto end = start;
//            size_t total = 0;
//            for (int i = 0; i < this->comm.size(); ++i) {
//              ::std::advance(end, recv_counts[i]);
//
//              // count results for process i
//              count_results.clear();
//              QueryProcessor::process(c, start, end, count_emplace_iter, count_element, sorted_input, pred);
//              send_counts[i] =
//                  ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
//                                    [](size_t v, ::std::pair<Key, size_t> const & x) {
//                return v + x.second;
//              });
//              //            for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//              //              send_counts[i] += it->second;
//              //            }
//              total += send_counts[i];
//              start = end;
//              //printf("Rank %d local count for src rank %d:  recv %d send %d\n", this->comm.rank(), i, recv_counts[i], send_counts[i]);
//            }
//            ::std::vector<::std::pair<Key, size_t> >().swap(count_results);
//            BL_BENCH_END(find, "local_count", total);
//
//
//            BL_BENCH_COLLECTIVE_START(find, "a2a_count", this->comm);
//            std::vector<size_t> resp_counts = mxx::all2all(send_counts, this->comm);  // compute counts of response to receive
//            BL_BENCH_END(find, "a2a_count", keys.size());
//
//
//            //==== reserve
//            BL_BENCH_START(find);
//            auto resp_displs = mxx::impl::get_displacements(resp_counts);  // compute response displacements.
//
//            auto resp_total = resp_displs[this->comm.size() - 1] + resp_counts[this->comm.size() - 1];
//            auto max_send_count = *(::std::max_element(send_counts.begin(), send_counts.end()));
//            results.resize(resp_total);   // allocate, not just reserve
//            local_results[0].resize(max_send_count);
//            local_results[1].resize(max_send_count);
//
//            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
//            BL_BENCH_END(find, "reserve", resp_total);
//
//            //=== process queries and send results.  O(p) iterations
//            BL_BENCH_START(find);
//            auto recv_displs = mxx::impl::get_displacements(recv_counts);  // compute response displacements.
//            int recv_from, send_to;
//            size_t found;
//            total = 0;
//            std::vector<MPI_Request> reqs(this->comm.size());
//
//            mxx::datatype dt = mxx::get_datatype<::std::pair<Key, T> >();
//
//            for (int i = 0; i < this->comm.size(); ++i) {
//              recv_from = (this->comm.rank() + (this->comm.size() - i)) % this->comm.size(); // rank to recv data from
//
//              // set up receive.
//              MPI_Irecv(&(results[resp_displs[recv_from]]), resp_counts[recv_from], dt.type(),
//                        recv_from, i, this->comm, &(reqs[i]));
//            }
//            BL_BENCH_END(find, "recv_request", this->comm.size());
//
//            BL_BENCH_START(find);
//
//            MPI_Request send_req = MPI_REQUEST_NULL;
//
//            for (int i = 0; i < this->comm.size(); ++i) {
//              send_to = (this->comm.rank() + i) % this->comm.size();    // rank to send data to
//
//              //== get data for the dest rank
//              start = keys.begin();                                   // keys for the query for the dest rank
//              ::std::advance(start, recv_displs[send_to]);
//              end = start;
//              ::std::advance(end, recv_counts[send_to]);
//
//              //local_results[i % 2].clear();
//              // work on query from process i.
//              auto oiter = local_results[i%2].begin();
//              found = QueryProcessor::process(c, start, end,
//                                              oiter, // local_emplace_iters[i % 2],
//                                              find_element, sorted_input, pred);
//              // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);
//              total += found;
//              //== now send the results immediately - minimizing data usage so we need to wait for both send and recv to complete right now.
//
//              // wait for previous to finish
//              if (send_req != MPI_REQUEST_NULL) MPI_Wait(&send_req, MPI_STATUS_IGNORE);
//
//              // then send current.
//              MPI_Isend(&(local_results[i%2][0]), found, dt.type(), send_to,
//                        i, this->comm, &send_req);
//            }
//
//            // wait for last send to finish
//            if (send_req != MPI_REQUEST_NULL) MPI_Wait(&send_req, MPI_STATUS_IGNORE);
//            BL_BENCH_END(find, "sent", total);
//
//
//            // wait for all recv requests to complete.
//            BL_BENCH_START(find);
//           MPI_Waitall(reqs.size(), &(reqs[0]), MPI_STATUSES_IGNORE);
//
//
//              // verify correct? done by comparing to previous code.
//
//            //printf("Rank %d total find %lu\n", this->comm.rank(), total);
//            BL_BENCH_END(find, "recved", results.size());
//
//          } else {
//
//            BL_BENCH_START(find);
//            // keep unique keys
//            ::fsc::unique(keys, sorted_input,
//    				typename Base::StoreTransformedFunc(),
//    				typename Base::StoreTransformedEqual());
//            BL_BENCH_END(find, "uniq1", keys.size());
//
//            // memory is constrained.  find EXACT count.
//            BL_BENCH_START(find);
//            ::std::vector<::std::pair<Key, size_t> > count_results;
//            count_results.reserve(keys.size());
//            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);
//
//            // count now.
//            QueryProcessor::process(c, keys.begin(), keys.end(), count_emplace_iter, count_element, sorted_input, pred);
//            size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
//                                             [](size_t v, ::std::pair<Key, size_t> const & x) {
//              return v + x.second;
//            });
//            //          for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//            //            count += it->second;
//            //          }
//            BL_BENCH_END(find, "local_count", count);
//
//            BL_BENCH_START(find);
//
//            results.reserve(count);                   // TODO:  should estimate coverage.
//            //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
//            BL_BENCH_END(find, "reserve", results.capacity());
//
//            BL_BENCH_START(find);
//            QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, sorted_input, pred);
//            BL_BENCH_END(find, "local_find", results.size());
//          }
//
//          BL_BENCH_REPORT_MPI_NAMED(find, "base_hashmap:find_isend", this->comm);
//
//          return results;
//
//      }


      template <class LocalFind, typename Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(LocalFind & find_element, Predicate const& pred = Predicate()) const {
          ::std::vector<::std::pair<Key, T> > results;

          if (this->local_empty()) return results;

          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

          auto keys = this->keys();

          ::std::vector<::std::pair<Key, size_t> > count_results;
          count_results.reserve(keys.size());
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

          // count now.
          QueryProcessor::process(c, keys.begin(), keys.end(), count_emplace_iter, count_element, false, pred);
          size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                           [](size_t v, ::std::pair<Key, size_t> const & x) {
            return v + x.second;
          });

          // then reserve
          results.reserve(count);                   // TODO:  should estimate coverage.

          QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, false, pred);

          return results;
      }

      template <class LocalErase, typename Predicate = ::fsc::TruePredicate>
      size_t erase(LocalErase & erase_element, ::std::vector<Key>& keys, bool sorted_input, Predicate const& pred) {
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;
          size_t before = this->c.size();
          BL_BENCH_INIT(erase);

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(erase, "base_hashed_vec:erase", this->comm);
            return 0;
          }

          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(erase, "base_hashed_vec:erase", this->comm);
            return 0;
          }

          BL_BENCH_START(erase);
          this->transform_input(keys);
          BL_BENCH_END(erase, "transform_intput", keys.size());

          if (this->comm.size() > 1) {

            BL_BENCH_START(erase);
            auto recv_counts(::dsc::distribute(keys, this->key_to_rank, sorted_input, this->comm));
            BLISS_UNUSED(recv_counts);
            BL_BENCH_END(erase, "dist_query", keys.size());

            // don't try to run unique further - have to use a set so might as well just have erase_element handle it.
            sorted_input = false;
          }

          BL_BENCH_START(erase);
          // then call local remove.
          ::fsc::unique(keys, sorted_input,
                                                  typename Base::StoreTransformedFunc(),
                                                  typename Base::StoreTransformedEqual());
          BL_BENCH_END(erase, "unique", keys.size());


          BL_BENCH_START(erase);
          // then call local remove.
          auto dummy_iter = keys.end();  // process requires a reference.
          QueryProcessor::process(this->c, keys.begin(), keys.end(), dummy_iter, erase_element, sorted_input, pred);
          BL_BENCH_END(erase, "erase", keys.size());

          BL_BENCH_REPORT_MPI_NAMED(erase, "base_hashmap:erase", this->comm);

          if (before != this->c.size()) local_changed = true;

          return before - this->c.size();
      }


      template <class LocalErase, typename Predicate>
      size_t erase(LocalErase & erase_element, Predicate const& pred) {
          size_t count = 0;

          if (this->local_empty()) return 0;


          if (!::std::is_same<Predicate, ::fsc::TruePredicate>::value) {

            auto keys = this->keys();  // already unique

            auto dummy_iter = keys.end();  // process requires a reference.
            count = QueryProcessor::process(c, keys.begin(), keys.end(), dummy_iter, erase_element, false, pred);
          } else {
            count = this->local_size();
            this->local_clear();
          }

          if (count > 0) local_changed = true;

          if (this->comm.size() > 1) this->comm.barrier();

          return count;
      }

      hashed_vec_base(const mxx::comm& _comm) : Base(_comm),
          key_to_rank(_comm.size()), local_changed(false) {}


      // ================ local overrides

      /// clears the hashed_vec
      virtual void local_clear() noexcept {
        c.clear();
      }

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      virtual void local_reserve( size_t n) {
        size_t buckets = std::ceil(static_cast<float>(n) / this->c.max_load_factor());

        if (this->c.bucket_count() < buckets) this->c.rehash(buckets);
      }



    public:

      virtual ~hashed_vec_base() {};

      /// returns the local storage.  please use sparingly.
      local_container_type& get_local_container() { return c; }

      const_iterator cbegin() const
      {
        return c.cbegin();
      }

      const_iterator cend() const {
        return c.cend();
      }

      using Base::size;
      using Base::unique_size;
      using Base::get_multiplicity;
      using Base::local_size;

      /// convert the map to a vector
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const {
        result.clear();
        if (c.empty()) return;
        result.reserve(c.size());
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(result);
        ::std::copy(c.begin(), c.end(), emplace_iter);
      }
      /// extract the unique keys of a map.
      virtual void keys(std::vector<Key> & result) const {
        result.clear();
        if (c.empty()) return;

        typename Base::template UniqueKeySetUtilityType<Key> temp(c.size());
        auto end = c.end();
        for (auto it = c.begin(); it != end; ++it) {
          temp.emplace((*it).first);
        }
        result.assign(temp.begin(), temp.end());
      }



      /**
       * @brief count elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys, bool sorted_input = false,
                                                        Predicate const& pred = Predicate() ) const {
          BL_BENCH_INIT(count);
          ::std::vector<::std::pair<Key, size_type> > results;

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(count, "base_hashed_vec:count", this->comm);
            return results;
          }

          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(count, "base_hashed_vec:count", this->comm);
            return results;
          }



          BL_BENCH_START(count);
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          this->transform_input(keys);
          BL_BENCH_END(count, "transform_intput", keys.size());


          if (this->comm.size() > 1) {

            BL_BENCH_START(count);
            // distribute (communication part)
            std::vector<size_t> recv_counts(
            		::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
            				typename Base::StoreTransformedFunc(),
            				typename Base::StoreTransformedEqual()));
            BL_BENCH_END(count, "dist_query", keys.size());


            // local count. memory utilization a potential problem.
            // do for each src proc one at a time.
            BL_BENCH_START(count);
            results.reserve(keys.size() );                   // TODO:  should estimate coverage.
            BL_BENCH_END(count, "reserve", results.capacity());

            BL_BENCH_START(count);
            auto start = keys.begin();
            auto end = start;
            for (int i = 0; i < this->comm.size(); ++i) {
              ::std::advance(end, recv_counts[i]);

              // within start-end, values are unique, so don't need to set unique to true.
              QueryProcessor::process(c, start, end, emplace_iter, count_element, sorted_input, pred);

              if (this->comm.rank() == 0)
                BL_DEBUGF("R %d added %lu results for %lu queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);

              start = end;
            }
            BL_BENCH_END(count, "local_count", results.size());

            // send back using the constructed recv count
            BL_BENCH_COLLECTIVE_START(count, "a2a2", this->comm);
            mxx::all2allv(results, recv_counts, this->comm).swap(results);
            BL_BENCH_END(count, "a2a2", results.size());
          } else {

            BL_BENCH_START(count);
            // keep unique keys
            ::fsc::unique(keys, sorted_input,
    				typename Base::StoreTransformedFunc(),
    				typename Base::StoreTransformedEqual());
            BL_BENCH_END(count, "uniq1", keys.size());


            BL_BENCH_START(count);
            results.reserve(keys.size());                   // TODO:  should estimate coverage.
            BL_BENCH_END(count, "reserve", results.capacity());


            BL_BENCH_START(count);
            // within start-end, values are unique, so don't need to set unique to true.
            QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, count_element, sorted_input, pred);
            BL_BENCH_END(count, "local_count", results.size());
          }

          BL_BENCH_REPORT_MPI_NAMED(count, "base_hashmap:count", this->comm);

          return results;

      }


      template <typename Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, size_type> > count(Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;

        if (this->local_empty()) return results;


        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size());

        QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, count_element, false, pred);

        if (this->comm.size() > 1) this->comm.barrier();
        return results;
      }



      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = ::fsc::TruePredicate>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate() ) {
          return this->erase(erase_element, keys, sorted_input, pred);
      }
      template <typename Predicate>
      size_t erase(Predicate const & pred = Predicate()) {
        return this->erase(erase_element, pred);
      }

      // ================  overrides

      // note that for each method, there is a local version of the operartion.
      // this is for use by the asynchronous version of communicator as callback for any messages received.
      /// check if empty.
      virtual bool local_empty() const {
        return this->c.empty();
      }

      /// get size of local container
      virtual size_t local_size() const {
        if (this->comm.rank() == 0) printf("rank %d hashmap_base local size %lu\n", this->comm.rank(), this->c.size());

        return this->c.size();
      }

      /// get size of local container
      virtual size_t local_unique_size() const {
        return this->local_size();
      }



  };



  /**
   * @brief  distributed unordered multimap following std unordered multimap's interface.
   * @details   This class is modeled after the std::unordered_multimap.
   *         it does not have all the methods of std::unordered_multimap.  Whatever methods that are present considers the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Iterators are assumed to be local rather than distributed, so methods that returns iterators are not provided.
   *         as an alternative, vectors are returned.
   *         methods that accept iterators as input assume that the input data is local.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  template <typename> class MapParams,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_multimap_hashvec : public hashed_vec_base<Key, T, ::fsc::hashed_vecmap, MapParams, Alloc> {

    protected:
      using Base = hashed_vec_base<Key, T, ::fsc::hashed_vecmap, MapParams, Alloc>;


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

      struct LocalFind {
        // unfiltered.
        template<class DB, typename Query, class OutputIter>
        size_t operator()(DB &db, Query const &v, OutputIter &output) const {
            auto range = db.equal_range(v);

            output = ::std::copy(range.first, range.second, output);  // tons faster to emplace - almost 3x faster
            return db.count(v);
        }
        // filtered element-wise.
        template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
        size_t operator()(DB &db, Query const &v, OutputIter &output,
                          Predicate const& pred) const {
            auto range = db.equal_range(v);

            // add the output entry.
            size_t count = 0;

            if (pred(range.first, range.second)) {
              for (auto it2 = range.first; it2 != range.second; ++it2) {

                if (pred(*it2)) {
                  *output = *it2;
                  ++output;
                  ++count;
                }
              }
            }
            return count;
        }
        // no filter by range AND elemenet for now.
      } find_element;


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_insert(InputIterator first, InputIterator last) {
        BL_BENCH_INIT(local_insert);

          if (first == last) return 0;

          size_t before = this->c.size();

          BL_BENCH_START(local_insert);
          this->c.insert(first, last);
          BL_BENCH_END(local_insert, "insert", this->c.size());

          if (this->c.size() != before) this->local_changed = true;


          BL_BENCH_REPORT_MPI_NAMED(local_insert, "hashvec:local_insert", this->comm);

          //          c.insert(first, last);  // mem usage?
          return this->c.size() - before;
      }



      /**
       * @brief insert new elements in the distributed unordered_multimap.  example use: stop inserting if more than x entries.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      size_t local_insert(InputIterator first, InputIterator last, Predicate const &pred) {

          auto new_end = std::partition(first, last, pred);
//
//          this->local_reserve(c.size() + ::std::distance(first, last));   // before branching, because reserve calls collective "empty()"
//
          if (first == last) return 0;

          size_t before = this->c.size();

          this->local_insert(first, new_end);
//
//          for (auto it = first; it != last; ++it) {
//            if (pred(*it)) c.emplace(*it);
//          }
//
          if (this->c.size() != before) this->local_changed = true;

          return this->c.size() - before;
      }


    public:


      unordered_multimap_hashvec(const mxx::comm& _comm) : Base(_comm) { }

      virtual ~unordered_multimap_hashvec() {}

      using Base::count;
      using Base::unique_size;


      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_overlap(::std::vector<Key>& keys, bool sorted_input = false,
                                               Predicate const& pred = Predicate()) const {
          return Base::find_overlap(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_collective(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find_a2a(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_sendrecv(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find_sendrecv(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }

      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = ::fsc::TruePredicate>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) {
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;
          size_t before = this->c.size();
          BL_BENCH_INIT(erase);

          if (::dsc::empty(keys, this->comm)) {
            BL_BENCH_REPORT_MPI_NAMED(erase, "hashvec:erase", this->comm);
            return 0;
          }

          if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(erase, "hashvec:erase", this->comm);
            return 0;
          }

          BL_BENCH_START(erase);
          this->transform_input(keys);
          BL_BENCH_END(erase, "transform_intput", keys.size());

          if (this->comm.size() > 1) {

            BL_BENCH_START(erase);
            auto recv_counts(::dsc::distribute(keys, this->key_to_rank, sorted_input, this->comm));
            BLISS_UNUSED(recv_counts);
            BL_BENCH_END(erase, "dist_query", keys.size());

            // don't try to run unique further - have to use a set so might as well just have erase_element handle it.
            sorted_input = false;
          }

          BL_BENCH_START(erase);
          // then call local remove.
          ::fsc::unique(keys, sorted_input,
                                                  typename Base::StoreTransformedFunc(),
                                                  typename Base::StoreTransformedEqual());
          BL_BENCH_END(erase, "unique", keys.size());


          BL_BENCH_START(erase);
          if (std::is_same<Predicate, ::fsc::TruePredicate>::value)
            this->c.erase(keys.begin(), keys.end());
          else
            this->c.erase(keys.begin(), keys.end(), pred);
          BL_BENCH_END(erase, "erase", keys.size());

          BL_BENCH_REPORT_MPI_NAMED(erase, "hashvec:erase", this->comm);

          if (before != this->c.size()) this->local_changed = true;

          return before - this->c.size();
      }



      template <typename Predicate>
      size_t erase(Predicate const & pred = Predicate()) {
        if (std::is_same<Predicate, ::fsc::TruePredicate>::value) {
          size_t s = this->c.size();
          this->c.clear();
          return s;
        }
        else
          return this->c.erase(pred);
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = ::fsc::TruePredicate>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(insert);

        if (::dsc::empty(input, this->comm)) {
          BL_BENCH_REPORT_MPI_NAMED(insert, "hashvec:insert", this->comm);
          return 0;
        }


        BL_BENCH_START(insert);
        this->transform_input(input);
        BL_BENCH_END(insert, "transform_intput", input.size());


        //        printf("r %d key size %lu, val size %lu, pair size %lu, tuple size %lu\n", this->comm.rank(), sizeof(Key), sizeof(T), sizeof(::std::pair<Key, T>), sizeof(::std::tuple<Key, T>));
        //        count_unique(input);
        //        count_unique(bucketing(input, this->key_to_rank, this->comm));

        // communication part
        if (this->comm.size() > 1) {
          BL_BENCH_START(insert);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          auto recv_counts = ::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm);
          BLISS_UNUSED(recv_counts);
          BL_BENCH_END(insert, "dist_data", input.size());

        }

        // reserve should be done outside of insert.

        BL_BENCH_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, ::fsc::TruePredicate>::value)
          count = this->local_insert(input.begin(), input.end(), pred);
        else
          count = this->local_insert(input.begin(), input.end());
        BL_BENCH_END(insert, "insert", this->c.size());

        BL_BENCH_REPORT_MPI_NAMED(insert, "hashvec:insert", this->comm);
        return count;

      }



      /// access the current the multiplicity.  only multimap needs to override this.
      virtual float get_multiplicity() const {
        // multimaps would add a collective function to change the multiplicity
       // if (this->comm.rank() == 0) printf("rank %d vec get_multiplicity called\n", this->comm.rank());

        // one approach is to add up the number of repeats for the key of each entry, then divide by total count.
        //  sum(count per key) / c.size.
        // problem with this approach is that for unordered map, to get the count for a key is essentially O(count), so we get quadratic time.
        // The approach is VERY SLOW for large repeat count.  - (0.0078125 human: 52 sec, synth: FOREVER.)

        // a second approach is to count the number of unique key then divide the map size by that.
        //  c.size / #unique.  requires unique set
        // To find unique set, we take each bucket, copy to vector, sort it, and then count unique.
        // This is precise, and is faster than the approach above.  (0.0078125 humanMaM MaMaÿ?: 54 sec.  synth: 57sec.)
        // but the n log(n) sort still grows with the duplicate count
//        size_t unique_count = this->unique_size();
//        float multiplicity = 1.0f;
//        if (unique_count > 0) {
//          // local unique
//          multiplicity =
//              static_cast<float>(this->size()) /
//              static_cast<float>(unique_count);
//        }
//          if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity max %lu\n", this->comm.rank(), this->c.get_max_multiplicity());
//          if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity mean %f\n", this->comm.rank(), this->c.get_mean_multiplicity());

          BL_BENCH_INIT(multiplicity);

          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "container", this->c.size());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "unique", this->c.unique_size());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "max_multi", this->c.get_max_multiplicity());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "mean_multi", this->c.get_mean_multiplicity());

          BL_BENCH_REPORT_MPI_NAMED(multiplicity, "hashvec:multiplicity", this->comm);


        return this->c.get_mean_multiplicity();
//        return 1.0f;
      }

      /// get the size of unique keys in the current local container.
      virtual size_t local_unique_size() const {
        return this->c.unique_size();
      }

  };


  /**
   * @brief  distributed unordered multimap following std unordered multimap's interface.
   * @details   This class is modeled after the std::unordered_multimap.
   *         it does not have all the methods of std::unordered_multimap.  Whatever methods that are present considers the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Iterators are assumed to be local rather than distributed, so methods that returns iterators are not provided.
   *         as an alternative, vectors are returned.
   *         methods that accept iterators as input assume that the input data is local.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  template <typename> class MapParams,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_multimap_vec : public hashed_vec_base<Key, T, ::fsc::unordered_vecmap, MapParams, Alloc> {

    protected:
      using Base = hashed_vec_base<Key, T, ::fsc::unordered_vecmap, MapParams, Alloc>;


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

      struct LocalFind {
        // unfiltered.
        template<class DB, typename Query, class OutputIter>
        size_t operator()(DB &db, Query const &v, OutputIter &output) const {
            auto range = db.equal_range_value_only(v);

            output = ::std::copy(range.first, range.second, output);  // tons faster to emplace - almost 3x faster
            return db.count(v);
        }
        // filtered element-wise.
        template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
        size_t operator()(DB &db, Query const &v, OutputIter &output,
                          Predicate const& pred) const {
            auto range = db.equal_range_value_only(v);

            // add the output entry.
            size_t count = 0;

            typename std::iterator_traits<decltype(range.first)>::value_type w;

            if (pred(range.first, range.second)) {
              for (auto it2 = range.first; it2 != range.second; ++it2) {
                w = *it2;
                if (pred(w)) {
                  *output = w;
                  ++output;
                  ++count;
                }
              }
            }
            return count;
        }
        // no filter by range AND elemenet for now.
      } find_element;

      // need to deal with fact that vector's iterators are invalidated.
      struct LocalErase {
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &) {
              size_t before = db.size();
              db.erase(v);
              return before - db.size();
          }
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &,
                            Predicate const & pred) {
              auto range = db.equal_range(v);  // get the range as values

              // check range first.  then erase.  only removed iterators are invalidated.
              // order of remaining elements preserved.  true for map/multimap, unordered or not.
              size_t count = 0;
              if (pred(range.first, range.second)) { // operator to decide if range matches.
                count = db.erase(v, pred);
              }

              return count;
          }
          // no filter by range AND elemenet for now.
      } erase_element;

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_insert(InputIterator first, InputIterator last) {
        BL_BENCH_INIT(local_insert);

        size_t dist = ::std::distance(first, last);

        BL_BENCH_START(local_insert);
          this->local_reserve(this->c.size() + dist);  // before branching, because reserve calls collective "empty()"
          BL_BENCH_END(local_insert, "reserve", this->c.size() + dist);

          if (first == last) return 0;


          BL_BENCH_START(local_insert);
//          ::std::sort(first, last, typename Base::Base::StoreTransformedLess());
          ::std::sort(first, last, ::fsc::TransformedComparator<Key, ::std::less, ::bliss::kmer::transform::identity>());

          BL_BENCH_END(local_insert, "sort", dist);


          size_t before = this->c.size();

          BL_BENCH_START(local_insert);
          this->c.insert_sorted(first, last);
          BL_BENCH_END(local_insert, "insert_sorted", this->c.size());

          if (this->c.size() != before) this->local_changed = true;


          BL_BENCH_REPORT_MPI_NAMED(local_insert, "vecmap:local_insert", this->comm);

          //          c.insert(first, last);  // mem usage?
          return this->c.size() - before;
      }

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      virtual void local_reserve( size_t n) {
        if (!this->empty()) { // nothing is available.   so assume multiplicity of 1.
          // has some data. so try to reserve buckets based on local information
          float multiplicity = static_cast<float>(this->c.size()) / static_cast<float>(this->c.unique_size());

          n = std::ceil(n / multiplicity);
        }

        this->Base::local_reserve(n);
      }

    public:


      unordered_multimap_vec(const mxx::comm& _comm) : Base(_comm) { }

      virtual ~unordered_multimap_vec() {}

      using Base::count;
      using Base::unique_size;


      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_overlap(::std::vector<Key>& keys, bool sorted_input = false,
                                               Predicate const& pred = Predicate()) const {
          /*


          // DEBUG

          auto keys2 = keys;
          auto result_a2a =  Base::find_a2a(find_element, keys2, sorted_input, pred);
          auto result = Base::find(find_element, keys, sorted_input, pred);


          // DEBUG
          if (result.size() != result_a2a.size()) {
              throw ::std::logic_error("ERROR: not same size");

          }

          ::std::stable_sort(result.begin(), result.end(), typename Base::Base::TransformedLess());
          ::std::stable_sort(result_a2a.begin(), result_a2a.end(), typename Base::Base::TransformedLess());

          typename Base::Base::TransformedEqual eq;
          for (int i = 0; i < result_a2a.size(); ++i) {
            if (!eq(result[i], result_a2a[i])) {
              printf("failing at %d:  result: %s, result_a2a: %s\n", i, result[i].first.toAlphabetString().c_str(), result_a2a[i].first.toAlphabetString().c_str());
              printf("  before   %d:  result: %s, result_a2a: %s\n", i-1, result[i-1].first.toAlphabetString().c_str(), result_a2a[i-1].first.toAlphabetString().c_str());
              printf("  after    %d:  result: %s, result_a2a: %s\n", i+1, result[i+1].first.toAlphabetString().c_str(), result_a2a[i+1].first.toAlphabetString().c_str());
              throw ::std::logic_error("ERROR: not same.");
            }
          }
//          bool same = ::std::equal(result.begin(), result.end(), result_a2a.begin(), typename Base::Base::TransformedEqual());
//          if (!same) throw ::std::logic_error("ERROR: not same.");

          return result;
           */
          return Base::find_overlap(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_collective(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find_a2a(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_sendrecv(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find_sendrecv(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }

      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = ::fsc::TruePredicate>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate() ) {
          return Base::erase(erase_element, keys, sorted_input, pred);
      }


      template <typename Predicate>
      size_t erase(Predicate const & pred = Predicate()) {
        return Base::erase(erase_element, pred);
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = ::fsc::TruePredicate>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(insert);

        if (::dsc::empty(input, this->comm)) {
          BL_BENCH_REPORT_MPI_NAMED(insert, "vecmap:insert", this->comm);
          return 0;
        }


        BL_BENCH_START(insert);
        this->transform_input(input);
        BL_BENCH_END(insert, "transform_intput", input.size());


        //        printf("r %d key size %lu, val size %lu, pair size %lu, tuple size %lu\n", this->comm.rank(), sizeof(Key), sizeof(T), sizeof(::std::pair<Key, T>), sizeof(::std::tuple<Key, T>));
        //        count_unique(input);
        //        count_unique(bucketing(input, this->key_to_rank, this->comm));

        // communication part
        if (this->comm.size() > 1) {
          BL_BENCH_START(insert);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          auto recv_counts = ::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm);
          BLISS_UNUSED(recv_counts);
          BL_BENCH_END(insert, "dist_data", input.size());

        }

        // reserve should be done outside of insert.

        BL_BENCH_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, ::fsc::TruePredicate>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->local_insert(input.begin(), input.end());
        BL_BENCH_END(insert, "insert", this->c.size());

        BL_BENCH_REPORT_MPI_NAMED(insert, "vecmap:insert", this->comm);
        return count;

      }



      /// access the current the multiplicity.  only multimap needs to override this.
      virtual float get_multiplicity() const {
        // multimaps would add a collective function to change the multiplicity
       // if (this->comm.rank() == 0) printf("rank %d vec get_multiplicity called\n", this->comm.rank());

        // one approach is to add up the number of repeats for the key of each entry, then divide by total count.
        //  sum(count per key) / c.size.
        // problem with this approach is that for unordered map, to get the count for a key is essentially O(count), so we get quadratic time.
        // The approach is VERY SLOW for large repeat count.  - (0.0078125 human: 52 sec, synth: FOREVER.)

        // a second approach is to count the number of unique key then divide the map size by that.
        //  c.size / #unique.  requires unique set
        // To find unique set, we take each bucket, copy to vector, sort it, and then count unique.
        // This is precise, and is faster than the approach above.  (0.0078125 humanMaM MaMaÿ?: 54 sec.  synth: 57sec.)
        // but the n log(n) sort still grows with the duplicate count
//        size_t unique_count = this->unique_size();
//        float multiplicity = 1.0f;
//        if (unique_count > 0) {
//          // local unique
//          multiplicity =
//              static_cast<float>(this->size()) /
//              static_cast<float>(unique_count);
//        }
//          if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity max %lu\n", this->comm.rank(), this->c.get_max_multiplicity());
//          if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity mean %f\n", this->comm.rank(), this->c.get_mean_multiplicity());

          BL_BENCH_INIT(multiplicity);

          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "container", this->c.size());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "unique", this->c.unique_size());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "max_multi", this->c.get_max_multiplicity());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "mean_multi", this->c.get_mean_multiplicity());

          BL_BENCH_REPORT_MPI_NAMED(multiplicity, "compact_vecmap:multiplicity", this->comm);


        return this->c.get_mean_multiplicity();
//        return 1.0f;
      }

      /// get the size of unique keys in the current local container.
      virtual size_t local_unique_size() const {
        return this->c.unique_size();
      }

  };


  /**
   * @brief  distributed unordered multimap following std unordered multimap's interface.
   * @details   This class is modeled after the std::unordered_multimap.
   *         it does not have all the methods of std::unordered_multimap.  Whatever methods that are present considers the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Iterators are assumed to be local rather than distributed, so methods that returns iterators are not provided.
   *         as an alternative, vectors are returned.
   *         methods that accept iterators as input assume that the input data is local.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  template <typename> class MapParams,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_multimap_compact_vec : public hashed_vec_base<Key, T, ::fsc::unordered_compact_vecmap, MapParams, Alloc> {
    protected:
      using Base = hashed_vec_base<Key, T, ::fsc::unordered_compact_vecmap, MapParams, Alloc>;


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

      struct LocalFind {
        // unfiltered.
        template<class DB, typename Query, class OutputIter>
        size_t operator()(DB &db, Query const &v, OutputIter &output) const {
            auto range = db.equal_range(v);

            //              // DEBUG
            //              bool equal = true;
            //              typename Base::TransformedEqual eq;
            //              for (auto it = range.first; it != range.second; ++it) {
            //                equal &= eq(*it, v);
            //              }
            //              if (!equal)  printf("ERROR: NOT EQUAL!");

            output = ::std::copy(range.first, range.second, output);  // tons faster to emplace - almost 3x faster
            return db.count(v);
        }
        // filtered element-wise.
        template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
        size_t operator()(DB &db, Query const &v, OutputIter &output,
                          Predicate const& pred) const {
            auto range = db.equal_range(v);

            // add the output entry.
            size_t count = 0;
            if (pred(range.first, range.second)) {
              for (auto it2 = range.first; it2 != range.second; ++it2) {
                if (pred(*it2)) {
                  *output = *it2;
                  ++output;
                  ++count;
                }
              }
            }
            return count;
        }
        // no filter by range AND elemenet for now.
      } find_element;

      struct LocalErase {
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &) {
              size_t before = db.size();
              db.erase(v);
              return before - db.size();
          }
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter, class Predicate = ::fsc::TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &,
                            Predicate const & pred) {
              auto range = db.equal_range(v);

              // check range first.  then erase.  only removed iterators are invalidated.
              // order of remaining elements preserved.  true for map/multimap, unordered or not.
              size_t count = 0;
              if (pred(range.first, range.second)) { // operator to decide if range matches.
                db.erase(v, pred);
              }

              return count;
          }
          // no filter by range AND elemenet for now.
      } erase_element;


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_insert(InputIterator first, InputIterator last) {
        BL_BENCH_INIT(local_insert);

        size_t dist = ::std::distance(first, last);

        BL_BENCH_START(local_insert);
          this->local_reserve(this->c.size() + dist);  // before branching, because reserve calls collective "empty()"
          BL_BENCH_END(local_insert, "reserve", this->c.size() + dist);

          if (first == last) return 0;

          BL_BENCH_START(local_insert);
//          ::std::sort(first, last, typename Base::Base::StoreTransformedLess());
          ::std::sort(first, last, ::fsc::TransformedComparator<Key, ::std::less, ::bliss::kmer::transform::identity>());
          BL_BENCH_END(local_insert, "sort", dist);


          size_t before = this->c.size();

          BL_BENCH_START(local_insert);
          this->c.insert_sorted(first, last);
          BL_BENCH_END(local_insert, "insert_sorted", this->c.size());

          if (this->c.size() != before) this->local_changed = true;


          BL_BENCH_REPORT_MPI_NAMED(local_insert, "compactvecmap:local_insert", this->comm);

          //          c.insert(first, last);  // mem usage?
          return this->c.size() - before;
      }


      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      virtual void local_reserve( size_t n) {
        if (!this->empty()) {// nothing is available.   so assume multiplicity of 1.
          // has some data. so try to reserve buckets based on local information
          float multiplicity = static_cast<float>(this->c.size()) / static_cast<float>(this->c.unique_size());

          n = std::ceil(n / multiplicity);
        }

        this->Base::local_reserve(n);

      }


    public:


      unordered_multimap_compact_vec(const mxx::comm& _comm) : Base(_comm) {}

      virtual ~unordered_multimap_compact_vec() {}

      using Base::count;
      using Base::unique_size;
      using Base::size;
      using Base::local_size;


      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_overlap(::std::vector<Key>& keys, bool sorted_input = false,
                                               Predicate const& pred = Predicate()) const {
          /*
          // DEBUG

          auto keys2 = keys;
          auto result_a2a =  Base::find_a2a(find_element, keys2, sorted_input, pred);
          auto result = Base::find(find_element, keys, sorted_input, pred);


          // DEBUG
          if (result.size() != result_a2a.size()) {
              throw ::std::logic_error("ERROR: not same size");

          }

          ::std::stable_sort(result.begin(), result.end(), typename Base::Base::TransformedLess());
          ::std::stable_sort(result_a2a.begin(), result_a2a.end(), typename Base::Base::TransformedLess());

          typename Base::Base::TransformedEqual eq;
          for (int i = 0; i < result_a2a.size(); ++i) {
            if (!eq(result[i], result_a2a[i])) {
              printf("failing at %d:  result: %s, result_a2a: %s\n", i, result[i].first.toAlphabetString().c_str(), result_a2a[i].first.toAlphabetString().c_str());
              printf("  before   %d:  result: %s, result_a2a: %s\n", i-1, result[i-1].first.toAlphabetString().c_str(), result_a2a[i-1].first.toAlphabetString().c_str());
              printf("  after    %d:  result: %s, result_a2a: %s\n", i+1, result[i+1].first.toAlphabetString().c_str(), result_a2a[i+1].first.toAlphabetString().c_str());
              throw ::std::logic_error("ERROR: not same.");
            }
          }
//          bool same = ::std::equal(result.begin(), result.end(), result_a2a.begin(), typename Base::Base::TransformedEqual());
//          if (!same) throw ::std::logic_error("ERROR: not same.");

          return result;
           */
          return Base::find_overlap(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_collective(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find_a2a(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_sendrecv(::std::vector<Key>& keys, bool sorted_input = false,
                                                          Predicate const& pred = Predicate()) const {
          return Base::find_sendrecv(find_element, keys, sorted_input, pred);
      }


      template <class Predicate = ::fsc::TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }

      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = ::fsc::TruePredicate>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate() ) {
          return Base::erase(erase_element, keys, sorted_input, pred);
      }


      template <typename Predicate>
      size_t erase(Predicate const & pred = Predicate()) {
        return Base::erase(erase_element, pred);
      }




      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = ::fsc::TruePredicate>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(insert);

        if (::dsc::empty(input, this->comm)) {
          BL_BENCH_REPORT_MPI_NAMED(insert, "compact_vecmap:insert", this->comm);
          return 0;
        }


        BL_BENCH_START(insert);
        this->transform_input(input);
        BL_BENCH_END(insert, "transform_intput", input.size());

        //        printf("r %d key size %lu, val size %lu, pair size %lu, tuple size %lu\n", this->comm.rank(), sizeof(Key), sizeof(T), sizeof(::std::pair<Key, T>), sizeof(::std::tuple<Key, T>));
        //        count_unique(input);
        //        count_unique(bucketing(input, this->key_to_rank, this->comm));

        // communication part
        if (this->comm.size() > 1) {
          BL_BENCH_START(insert);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          auto recv_counts = ::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm);
          BLISS_UNUSED(recv_counts);
          BL_BENCH_END(insert, "dist_data", input.size());
        }


        BL_BENCH_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, ::fsc::TruePredicate>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->local_insert(input.begin(), input.end());
        BL_BENCH_END(insert, "insert", this->c.size());

        BL_BENCH_REPORT_MPI_NAMED(insert, "compact_vecmap:insert", this->comm);
        return count;

      }


      /// update the multiplicity.  only multimap needs to do this.
      virtual float get_multiplicity() const {
          BL_BENCH_INIT(multiplicity);

          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "container", this->c.size());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "unique", this->c.unique_size());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "max_multi", this->c.get_max_multiplicity());
          BL_BENCH_START(multiplicity);
          BL_BENCH_END(multiplicity, "mean_multi", this->c.get_mean_multiplicity());

          BL_BENCH_REPORT_MPI_NAMED(multiplicity, "compact_vecmap:multiplicity", this->comm);

        //        if (this->comm.rank() == 0) printf("rank %d compactvec get_multiplicity called\n", this->comm.rank());
//
//        if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity min %lu\n", this->comm.rank(), this->c.get_min_multiplicity());
//        if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity max %lu\n", this->comm.rank(), this->c.get_max_multiplicity());
//        if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity mean %f\n", this->comm.rank(), this->c.get_mean_multiplicity());
//        if (this->comm.rank() == 0) printf("rank %d compactvec local multiplicity stdev %f\n", this->comm.rank(), this->c.get_stdev_multiplicity());



        // one approach is to add up the number of repeats for the key of each entry, then divide by total count.
        //  sum(count per key) / c.size.
        // problem with this approach is that for unordered map, to get the count for a key is essentially O(count), so we get quadratic time.
        // The approach is VERY SLOW for large repeat count.  - (0.0078125 human: 52 sec, synth: FOREVER.)

        // a second approach is to count the number of unique key then divide the map size by that.
        //  c.size / #unique.  requires unique set
        // To find unique set, we take each bucket, copy to vector, sort it, and then count unique.
        // This is precise, and is faster than the approach above.  (0.0078125 human: 54 sec.  synth: 57sec.)
        // but the n log(n) sort still grows with the duplicate count
//        size_t unique_count = this->unique_size();
//        size_t total = this->size();
//        if (this->comm.rank() == 0) printf("rank %d compactvec unique_count is %lu, size is %lu\n", this->comm.rank(), unique_count, total);
//        float multiplicity = 1.0f;
//        if (unique_count > 0) {
//          // local unique
//          multiplicity =
//              static_cast<float>(total) /
//              static_cast<float>(unique_count);
//        }


        return this->c.get_mean_multiplicity();
      }


      /// get the size of unique keys in the current local container.
      virtual size_t local_unique_size() const {
        return this->c.unique_size();
      }
  };



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_HASHED_VEC_HPP
