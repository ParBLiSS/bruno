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
 * debruijn_biedge_filters.hpp
 *
 *  Created on: Sept 26, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_BIEDGE_FILTERS_HPP_
#define DEBRUIJN_BIEDGE_FILTERS_HPP_

#include <vector>
#include "iterators/zip_iterator.hpp"
#include "io/kmer_parser.hpp"
#include "utils/benchmark_utils.hpp"
#include "utils/filter_utils.hpp"

namespace bliss {
  namespace debruijn {
    namespace biedge {
      namespace filter {


        /// extract the biedge elements from the nodes.
        template <typename KmerType>
        std::vector<::bliss::debruijn::biedge::compact_simple_biedge>
        extract_biedges(std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > const & biedges) {
          std::vector<::bliss::debruijn::biedge::compact_simple_biedge> results(biedges.size());

          ::std::transform(biedges.begin(), biedges.end(), results.begin(), ::bliss::debruijn::biedge::transform::extract_biedge());

          return results;
        }

        /// merge kmer and biedges back into a contiguous memory data structure.
        template <typename KmerType>
        std::vector<::std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> >
        reconstruct_nodes(std::vector<KmerType> const & kmers, std::vector<::bliss::debruijn::biedge::compact_simple_biedge> const & biedges) {
          std::vector<::std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > results(kmers.size());

          ::std::transform(kmers.begin(), kmers.end(), biedges.begin(), results.begin(), ::bliss::debruijn::biedge::transform::merge_kmer_biedge());

          return results;
        }

        /// select only non-isolated nodes.
        template <typename KmerType>
        void filter_isolated_nodes(std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > & biedges) {
          // first partition.

          auto new_end = std::partition(biedges.begin(), biedges.end(), [](std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> const & x){
            return x.second.getData()[0] != 0;  // at least one edge char is not 0, i.e. kmer is not isolated.
          });

          // next erase
          biedges.erase(new_end, biedges.end());
        }



        template <typename KmerType>
        std::vector<::std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> >
        reconstruct_filter_nodes(std::vector<KmerType> const & kmers, std::vector<::bliss::debruijn::biedge::compact_simple_biedge> const & biedges) {

          using merge_iterator = ::bliss::iterator::ZipIterator<typename std::vector<KmerType>::const_iterator,
                typename std::vector<::bliss::debruijn::biedge::compact_simple_biedge>::const_iterator>;


          std::vector<::std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > results;
          results.reserve(kmers.size());
          ::fsc::back_emplace_iterator<std::vector<::std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > > emplacer(results);


          ::std::copy_if(merge_iterator(kmers.begin(), biedges.begin()), merge_iterator(kmers.end(), biedges.end()),
                         emplacer, [](std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> const & x){
                           return x.second.getData()[0] != 0;  // at least one edge char is not 0, i.e. kmer is not isolated.
                         });

          return results;
        }

        /// reconstruct filtered nodes by replacing the edges, and copy only if filtered
        template <typename KmerType>
        std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> >
        reconstruct_filter_nodes(
        		std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > const & nodes,
				std::vector<::bliss::debruijn::biedge::compact_simple_biedge> const & biedges)
		{

          assert(nodes.size() == biedges.size());

          std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > results;
          results.reserve(nodes.size());
          ::fsc::back_emplace_iterator<std::vector<::std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > > emplacer(results);

          auto n_it = nodes.begin();
          auto e_it = biedges.begin();
          for (; n_it != nodes.end(); ++n_it, ++e_it ) {
        	  if ((*e_it).getData()[0] != 0) {
        		  // edge is not empty.
        		  results.emplace_back((*n_it).first, *e_it);
        	  }
          }
          return results;
        }

#define USE_MPI
#if defined(USE_MPI)

        /// convenience function to generate k+1mer frequency map.
        template <template <typename> class SeqParser, template <typename, template <typename> class> class SeqIterType, typename CounterType>
        void compute_edge_frequency(::std::vector<::bliss::io::file_data> const & file_data, CounterType & counter, mxx::comm const & comm,
                                    typename CounterType::mapped_type const & lower_thresh = 0,
                                    typename CounterType::mapped_type const & upper_thresh = ::std::numeric_limits<typename CounterType::mapped_type>::max()) {

          // k+1-mer count map.  note that it should use the same hash function as Index.
          using K1merType = typename CounterType::key_type;
          using K1merParser = ::bliss::index::kmer::KmerParser<K1merType>;

          BL_BENCH_INIT(compute_edge_frequency);

          // ========  count the k+1-mers.
          BL_BENCH_START(compute_edge_frequency);
          {
            ::std::vector<K1merType> temp;
            for (auto x : file_data) {
              temp.clear();
              // the parser needs to collectively parse the records.
              ::bliss::io::KmerFileHelper::template parse_file_data<K1merParser, SeqParser, SeqIterType>(x, temp, comm);
              counter.insert(temp);  // this distributes the counts according to k-mer hash.
              // TODO: build from k+2 mer, so that overlap region does not become an issue for fasta files.
            }
          }  // erase temp
          BL_BENCH_COLLECTIVE_END(compute_edge_frequency, "count edge", counter.local_size(), comm);


          BL_BENCH_START(compute_edge_frequency);
          if ((lower_thresh > 0) || (upper_thresh < ::std::numeric_limits<typename CounterType::mapped_type>::max())) {
            // ======= filter k+1-mers
            // now filter out the low frequency (and high frequency) ones.
            counter.erase([&lower_thresh, &upper_thresh](typename CounterType::value_type const & x) {
              return ((x.second < lower_thresh) || (x.second >= upper_thresh));
            });
            counter.reserve(0);  // compact the counter.
          }
          BL_BENCH_COLLECTIVE_END(compute_edge_frequency, "filter counter", counter.local_size(), comm);

          BL_BENCH_REPORT_MPI_NAMED(compute_edge_frequency, "compute_edge_frequency", comm);
        }


        /**
         * @brief transform biedge by edge frequency
         * @details   compact_simple_biedge are zeroed out (char by char) if the associated edge has invalid frequency (not present in k1mer_counter.
         *            the biedge is modified.
         *
         *            process:   get out edges of the <kmer, biedge> nodes, then query for frequency
         *                       frequencies stored in local count index (done in blocks)
         *                       the <kmer, biedge> nodes are then scanned, if edge does not exist in local count index,
         *                       then corresponding compact_simple_biedge character is zeroed.
         *
         *            order of biedges does not change.  no entry of biedges are deleted here.
         * @param biedges  <kmer, biedge> nodes.  may be ordered, e.g. same as file reading order.
         * @param k1mer_counter   edge frequencies, distributed hash table.
         * @param comm     MPI communicator
         */
        template <typename KmerType, typename Counter>
        void transform_biedges_by_frequency(std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > & biedges,
                                                              Counter const & k1mer_counter, mxx::comm const & comm) {
          static_assert(Counter::key_type::size > 0, "counter kmer type should not be zero-mers.  possible?");
          static_assert(KmerType::size + 1 == Counter::key_type::size,
                        "counter kmer type should have length 1 less than the input kmer type to filter for implicit debruijn chain nodes");

          BL_BENCH_INIT(filter_biedge_by_frequency);

          // define k1mer type
          using K1merType = typename Counter::key_type;

          size_t changed_count = 0;
          bool edge_changed = false;

          // create a mask.
          BL_BENCH_START(filter_biedge_by_frequency);

          // calculate the total number of iterations globally
          size_t step_size = 8000000;   // this is arbitrary choice.
          size_t iterations = (biedges.size() + step_size - 1) / step_size;
          iterations = ::mxx::allreduce(iterations, [](size_t const & x, size_t const & y){
            return std::max(x, y);
          }, comm);

          // local storage of query results.
          ::std::vector<K1merType> query;
          query.reserve(std::min(step_size, biedges.size()));  // preallocate space
          ::std::vector<unsigned char> remote_counts;
          remote_counts.reserve(query.capacity());
          BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "init", biedges.size(), comm);

          // do left and right together, in batches of step_size.
          BL_BENCH_START(filter_biedge_by_frequency);

          size_t jmin, jmax;
          K1merType k1;
          ::bliss::kmer::transform::lex_less<K1merType> canonical;

          for (size_t i = 0; i < iterations; ++i) {
            // biedges indices for this iteration.
            jmin = std::min(biedges.size(), i * step_size);
            jmax = std::min(biedges.size(), jmin + step_size);

            // clear query from previous iteration.
            query.clear();
            remote_counts.clear();

            //===== add the very first left query (for reads that are split between partitions).
            //  when read is split between partitions, the second half gets the first node at offset of 1 from the partition start.
            //  we need to query for this edge.  this is at i == 0.
            if (jmin != jmax) query.emplace_back(canonical(::bliss::debruijn::biedge::get_in_edge_k1mer(biedges[jmin])));

            //===== populate right query.  only right is needed since left edge of the next node is the same.
            // NOTE: do both and not checking local counts from prev iteration.
            // next time, do check and update the bit vec here and later.
            for (size_t j = jmin; j < jmax; ++j) {
              // get left as canonical.
              k1 = canonical(::bliss::debruijn::biedge::get_out_edge_k1mer(biedges[j]));  // if empty edge, would not be a valid k+1 mer in counter.

              // always insert.
              query.emplace_back(k1);
            }

            //===== query right and insert into local map.
            // one to one because mxx::bucketing is stable.
            size_t query_size = query.size();
	    k1mer_counter.exists(query, false).swap(remote_counts);
            std::cout << "rank " << comm.rank() << " query size " << query_size << " result " << remote_counts.size() << std::endl;

            assert(remote_counts.size() == query_size);


            //====== NOW go through k2mers again and modify the biedges based on frequency.
            // we check left and right edges here so that all edges are consistent.
            size_t k = 0;
            for (size_t j = jmin; j < jmax; ++j, ++k) {
              // check each for query result
            	edge_changed = false;  // DEBUG ONLY

              if (remote_counts[k] == 0) {
                // did not find, so clear the in edge (upper 4 bits of the byte).
                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0x0F));
//                if (edge_changed) std::cout << "rank " << comm.rank() << " missing count locally for in edge " << k1 << std::endl;

                biedges[j].second.getDataRef()[0] &= 0x0F;

              }

              if (remote_counts[k+1] == 0) {
                // did not find, so clear the out edge (lower 4 bits of the byte).
                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0xF0));
//                if (edge_changed) std::cout << "rank " << comm.rank() << " missing count locally for out edge " << k1 << std::endl;

                biedges[j].second.getDataRef()[0] &= 0xF0;
              }

              if (edge_changed) ++changed_count;

            } // end scan of current step to create the bit vector.

          }  // end iterations
          BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "query_filter", changed_count, comm);

//          std::cout << "rank " << comm.rank() << " transformed " << changed_count << " by kmer frequency. " << std::endl;


          BL_BENCH_REPORT_MPI_NAMED(filter_biedge_by_frequency, "filter_biedge_by_frequency", comm);

        }


//        /**
//         * @brief transform biedge by edge frequency
//         * @details   compact_simple_biedge are zeroed out (char by char) if the associated edge has invalid frequency (not present in k1mer_counter.
//         *            the biedge is modified.
//         *
//         *            process:   get out edges of the <kmer, biedge> nodes, then query for frequency
//         *                       frequencies stored in local count index (done in blocks)
//         *                       the <kmer, biedge> nodes are then scanned, if edge does not exist in local count index,
//         *                       then corresponding compact_simple_biedge character is zeroed.
//         *
//         *            order of biedges does not change.  no entry of biedges are deleted here.
//         * @param biedges  <kmer, biedge> nodes.  may be ordered, e.g. same as file reading order.
//         * @param k1mer_counter   edge frequencies, distributed hash table.
//         * @param comm     MPI communicator
//         */
//        template <typename KmerType, typename Counter>
//        void transform_biedges_by_frequency(std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > & biedges,
//                                                              Counter const & k1mer_counter, mxx::comm const & comm) {
//          static_assert(Counter::key_type::size > 0, "counter kmer type should not be zero-mers.  possible?");
//          static_assert(KmerType::size + 1 == Counter::key_type::size,
//                        "counter kmer type should have length 1 less than the input kmer type to filter for implicit debruijn chain nodes");
//
//          BL_BENCH_INIT(filter_biedge_by_frequency);
//
//          // define k1mer type
//          using K1merType = typename Counter::key_type;
//
//          size_t changed_count = 0;
//          bool edge_changed = false;
//
//          // create a mask.
//          BL_BENCH_START(filter_biedge_by_frequency);
//
//          // calculate the total number of iterations globally
//          size_t step_size = 8000000;   // this is arbitrary choice.
//          size_t iterations = (biedges.size() + step_size - 1) / step_size;
//          iterations = ::mxx::allreduce(iterations, [](size_t const & x, size_t const & y){
//            return std::max(x, y);
//          }, comm);
//
//          // local storage of query results.
//          ::std::vector<K1merType> query;
//          query.reserve(std::min(step_size, biedges.size()));  // preallocate space
//
//          using K1merCountMap = typename Counter::local_container_type;
//          using K1merCountMapIter = typename K1merCountMap::const_iterator;
//          K1merCountMap local_counts;  // keep it small, approximately size of query.
//          local_counts.resize(query.capacity());
//
//          K1merCountMapIter count_iter;
//          BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "init", biedges.size(), comm);
//
//          // do left and right together, in batches of step_size.
//          BL_BENCH_START(filter_biedge_by_frequency);
//
//          size_t jmin, jmax;
//          K1merType k1;
//          ::bliss::kmer::transform::lex_less<K1merType> canonical;
//
//          std::vector<::std::pair<K1merType, size_t> > remote_counts;
//          remote_counts.reserve(query.capacity());
//          for (size_t i = 0; i < iterations; ++i) {
//            // biedges indices for this iteration.
//            jmin = std::min(biedges.size(), i * step_size);
//            jmax = std::min(biedges.size(), jmin + step_size);
//
//            // clear query from previous iteration.
//            query.clear();
//
//            //===== add the very first left query (for reads that are split between partitions).
//            //  when read is split between partitions, the second half gets the first node at offset of 1 from the partition start.
//            //  we need to query for this edge.  this is at i == 0.
//            if (jmin != jmax) query.emplace_back(canonical(::bliss::debruijn::biedge::get_in_edge_k1mer(biedges[jmin])));
//
//            //===== populate right query.  only right is needed since left edge of the next node is the same.
//            // NOTE: do both and not checking local counts from prev iteration.
//            // next time, do check and update the bit vec here and later.
//            for (size_t j = jmin; j < jmax; ++j) {
//              // get left as canonical.
//              k1 = canonical(::bliss::debruijn::biedge::get_out_edge_k1mer(biedges[j]));  // if empty edge, would not be a valid k+1 mer in counter.
//
//              // always insert.
//              query.emplace_back(k1);
//            }
//
////            for (auto q : query) {
////              std::cout << "query k1mer: " << bliss::utils::KmerUtils::toASCIIString(q) << ::std::endl;
////            }
//
//            //===== query right and insert into local map.
//            {
//              local_counts.clear();
//              // not doing unique since unique requires a high AVERAGE repeat rate to make sense.  this reduces running time and space usage.
//              // using find would reduce return data size.
//              auto remote_counts = k1mer_counter.template count<false>(query);
//
//              // results have canonical kmers.  insert directly into local_counts.
//              // insert into local count map.
//              for (auto it = remote_counts.begin(); it != remote_counts.end(); ++it) {
//                if ((*it).second > 0) local_counts.insert(*it);
//              }
//            }
//
////            std::cout << "rank " << comm.rank() << " iter " << i << " query count " << query.size() << " result " << local_counts.size() << std::endl;
//
//            //====== NOW go through k2mers again and modify the biedges based on frequency.
//            // we check left and right edges here so that all edges are consistent.
//            for (size_t j = jmin; j < jmax; ++j) {
//              // check each for query result
//              edge_changed = false;  // DEBUG ONLY
//
//              // get left as canonical.
//              k1 = canonical(::bliss::debruijn::biedge::get_in_edge_k1mer(biedges[j]));
//              // check local count
//              count_iter = local_counts.find(k1);
//              if (count_iter == local_counts.end()) {
//                // did not find, so clear the in edge (upper 4 bits of the byte).
//                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0x0F));
////                if (edge_changed) std::cout << "rank " << comm.rank() << " missing count locally for in edge " << k1 << std::endl;
//
//                biedges[j].second.getDataRef()[0] &= 0x0F;
//
//              }
//
//              // get right as canonical.
//              k1 = canonical(::bliss::debruijn::biedge::get_out_edge_k1mer(biedges[j]));
//              // check local count
//              count_iter = local_counts.find(k1);
//              if (count_iter == local_counts.end()) {
//                // did not find, so clear the out edge (lower 4 bits of the byte).
//                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0xF0));
////                if (edge_changed) std::cout << "rank " << comm.rank() << " missing count locally for out edge " << k1 << std::endl;
//
//                biedges[j].second.getDataRef()[0] &= 0xF0;
//              }
//
//              if (edge_changed) ++changed_count;
//
//            } // end scan of current step to create the bit vector.
//
//          }  // end iterations
//          BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "query_filter", changed_count, comm);
//
////          std::cout << "rank " << comm.rank() << " transformed " << changed_count << " by kmer frequency. " << std::endl;
//
//
//          BL_BENCH_REPORT_MPI_NAMED(filter_biedge_by_frequency, "filter_biedge_by_frequency", comm);
//
//        }

#endif


      } //namespace filter
    } // namespace biedge
  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_GRAPH_FILTERS_HPP_
