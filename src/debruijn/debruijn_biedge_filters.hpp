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
#include "utils/sys_utils.hpp"

#include <mxx/reduction.hpp>

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
          {
            BL_BENCH_LOOP_START(compute_edge_frequency, 0);
            BL_BENCH_LOOP_START(compute_edge_frequency, 1);
            size_t total = 0;
            ::std::vector<K1merType> temp;
            for (auto x : file_data) {
              BL_BENCH_LOOP_RESUME(compute_edge_frequency, 0);
              temp.clear();
              // the parser needs to collectively parse the records.
              ::bliss::io::KmerFileHelper::template parse_file_data<K1merParser, SeqParser, SeqIterType>(x, temp, comm);
              total += temp.size();
              BL_BENCH_LOOP_PAUSE(compute_edge_frequency, 0);

              BL_BENCH_LOOP_RESUME(compute_edge_frequency, 1);
              counter.insert(temp);  // this distributes the counts according to k-mer hash.
              // TODO: build from k+2 mer, so that overlap region does not become an issue for fasta files.
              BL_BENCH_LOOP_PAUSE(compute_edge_frequency, 1);
            }
            BL_BENCH_LOOP_END(compute_edge_frequency, 0, "parse_file", total);
            BL_BENCH_LOOP_END(compute_edge_frequency, 1, "idx_insert", counter.local_size());
          }  // erase temp


          BL_BENCH_START(compute_edge_frequency);
          if ((lower_thresh > 0) || (upper_thresh < ::std::numeric_limits<typename CounterType::mapped_type>::max())) {
            // ======= filter k+1-mers
            // now filter out the low frequency (and high frequency) ones.
            counter.erase([&lower_thresh, &upper_thresh](typename CounterType::value_type const & x) {
              return ((x.second < lower_thresh) || (x.second >= upper_thresh));
            });
            counter.reserve(0);  // compact the counter.
          }
          BL_BENCH_END(compute_edge_frequency, "filter counter", counter.local_size());

          BL_BENCH_REPORT_MPI_NAMED(compute_edge_frequency, "compute_edge_frequency", comm);
        }


        /// convenience function to generate k+1mer frequency map.
        template <template <typename> class SeqParser, template <typename, template <typename> class> class SeqIter, typename CounterType>
        void compute_edge_frequency_incremental(::std::vector<::bliss::io::file_data> const & file_data, CounterType & counter, mxx::comm const & comm,
                                    typename CounterType::mapped_type const & lower_thresh = 0,
                                    typename CounterType::mapped_type const & upper_thresh = ::std::numeric_limits<typename CounterType::mapped_type>::max()) {


        	BL_BENCH_INIT(compute_edge_frequency);

        	// ========  count the k+1-mers.
        	if (comm.rank() == 0) printf("Compute Edge Frequency\n");

        	// k+1-mer count map.  note that it should use the same hash function as Index.
        	using K1merType = typename CounterType::key_type;

        	using CharIterType = typename ::bliss::io::file_data::const_iterator;
        	using SeqParserType = SeqParser<CharIterType>;
        	using SeqIterType = SeqIter<CharIterType, SeqParser>;
        	using KmerParser = ::bliss::index::kmer::KmerParser<K1merType>;
        	using Iter = typename ::bliss::iterator::ContainerConcatenatingIterator<SeqIterType, KmerParser>;

        	BL_BENCH_START(compute_edge_frequency);

        	// estimate the largest amount of memory to use.
        	unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);

        	// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
        	size_t block_size = (free_mem / (8 * sizeof(typename KmerParser::value_type)));  // number of elements that can be held in freemem

        	if (comm.rank() == 0) std::cout << "estimate num elements=" << block_size << ", value_type size=" <<
        			sizeof(typename KmerParser::value_type) << " bytes" << std::endl;

        	::std::vector<typename KmerParser::value_type> temp2;
        	temp2.reserve(block_size);
        	::fsc::back_emplace_iterator<std::vector<typename KmerParser::value_type> > emplace_iter(temp2);
        	// TESTING END;
        	BL_BENCH_END(compute_edge_frequency, "reserve", block_size);

        	BL_BENCH_LOOP_START(compute_edge_frequency, 0);  // for init
        	BL_BENCH_LOOP_START(compute_edge_frequency, 1);  // for parse
        	BL_BENCH_LOOP_START(compute_edge_frequency, 2);  // for insert
        	size_t count = 0, i;
        	bool all_done = false;

        	for (auto x : file_data) {

        		// initialization
        		BL_BENCH_LOOP_RESUME(compute_edge_frequency, 0);

        		// not reusing the SeqParser in loader.  instead, reinitializing one.
        		SeqParserType seq_parser;
        		seq_parser.init_parser(x.in_mem_cbegin(), x.parent_range_bytes, x.in_mem_range_bytes, x.getRange());

        		//==  and wrap the chunk inside an iterator that emits Reads.
        		SeqIterType seqs_start(seq_parser, x.cbegin(), x.in_mem_cend(), x.getRange().start);
        		SeqIterType seqs_end(x.in_mem_cend());

        		//== sequence parser type
        		KmerParser kmer_parser(x.valid_range_bytes);

        		// now make the concatenated iterators
        		Iter start(kmer_parser, seqs_start, seqs_end);
        		Iter endd(kmer_parser, seqs_end);

        		BL_BENCH_LOOP_PAUSE(compute_edge_frequency, 0);

        		//=== copy into index incrementally
        		while (! all_done) {
            		temp2.clear();

        			//== process the chunk of data
        			BL_BENCH_LOOP_RESUME(compute_edge_frequency, 1);

        			for (i = 0; (i < block_size) && (start != endd); ++i) {
        				*emplace_iter = *start;
        				++start;
        				++emplace_iter;
        			}
        			count += i;

        			BL_BENCH_LOOP_PAUSE(compute_edge_frequency, 1);


        			BL_BENCH_LOOP_RESUME(compute_edge_frequency, 2);
        			all_done = i == 0;
        			all_done = mxx::all_of(all_done, comm);

        			counter.insert(temp2);
        			BL_BENCH_LOOP_PAUSE(compute_edge_frequency, 2);

        		}
        	}
        	BL_BENCH_LOOP_END(compute_edge_frequency, 0, "setup", temp2.capacity());  // for init
        	BL_BENCH_LOOP_END(compute_edge_frequency, 1, "parse", count);  // for parse
        	BL_BENCH_LOOP_END(compute_edge_frequency, 2, "count", counter.local_size());  // for insert

        	size_t total = counter.size();
        	if (comm.rank() == 0) printf("KmerCounter DONE: total size after insert/rehash is %lu\n", total);


			BL_BENCH_START(compute_edge_frequency);
			if ((lower_thresh > 0) || (upper_thresh < ::std::numeric_limits<typename CounterType::mapped_type>::max())) {
				// ======= filter k+1-mers
				// now filter out the low frequency (and high frequency) ones.
				counter.erase([&lower_thresh, &upper_thresh](typename CounterType::value_type const & x) {
				  return ((x.second < lower_thresh) || (x.second >= upper_thresh));
				});
				counter.reserve(0);  // compact the counter.
			}
			BL_BENCH_COLLECTIVE_END(compute_edge_frequency, "filter", counter.local_size(), comm);

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

          //BL_BENCH_INIT(filter_biedge_by_frequency);

          if (comm.rank() == 0) {
        	  std::cout << "SIZES compact simple biedge size: " << sizeof(::bliss::debruijn::biedge::compact_simple_biedge) <<
        			  " kmer size " << sizeof(KmerType) << " node size " << sizeof(std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge>) << std::endl;
          }

          // define k1mer type
          using K1merType = typename Counter::key_type;

          size_t changed_count = 0;
          bool edge_changed = false;

          // create a mask.
          // BL_BENCH_START(filter_biedge_by_frequency);

          // local storage of query results.
          ::std::vector<K1merType> query;
          query.reserve(biedges.size() + 1);  // preallocate space

          // no k1mertype here...
          ::std::vector<unsigned char> remote_exists;
          remote_exists.reserve(query.capacity());

          //BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "init", biedges.size(), comm);


          // do left and right together, in batches of step_size.
          //BL_BENCH_START(filter_biedge_by_frequency);

          K1merType k1;
          ::bliss::kmer::transform::lex_less<K1merType> canonical;

          //===== add the very first left query (for reads that are split between partitions).
          //  when read is split between partitions, the second half gets the first node at offset of 1 from the partition start.
          //  we need to query for this edge.  this is at i == 0.
          if (biedges.size() > 0) {
            k1 = canonical(::bliss::debruijn::biedge::get_in_edge_k1mer(biedges[0]));
            query.emplace_back(k1);
          }

            //===== populate right query.  only right is needed since left edge of the next node is the same.
            // NOTE: do both and not checking local counts from prev iteration.
            // next time, do check and update the bit vec here and later.
            for (size_t j = 0; j < biedges.size(); ++j) {
              // get left as canonical.
              k1 = canonical(::bliss::debruijn::biedge::get_out_edge_k1mer(biedges[j]));  // if empty edge, would not be a valid k+1 mer in counter.

              // always insert.
              query.emplace_back(k1);
            }

            //===== query right and insert into local map.
            // one to one because mxx::bucketing is stable.
            k1mer_counter.exists(query, false).swap(remote_exists);
            assert(remote_exists.size() == query.size());

            //====== NOW go through k2mers again and modify the biedges based on frequency.
            // we check left and right edges here so that all edges are consistent.
            for (size_t j = 0, k = 0; j < biedges.size(); ++j, ++k) {
              // check each for query result
            	edge_changed = false;  // DEBUG ONLY

              if (remote_exists[k] == 0) {
                // did not find, so clear the in edge (upper 4 bits of the byte).
                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0x0F));

                biedges[j].second.getDataRef()[0] &= 0x0F;
              }

              if (remote_exists[k+1]== 0) {
                // did not find, so clear the out edge (lower 4 bits of the byte).
                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0xF0));

                biedges[j].second.getDataRef()[0] &= 0xF0;
              }

              if (edge_changed) ++changed_count;

            } // end scan of current step to create the bit vector.

          //BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "query_filter", changed_count, comm);

//          std::cout << "rank " << comm.rank() << " transformed " << changed_count << " by kmer frequency. " << std::endl;


          //BL_BENCH_REPORT_MPI_NAMED(filter_biedge_by_frequency, "filter_biedge_by_frequency", comm);

        }


        // return the number of attempted insertion into the index
        template <typename KmerType, typename InputIter, typename Counter, typename OutputIndex, typename EdgeOutputIter>
        size_t freq_filter_insert_biedges(InputIter start, InputIter end,
        		Counter const & k1mer_counter, const size_t block_size,
        		OutputIndex & idx, EdgeOutputIter edge_start,
        		mxx::comm const & comm) {

        	// some verification
          using NodeType = typename std::iterator_traits<InputIter>::value_type;
          // define k1mer type
          using K1merType = typename Counter::key_type;

          static_assert(std::is_same<NodeType,
        		  std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> >::value,
        		  "ERROR: Input iterator DOES NOT have <KmerType,compact_simple_biedge> as value type.");
          static_assert(K1merType::size > 0, "counter kmer type should not be zero-mers.  possible?");
          static_assert(KmerType::size + 1 == K1merType::size,
                        "counter kmer type should have length 1 less than the input kmer type to filter for implicit debruijn chain nodes");


          static_assert(std::is_same<KmerType,
        		  typename OutputIndex::key_type >::value,
        		  "ERROR: Output index DOES NOT have KmerType as key type.");
          static_assert(std::is_same<typename std::iterator_traits<EdgeOutputIter>::value_type,
        		  ::bliss::debruijn::biedge::compact_simple_biedge >::value,
        		  "ERROR: Edge Output iterator DOES NOT have compact_simple_biedge as value type.");

          BL_BENCH_INIT(freq_filter_insert_biedges);

          if (comm.rank() == 0) {
        	  std::cout << "SIZES compact simple biedge size: " << sizeof(::bliss::debruijn::biedge::compact_simple_biedge) <<
        			  " kmer size " << sizeof(KmerType) << " node size " << sizeof(NodeType) <<
					  " k1mer size " << sizeof(K1merType) << std::endl;
          }

          size_t changed_count = 0;
          bool edge_changed = false;

          // create a mask.
          BL_BENCH_START(freq_filter_insert_biedges);

          // local storage of query results.
          ::std::vector<K1merType> query;
          query.reserve(block_size + 1);  // preallocate space

          // no k1mertype here...
          ::std::vector<unsigned char> remote_exists;
          remote_exists.reserve(block_size + 1);

          // storage for index insertion.
          ::std::vector<NodeType> nodes;
          nodes.reserve(block_size);  // preallocate space

          // intermediate stuff.
          ::bliss::kmer::transform::lex_less<K1merType> canonical;
          bool all_done = false;
          size_t i;
          auto it1 = start, it2 = start;
          K1merType k1;
          NodeType node;

          size_t total0 = 0, total1 = 0, total2 = 0, total3 = 0;
          //size_t before = idx.local_size();
          BL_BENCH_END(freq_filter_insert_biedges, "setup", block_size);


          BL_BENCH_LOOP_START(freq_filter_insert_biedges, 0);  // for reset
          BL_BENCH_LOOP_START(freq_filter_insert_biedges, 1);  // for parse
		  BL_BENCH_LOOP_START(freq_filter_insert_biedges, 2);  // for query
		  BL_BENCH_LOOP_START(freq_filter_insert_biedges, 3);  // for filter
		  BL_BENCH_LOOP_START(freq_filter_insert_biedges, 4);  // for insert

		  while (! all_done) {   // loop until all nodeds are done

			  BL_BENCH_LOOP_RESUME(freq_filter_insert_biedges, 0);
			  total0 += query.size();
			  query.clear();
			  remote_exists.clear();
			  nodes.clear();
			  BL_BENCH_LOOP_PAUSE(freq_filter_insert_biedges, 0);


			  //=====  get the k1mers for query
	          //===== add the very first left query (for reads that are split between partitions).
	          //  when read is split between partitions, the second half gets the first node at offset of 1 from the partition start.
	          //  we need to query for this edge.  this is at i == 0.
			  BL_BENCH_LOOP_RESUME(freq_filter_insert_biedges, 1);
	          if (it1 != end) {
	            k1 = canonical(::bliss::debruijn::biedge::get_in_edge_k1mer(*it1));
	            query.emplace_back(k1);
	          }
			    //===== populate right query.  only right is needed since left edge of the next node is the same.
	            // NOTE: do both and not checking local counts from prev iteration.
	            // next time, do check and update the bit vec here and later.
			  for (i = 0; (i < block_size) && (it1 != end); ++i, ++it1) {
	              // get left as canonical.
	              k1 = canonical(::bliss::debruijn::biedge::get_out_edge_k1mer(*it1));  // if empty edge, would not be a valid k+1 mer in counter.

	              // always insert.
	              query.emplace_back(k1);
			  }
			  total1 += query.size();
			  BL_BENCH_LOOP_PAUSE(freq_filter_insert_biedges, 1);

              //===== query right and insert into local map.
              // one to one because mxx::bucketing is stable.
			  BL_BENCH_LOOP_RESUME(freq_filter_insert_biedges, 2);
  			  //std::cout << "rank " << comm.rank() << " query size " << query_size << " i " << i << std::endl;
  			  assert((query.size() == 0) || (query.size() == i + 1));
  			  k1mer_counter.exists(query, false).swap(remote_exists);
              assert(remote_exists.size() == query.size());
			  total2 += remote_exists.size();
			  BL_BENCH_LOOP_PAUSE(freq_filter_insert_biedges, 2);


              //====== NOW go through k2mers again and modify the biedges based on frequency.
              // we check left and right edges here so that all edges are consistent.
			  // co-iterate between input and output.
			  BL_BENCH_LOOP_RESUME(freq_filter_insert_biedges, 3);
			  for (i = 0; (i < block_size) && (it2 != end); ++i, ++it2) {

				  edge_changed = false;  // DEBUG ONLY

				  node = *it2;

				  //==== transform the node's edge based on existence (i.e by frequency)
	              if (remote_exists[i] == 0) {
	                // did not find, so clear the in edge (upper 4 bits of the byte).
	                edge_changed = (node.second.getDataRef()[0] != (node.second.getDataRef()[0] & 0x0F));

	                node.second.getDataRef()[0] &= 0x0F;
	              }

	              if (remote_exists[i+1]== 0) {
	                // did not find, so clear the out edge (lower 4 bits of the byte).
	                edge_changed = (node.second.getDataRef()[0] != (node.second.getDataRef()[0] & 0xF0));

	                node.second.getDataRef()[0] &= 0xF0;
	              }

	              if (edge_changed) ++changed_count;

				  //==== unconditionally insert biedge into EdgeOutputIter
				  *edge_start = node.second;
				  ++edge_start;

				  //==== conditionally insert into nodes.
	              if (node.second.getDataRef()[0] != 0)
	            	  nodes.emplace_back(node);

			  }
			  total3 += nodes.size();
			  BL_BENCH_LOOP_PAUSE(freq_filter_insert_biedges, 3);

			  //==== index insert
			  BL_BENCH_LOOP_RESUME(freq_filter_insert_biedges, 4);
			  all_done = (i == 0);
			  all_done = mxx::all_of(all_done, comm);

			  idx.insert(nodes);
			  BL_BENCH_LOOP_PAUSE(freq_filter_insert_biedges, 4);
		  }

          BL_BENCH_LOOP_END(freq_filter_insert_biedges, 0, "reset", total0);  // for reset  // amount cleared
          BL_BENCH_LOOP_END(freq_filter_insert_biedges, 1, "parse", total1);  // for parse  // num queries
		  BL_BENCH_LOOP_END(freq_filter_insert_biedges, 2, "query", total2);  // for query  // num query results
		  BL_BENCH_LOOP_END(freq_filter_insert_biedges, 3, "filter", total3);  // for filter  // insert attempts
		  BL_BENCH_LOOP_END(freq_filter_insert_biedges, 4, "insert", idx.local_size());  // for insert  // curr index size


          BL_BENCH_REPORT_MPI_NAMED(freq_filter_insert_biedges, "filter_insert_biedges", comm);

          return total3;

        }


//
//
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
//        void transform_biedges_by_frequency3(std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > & biedges,
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
////          ::std::vector<K1merType> query2;
////          ::std::vector<K1merType> query3;
////          query2.reserve(std::min(step_size, biedges.size()));  // preallocate space
////          query3.reserve(std::min(step_size, biedges.size()));  // preallocate space
//
//          // no k1mertype here...
//          ::std::vector<unsigned char> remote_exists;
////			::std::vector<std::pair<K1merType, unsigned char> > remote_exists;
//			remote_exists.reserve(query.capacity());
//
////          ::std::vector<std::pair<K1merType, size_t> > remote_counts;
////          remote_counts.reserve(query.capacity());
//
//          BL_BENCH_COLLECTIVE_END(filter_biedge_by_frequency, "init", biedges.size(), comm);
//
//          // do left and right together, in batches of step_size.
//          BL_BENCH_START(filter_biedge_by_frequency);
//
//          size_t jmin, jmax;
//          K1merType k1;
//          ::bliss::kmer::transform::lex_less<K1merType> canonical;
//
//          for (size_t i = 0; i < iterations; ++i) {
//            // biedges indices for this iteration.
//            jmin = std::min(biedges.size(), i * step_size);
//            jmax = std::min(biedges.size(), jmin + step_size);
//
//            // clear query from previous iteration.
//            query.clear();
////            query2.clear();
////            query3.clear();
//            remote_exists.clear();
////            remote_counts.clear();
//
//            //===== add the very first left query (for reads that are split between partitions).
//            //  when read is split between partitions, the second half gets the first node at offset of 1 from the partition start.
//            //  we need to query for this edge.  this is at i == 0.
//            if (jmin != jmax) {
//            	k1 = canonical(::bliss::debruijn::biedge::get_in_edge_k1mer(biedges[jmin]));
//            	query.emplace_back(k1);
////            	query2.emplace_back(k1);
////            	query3.emplace_back(k1);
//            }
//
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
////              query2.emplace_back(k1);
////              query3.emplace_back(k1);
//
//            }
//
//            size_t query_size = query.size();
////            if (query_size > 0) {
////            	K1merType first = query.front();
////            	K1merType last = query.back();
////            	std::cout << "rank " << comm.rank() << "n edges=" << query.size() << " first" << first << " last=" << last << std::endl;
////            } else {
////            	std::cout << "rank " << comm.rank() << "no edges." << std::endl;
////            }
//
//            //===== query right and insert into local map.
//            // one to one because mxx::bucketing is stable.
//            k1mer_counter.exists(query, false).swap(remote_exists);
//
////            k1mer_counter.template count<false>(query2, false).swap(remote_counts);
//
////			bool same = true;
////			bool same2 = true;
////			for (size_t l = 0; l < query_size; ++l) {
////				same = ((remote_exists[l].second == 1) == (remote_counts[l].second > 0));
////				if (!same) {
////					std::cout << "rank " << comm.rank() << " result differ at " << l << ": " << (remote_exists[l].second == 1 ? "y" : "n") << ", " <<  remote_counts[l].second << std::endl;
////				}
////			}
////			for (size_t l = 0; l < query_size; ++l) {
////				same = (remote_exists[l].first == remote_counts[l].first);
////				if (!same) {
////					std::cout << "rank " << comm.rank() << " kmer ordering differ at " << l << ": " << remote_exists[l].first << " VS " << remote_counts[l].first << std::endl;
////				}
////			}
////			for (size_t l = 0; l < query_size; ++l) {
////				same = (remote_exists[l].first == query3[l]);
////				if (!same) {
////					std::cout << "rank " << comm.rank() << " bucket stability differ at " << l << ": " << remote_exists[l].first << " VS query " << query3[l] << std::endl;
////				}
////			}
////			for (size_t l = 0; l < query.size(); ++l) {
////				same2 = query[l] == query2[l];
////				if (!same2) {
////					std::cout << "rank " << comm.rank() << " distribute consistency differ at " << l << ": " << query[l] << ",  " << query2[l] << std::endl;
////				}
////			}
////			std::cout << "rank " << comm.rank() << " query size " << query_size << " exists " << remote_exists.size() << " counts " << remote_counts.size() << " same? " << (same ? "y" : "n") << std::endl;
////
////			assert(remote_counts.size() == query_size);
//            assert(remote_exists.size() == query_size);
//
////            if (remote_exists.size() > 0)
////            std::cout << "rank " << comm.rank() << " first edge " << remote_exists.front() <<
////            		" last edge " << remote_exists.back() << " n edges = " << remote_exists.size() << std::endl;
//
//            //====== NOW go through k2mers again and modify the biedges based on frequency.
//            // we check left and right edges here so that all edges are consistent.
//            size_t k = 0;
//            for (size_t j = jmin; j < jmax; ++j, ++k) {
//              // check each for query result
//            	edge_changed = false;  // DEBUG ONLY
//
//              if (remote_exists[k] == 0) {
//                // did not find, so clear the in edge (upper 4 bits of the byte).
//                edge_changed = (biedges[j].second.getDataRef()[0] != (biedges[j].second.getDataRef()[0] & 0x0F));
////                if (edge_changed) std::cout << "rank " << comm.rank() << " missing count locally for in edge " << k1 << std::endl;
//
//                biedges[j].second.getDataRef()[0] &= 0x0F;
//              }
//
//              if (remote_exists[k+1]== 0) {
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
//         *
//         *            VERSION USING COUNTMAP's COUNT FUNCTION, REQUIRES A LOCAL MAP FOR NOW.
//         *
//         * @param biedges  <kmer, biedge> nodes.  may be ordered, e.g. same as file reading order.
//         * @param k1mer_counter   edge frequencies, distributed hash table.
//         * @param comm     MPI communicator
//         */
//        template <typename KmerType, typename Counter>
//        void transform_biedges_by_frequency_2(std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > & biedges,
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
