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
 * compact_debuijn_graph_block_construct.cpp
 *
 * goal: 	unified compact_debruijn_graph_refactor.cpp and compact_debruijn_graph_block_construct.cpp into one file.
 * 			minimize memory usage during kmer counting and node filtering phase.
 *
 *  Created on: Aug 17, 2015
 *      Author: yongchao
 *
 *  Rewrote: June 10, 2016
 *      Author: tony pan
 *
 */

#include "bliss-config.hpp"

#include <unistd.h>  // get hostname

#include <functional>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>  // for system("pause");
#include <fstream>  // ofstream
#include <utility>  // std::declval

#include "utils/logging.h"

#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"

#include "io/sequence_iterator.hpp"
#include "io/filtered_sequence_iterator.hpp"
//#include "io/sequence_id_iterator.hpp"
#include "io/kmer_file_helper.hpp"
#include "io/kmer_parser.hpp"

#include "iterators/transform_iterator.hpp"

//#include "common/kmer_iterators.hpp"
//#include "iterators/zip_iterator.hpp"
//#include "iterators/unzip_iterator.hpp"
//#include "index/quality_score_iterator.hpp"

#include "index/kmer_index.hpp"
#include "index/kmer_hash.hpp"

#include "debruijn/debruijn_common.hpp"
#include "debruijn/debruijn_biedge_loader.hpp"
#include "debruijn/debruijn_biedge_filters.hpp"

#include "debruijn/debruijn_graph.hpp"
#include "debruijn/debruijn_graph_node.hpp"
#include "debruijn/debruijn_graph_map.hpp"
#include "debruijn/debruijn_graph_loader.hpp"
#include "debruijn/debruijn_graph_filters.hpp"
#include "debruijn/debruijn_graph_operations.hpp"

#include "debruijn/debruijn_chain_filters.hpp"
#include "debruijn/debruijn_chain_node.hpp"
#include "debruijn/debruijn_chain_graph.hpp"
#include "debruijn/debruijn_chain_operations.hpp"


#include "debruijn/debruijn_stats.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"
#include "utils/tclap_utils.hpp"
#include "utils/sys_utils.hpp"


#include "mxx/env.hpp"
#include "mxx/comm.hpp"

#if (pDNA == 16)
	using Alphabet = bliss::common::DNA16;
	#define bDNA 4
#elif (pDNA == 5)
	using Alphabet = bliss::common::DNA5;
	#define bDNA 3
#elif (pDNA == 4)
	using Alphabet = bliss::common::DNA;
	#define bDNA 2
#else
	#define bDNA 0
#endif

// need to handle different bit depths.
//using CountType = uint16_t;  // by default, 16 bit to ensure most frequencies are captured.

#if defined(MIN_MEM)

	#if defined(pK)
	// compute the data type from bDNA and pK.  DEVIATION FROM PREV:  uint8_t default instead of uint16_t
		#if (((bDNA * pK) & 0x3F) == 0)
			using KmerWordType = uint64_t;
		#elif (((bDNA * pK) & 0x1F) == 0)
			using KmerWordType = uint32_t;
		#elif (((bDNA * pK) & 0x0F) == 0)
			using KmerWordType = uint16_t;
		#else
			using KmerWordType = uint8_t;
		#endif	

		using KmerType = bliss::common::Kmer<pK, Alphabet, KmerWordType>;
	#else   // pK not defined.  assume 31.  compute data type.
		#if (pDNA == 16)
			using KmerWordType = uint64_t;
		#elif (pDNA == 5)
			using KmerWordType = uint32_t;
		#elif (pDNA == 4)
			using KmerWordType = uint64_t;
		#endif
		using KmerType = bliss::common::Kmer<31, Alphabet, KmerWordType>;
	#endif
		using CountType = uint8_t;  // by default, 16 bit to ensure most frequencies are captured.
#else   // MIN_MEM not defined.  use uint64_t.  this is same as compact_debruijn_graph_refactor.cpp

	using KmerWordType = uint64_t;  // matches system architecture.

	#if defined(pK)
		using KmerType = bliss::common::Kmer<pK, Alphabet, KmerWordType>;
	#else
		using KmerType = bliss::common::Kmer<31, Alphabet, KmerWordType>;
	#endif
		using CountType = uint16_t;  // by default, 16 bit to ensure most frequencies are captured.
#endif

#define FASTA 1
#define FASTQ 0

//============== index input file format
#if (pPARSER == FASTA)
#define FileParser ::bliss::io::FASTAParser
#elif (pPARSER == FASTQ)
#define FileParser ::bliss::io::FASTQParser
#endif


using EdgeEncoding = Alphabet;

// sequence iterator for use in constructing the dbg.  here we throw away reads containing N.  technically, the splitting sequence iterator works here too,
// but the frequency calculated would include kmers in reads containing N
template <typename Iterator, template <typename> class SeqParser>
using SplitSeqIterType = bliss::io::NFilterSequencesIterator<Iterator, SeqParser>;

// sequence iterator for use in identifying the first and last valid kmers in a read.  we want 1 sequence per read, so just use non-filtering and non-splitting kmer reader.
template <typename Iterator, template <typename> class SeqParser>
using SeqIterType = bliss::io::SequencesIterator<Iterator, SeqParser>;

using FileReaderType = ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, FileParser >;

using DBGMapType = ::bliss::debruijn::graph::simple_hash_debruijn_graph_map<KmerType>;
using DBGType = ::bliss::debruijn::graph::debruijn_graph<KmerType, bool>;

using CountDBGMapType = ::bliss::debruijn::graph::count_hash_debruijn_graph_map<KmerType, CountType>;
using CountDBGType = ::bliss::debruijn::graph::debruijn_graph<KmerType, CountType>;

using ChainNodeType = ::bliss::debruijn::simple_biedge<KmerType>;

//template <typename K>
//using ChainMapParams = typename DBGType::map_params_template<K>;
//
//using ChainGraphType = ::dsc::densehash_map<KmerType, ChainNodeType,
//		ChainMapParams,
//		 ::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

using ChainGraphType = ::bliss::debruijn::graph::debruijn_chain_graph<KmerType>;

template <typename Key>
using FreqMapParams = ::bliss::index::kmer::CanonicalHashMapParams<Key>;

using CountMapType = ::dsc::saturating_counting_densehash_map<KmerType, CountType,
		FreqMapParams,
		::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

using CountIndexType = ::bliss::index::kmer::CountIndex2<CountMapType>;


using FreqSummaryType = std::tuple<size_t, size_t, CountType, CountType>;
using FreqMapType = ::dsc::reduction_densehash_map<KmerType, FreqSummaryType,
		FreqMapParams,
		::bliss::kmer::hash::sparsehash::special_keys<KmerType, true>,
		::bliss::debruijn::operation::chain::freq_summary<CountType> >;

using ListRankedChainNodeVecType = std::vector<::bliss::debruijn::chain::listranked_chain_node<KmerType> >;

using ChainVecType = ::std::vector<std::pair<KmerType, ChainNodeType> >;



  /// append to container via emplace.
  /// modified based on http://stackoverflow.com/questions/18724999/why-no-emplacement-iterators-in-c11-or-c14
  template<class Container>
  class discard_iterator : public std::iterator< std::output_iterator_tag,
                                                      typename::std::iterator_traits<decltype(std::declval<Container>().begin())>::value_type,
                                                      void, void, void >
  {
  public:
      typedef Container container_type;

      discard_iterator() {}

      template<class... Args>
      discard_iterator& operator=(Args&&... args) { return *this; }
      discard_iterator& operator*() { return *this; }
      discard_iterator& operator++() { return *this; }
      discard_iterator& operator++(int) { return *this; }
  };



template <typename Index>
void build_index(::std::vector<::bliss::io::file_data> const & file_data, Index & idx, mxx::comm const & comm) {
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) printf("PARSING and INSERT\n");

	// TESTING
	::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > temp2;
	// TESTING END;

	for (auto x : file_data) {
		temp2.clear();

		BL_BENCH_START(build);
		::bliss::io::KmerFileHelper::template parse_file_data<
		 ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>,
		 FileParser, SplitSeqIterType>(x, temp2, comm);
		BL_BENCH_COLLECTIVE_END(build, "parse", temp2.size(), comm);

		BL_BENCH_START(build);
		idx.insert(temp2);
		BL_BENCH_COLLECTIVE_END(build, "insert", idx.local_size(), comm);
	}

	size_t total = idx.size();
	if (comm.rank() == 0) printf("PARSING and INSERT DONE: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "construct", comm);
}

#if defined(MIN_MEM)

template <typename Index>
void build_index_incremental(::std::vector<::bliss::io::file_data> const & file_data, Index & idx, mxx::comm const & comm) {
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) printf("PARSING and INSERT\n");

    using CharIterType = typename ::bliss::io::file_data::const_iterator;
	using SeqParserType = FileParser<CharIterType>;
	using SeqIterType = SplitSeqIterType<CharIterType, FileParser>;
	using KmerParser = ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>;
	using Iter = typename ::bliss::iterator::ContainerConcatenatingIterator<SeqIterType, KmerParser>;


	// TESTING
    BL_BENCH_START(build);


//	// estimate the largest amount of memory to use.
//	unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);
//	// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
//	size_t block_size = (free_mem / (8 * sizeof(typename KmerParser::value_type)));  // number of elements that can be held in freemem
//
//	if (comm.rank() == 0) std::cout << "estimate num elements=" << block_size << ", value_type size=" <<
//			sizeof(typename KmerParser::value_type) << " bytes" << std::endl;


	::std::vector<typename KmerParser::value_type> temp2;
	 // ::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType>

//	temp2.reserve(block_size);

    ::fsc::back_emplace_iterator<std::vector<typename KmerParser::value_type> > emplace_iter(temp2);
    // TESTING END;
    BL_BENCH_END(build, "reserve", temp2.capacity());

	BL_BENCH_LOOP_START(build, 0);  // for init
	BL_BENCH_LOOP_START(build, 1);  // for parse
	BL_BENCH_LOOP_START(build, 2);  // for insert
	size_t count = 0, i;
	bool all_done = false;

	for (auto x : file_data) {

		// initialization
		BL_BENCH_LOOP_RESUME(build, 0);

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

    	// estimate the largest amount of memory to use.
    	unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);

    	// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
    	size_t block_size = (free_mem / (8 * sizeof(typename KmerParser::value_type)));  // number of elements that can be held in freemem
    	block_size = std::min(block_size, x.getRange().size());

    	if (comm.rank() == 0) std::cout << "estimate num elements=" << block_size << ", value_type size=" <<
    			sizeof(typename KmerParser::value_type) << " bytes" << std::endl;

    	if (block_size > temp2.capacity()) temp2.clear();
    	temp2.reserve(block_size);


		BL_BENCH_LOOP_PAUSE(build, 0);

    	//=== copy into array incrementally
		while (! all_done) {
			temp2.clear();

            //== process the chunk of data
			BL_BENCH_LOOP_RESUME(build, 1);

			for (i = 0; (i < block_size) && (start != endd); ++i) {
				*emplace_iter = *start;
				++start;
				++emplace_iter;
			}
			count += i;

        	BL_BENCH_LOOP_PAUSE(build, 1);


			BL_BENCH_LOOP_RESUME(build, 2);
			all_done = (i == 0);
			all_done = mxx::all_of(all_done, comm);

			idx.insert(temp2);
			BL_BENCH_LOOP_PAUSE(build, 2);
		}
	}
	BL_BENCH_LOOP_END(build, 0, "setup", temp2.capacity());  // for init
	BL_BENCH_LOOP_END(build, 1, "parse", count);  // for parse
	BL_BENCH_LOOP_END(build, 2, "insert", idx.local_size());  // for insert

	size_t total = idx.size();
	if (comm.rank() == 0) printf("PARSING and INSERT DONE: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "construct", comm);
}
#endif


// ==========  choices:
//  1. no filtering:  parse simple nodes directly, insert
//  2. de novo filtering:  parse simple nodes, parse edge k1mers, transform edges, copy edges, filter and erase, insert
//  3. reconstructive filtering:  parser kmers, add edges, filter and copy, insert.

// FILTERING ONLY WORKS WITH FASTQ files right now.  TODO: STILL NEEDED?
//#if (pPARSER == FASTQ)

/*
 * @brief  build an index with thresholded k+2-mers.  note that the threshold is specified for k+2-mers, not k-mers, and is exclusive.
 * @details		The goal is to identify erroneous edges and nodes.  k+2 mer satisfies this goal.  the only case that is
 * 				not dealt with is when there is 1 erroneous edge, and the other edge is good,
 * 				in which case the good edge's final count would be lower.
 * 				The hope is that the next k+2-mer would have enough count to make this boundary case insignificant.
 * @note		The k-mer counts are generated from the center k of k+2-mers
 * @note 		Filtering after the dbg is build loses context of biedge - can only operate on edges, and the central k-mer
 * 				(dbg node) count is based on the biedge.
 * @note		Finally, counting k+1-mer is tricky also because it does not have the symmetry, so canonical changes which side
 * 				of k+1-mer the owning k-mer sits on.  This means that easiest way to do the k-mer
 * 				counting is to count both ends.  however, due to high number of reads, not all k-mers are counted 2x,
 * 				 so we can't simply divide the count by 2 to get the true count.
 * @note		counting using k+2 mer and filter would result in missing nodes - prev node with valid edge pointing to (non-existent) node with invalid edge
 * 				have to use a hybrid approach - read in k+2-mers, but check if should insert using 2 k+1 mer.
 * @tparam Index		Type of the debruijn graph
 * @param file_data		Input raw file type
 * @param idx			debruijn graph to build
 * @param lower_thresh	inclusive lower threshold for INCLUDING a k+1 mers.  note this is not threshold of edge or node frequency.
 * @param upper_thresh	exclusive upper threshold for INCLUDING a k+1 mers.  note this is not threshold of edge or node frequency.
 * @param comm			mpi communicator.
 * @return      vector of vectors of bool indicating whether the left and right k+1mer of the k+2-mer nodes are within valid frequency ranges or not.
 */
template <typename Index>
::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> >
build_index_thresholded(::std::vector<::bliss::io::file_data> const & file_data, Index & idx,
		size_t const & lower_thresh, size_t const & upper_thresh,  mxx::comm const & comm, bool benchmark = true) {
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) printf("PARSING, FILTER, and INSERT\n");

	// k+1-mer count map.  note that it should use the same hash function as Index.
	using K1merType	= ::bliss::common::Kmer<KmerType::size + 1, typename KmerType::KmerAlphabet, typename KmerType::KmerWordType>;
	using CountMap1Type = ::dsc::saturating_counting_densehash_map<K1merType, CountType,
	    FreqMapParams,
			::bliss::kmer::hash::sparsehash::special_keys<K1merType, true> >;


	// ========  count the k+1-mers.
	BL_BENCH_START(build);
  CountMap1Type counter2(comm);
  ::bliss::debruijn::biedge::filter::compute_edge_frequency<FileParser, SplitSeqIterType>(file_data, counter2, comm, lower_thresh, upper_thresh);
	BL_BENCH_COLLECTIVE_END(build, "count_and_filter k1mer", counter2.local_size(), comm);


	// ======= filter k+2-mers.  incremental by file
	::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> > results;
	{
		::std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> > nodes2;

		for (auto x : file_data) {

			BL_BENCH_START(build);
			nodes2.clear();
			// construct biedges (nodes)
			::bliss::io::KmerFileHelper::template parse_file_data<
			  ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>,
			   FileParser, SplitSeqIterType>(x, nodes2, comm);
			BL_BENCH_COLLECTIVE_END(build, "parse", nodes2.size(), comm);

			// remove low frequency edges
			BL_BENCH_START(build);
			::bliss::debruijn::biedge::filter::transform_biedges_by_frequency<KmerType, CountMap1Type>(nodes2, counter2, comm);

			// extract the biedges. and save them
			if (!benchmark) results.emplace_back(::bliss::debruijn::biedge::filter::extract_biedges<KmerType>(nodes2));
			BL_BENCH_END(build, "transform_by_freq", nodes2.size());

			// now remove nodes with no edges.
			BL_BENCH_START(build);
			::bliss::debruijn::biedge::filter::filter_isolated_nodes<KmerType>(nodes2);
			BL_BENCH_END(build, "filter", nodes2.size());

			// build the index
			BL_BENCH_START(build);
			idx.insert(nodes2);
			BL_BENCH_COLLECTIVE_END(build, "insert", idx.local_size(), comm);
//			std::cout << "rank " << comm.rank() << " input size " << node_size << " idx size " << idx.local_size() << " buckets " << idx.get_map().get_local_container().bucket_count() << std::endl;

		}
	}

	size_t total = idx.size();
	if (comm.rank() == 0) printf("PARSING, FILTER, and INSERT: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "filtered_construct", comm);

	return results;
}

#if defined(MIN_MEM)

template <typename Index>
::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> >
build_index_thresholded_incremental(::std::vector<::bliss::io::file_data> const & file_data, Index & idx,
		size_t const & lower_thresh, size_t const & upper_thresh,  mxx::comm const & comm, bool benchmark = true) {
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) printf("PARSING, FILTER, and INSERT\n");


	// ========  count the k+1-mers.
	using K1merType	= ::bliss::common::Kmer<KmerType::size + 1, typename KmerType::KmerAlphabet, typename KmerType::KmerWordType>;
	// k+1-mer count map.  note that it should use the same hash function as Index.
	using CountMap1Type = ::dsc::saturating_counting_densehash_map<K1merType, CountType,
	    FreqMapParams,
			::bliss::kmer::hash::sparsehash::special_keys<K1merType, true> >;

	BL_BENCH_START(build);
	CountMap1Type counter2(comm);
	::bliss::debruijn::biedge::filter::compute_edge_frequency_incremental<FileParser, SplitSeqIterType>(
			file_data, counter2, comm, lower_thresh, upper_thresh);
	BL_BENCH_END(build, "count_and_filter k1mer", counter2.local_size());


	// ======== filter by edge frequency and insert

  using CharIterType = typename ::bliss::io::file_data::const_iterator;
	using SeqParserType = FileParser<CharIterType>;
	using SeqIterType = SplitSeqIterType<CharIterType, FileParser>;
	using KmerParser = ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>;
	using Iter = typename ::bliss::iterator::ContainerConcatenatingIterator<SeqIterType, KmerParser>;

	using EdgeType = ::bliss::debruijn::biedge::compact_simple_biedge;
	//using NodeType = std::pair<KmerType, EdgeType>;

	::std::vector< ::std::vector<EdgeType> > results;  // each inner vector corresponds to one of the data objects.


	size_t total0 = 0;
	size_t total1 = 0;
	size_t total2 = 0;
	size_t total3 = 0;

	BL_BENCH_LOOP_START(build, 0);  // for init
	BL_BENCH_LOOP_START(build, 1);  // for reserve
	BL_BENCH_LOOP_START(build, 2);  // for filter_insert
	BL_BENCH_LOOP_START(build, 3);  // for save_edges

	for (auto x : file_data) {

		// initialization
		BL_BENCH_LOOP_RESUME(build, 0);

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

    	total0 += x.getRange().size();
    BL_BENCH_LOOP_PAUSE(build, 0);


        //== process the chunk of data, and filter insert into index.
		BL_BENCH_LOOP_RESUME(build, 1);
		// initialize a edge vector.
		::std::vector<EdgeType> edges;
		edges.reserve(x.getRange().size() / 4);  // estimate.  may be small for fasta, and large for fastq.
		total1 += edges.capacity();
		BL_BENCH_LOOP_PAUSE(build, 1);

		BL_BENCH_LOOP_RESUME(build, 2);
	  // estimate the largest amount of memory to use.
	  unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);
	  // use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
	  size_t block_size = (free_mem / (8 * sizeof(typename CountMap1Type::key_type)));  // number of elements that can be held in freemem
	  block_size = std::min(block_size, x.getRange().size());

	  if (comm.rank() == 0) std::cout << "estimate num elements=" << block_size << ", value_type size=" <<
	      sizeof(typename CountMap1Type::key_type) << " bytes" << std::endl;



		// filter and insert
		if (!benchmark) 
			total2 += ::bliss::debruijn::biedge::filter::freq_filter_insert_biedges<
				KmerType, Iter, CountMap1Type, Index,
				::fsc::back_emplace_iterator<std::vector<EdgeType> > >(start, endd, counter2, block_size, idx, ::fsc::back_emplace_iterator<std::vector<EdgeType> >(edges), comm);
		else 
			total2 += ::bliss::debruijn::biedge::filter::freq_filter_insert_biedges<
				KmerType, Iter, CountMap1Type, Index,
				discard_iterator<std::vector<EdgeType> > >(start, endd, counter2, block_size, idx, discard_iterator<std::vector<EdgeType> >(), comm);
		BL_BENCH_LOOP_PAUSE(build, 2);

		BL_BENCH_LOOP_RESUME(build, 3);
		total3 += edges.size();
		// save the results.
		if (!benchmark) results.emplace_back(std::move(edges));
		BL_BENCH_LOOP_PAUSE(build, 3);

	}
	BL_BENCH_LOOP_END(build, 0, "setup", 		 total0);  // for init. size of input data
	BL_BENCH_LOOP_END(build, 1, "reserve", 		 total1);  // for parse.  reserved capacity for edges.
	BL_BENCH_LOOP_END(build, 2, "filter_insert", total2);  // total attempted insertion
	BL_BENCH_LOOP_END(build, 3, "save_edges", 	 total3);  // total edges == kmers.

	// for information only.
	BL_BENCH_START(build);
	BL_BENCH_END(build, "index", idx.local_size());


	size_t total = idx.size();
	if (comm.rank() == 0) printf("PARSING, FILTER, and INSERT: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "filtered_construct", comm);

	return results;
}


#endif

#include "../common/compact_dbg_io.hpp"

#include "../common/compact_dbg_stats.hpp"


/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {

	//////////////// init logging
	LOG_INIT();

	//////////////// initialize MPI and openMP

	mxx::env e(argc, argv);
	mxx::comm comm;

	if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);

	comm.barrier();


	//////////////// parse parameters

	//////////////// parse parameters
	std::vector<std::string> filenames;
	std::string filename;


	std::string out_prefix;
	out_prefix.assign("./output");

//#if (pPARSER == FASTQ)
	size_t lower, upper;
	lower = upper = ((sizeof(CountType) > 4) ? ::std::numeric_limits<uint32_t>::max() : ::std::numeric_limits<CountType>::max()) + 1;
//#endif

	bool thresholding = false;
	bool benchmark = false;
	bool LRoptimized = false;
	bool compress = false;
  bool mpiio = false;

	//  std::string queryname(filename);
	//  int sample_ratio = 100;

	// Wrap everything in a try block.  Do this every time,
	// because exceptions will be thrown for problems.
	try {

		// Define the command line object, and insert a message
		// that describes the program. The "Command description message"
		// is printed last in the help text. The second argument is the
		// delimiter (usually space) and the last one is the version number.
		// The CmdLine object parses the argv array based on the Arg objects
		// that it contains.
		TCLAP::CmdLine cmd("Parallel de bruijn graph compaction", ' ', "0.4");

		// MPI friendly commandline output.
		::bliss::utils::tclap::MPIOutput cmd_output(comm);
		cmd.setOutput(&cmd_output);

		// Define a value argument and add it to the command line.
		// A value arg defines a flag and a type of value that it expects,
		// such as "-n Bishop".
		TCLAP::ValueArg<std::string> outputArg("O", "output_prefix", "Prefix for output files, including directory", false, "", "string", cmd);

		// value type is uint16_t - to ensure that the values are parsed as integer instead of being treated as unsigned chars.
		TCLAP::SwitchArg threshArg("T", "thresholding", "on/off for thresholding", cmd, false);
		TCLAP::ValueArg<size_t> lowerThreshArg("L", "lower_thresh", "Lower Threshold for Kmer and Edge frequency", false, 0, "uint16", cmd);
		TCLAP::ValueArg<size_t> upperThreshArg("U", "upper_thresh", "Upper Threshold for Kmer and Edge frequency", false,
				upper, "uint16", cmd);

		TCLAP::SwitchArg benchmarkArg("B", "benchmark", "on/off for benchmarking (no file output)", cmd, false);

		TCLAP::SwitchArg lrOptimizeArg("R", "list_rank_opt", "on/off for list ranking optimization", cmd, false);

    TCLAP::SwitchArg compressArg("C", "compress", "on/off for in-mem compressed chains", cmd, false);

    TCLAP::SwitchArg mpiioArg("M", "mpiio", "on - use mpiio for input.  off - use posix io", cmd, false);


		//    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);
		TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names", false, "string", cmd);


		// Parse the argv array.
		cmd.parse( argc, argv );

		filenames = fileArg.getValue();
		out_prefix = outputArg.getValue();

//#if (pPARSER == FASTQ)
    lower = ::std::min(lowerThreshArg.getValue(), lower);
    upper = ::std::min(upperThreshArg.getValue(), upper);
//#endif

		thresholding = threshArg.getValue();
		benchmark = benchmarkArg.getValue();
		LRoptimized = lrOptimizeArg.getValue();
		compress = compressArg.getValue();
    mpiio = mpiioArg.getValue();


	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		exit(-1);
	}

	if (comm.rank() == 0)  std::cout << "compact_debruijn_graph_block_construct v0.4a" << std::endl;
  if (thresholding && (comm.rank() == 0)) {
    std::cout << "THRESHOLDS: ";
    std::cout << "k1:  " << lower << "-" << upper << std::endl;
  }


//#if (pPARSER == FASTA)
//    if (thresholding)
//      throw std::invalid_argument("ERROR: FASTA version of debruijn graph compaction does not currently support filtering.");
//#endif

	if (filenames.size() == 0) {
		filename.assign(PROJ_SRC_DIR);

#if (pPARSER == FASTA)
		filename.append("/test/data/test.debruijn.small.fasta");
#elif (pPARSER == FASTQ)
		filename.append("/test/data/test.debruijn.small.fastq");
#endif


		filenames.push_back(filename);
	}

	// filename for compacted chain strings
	// string starts with smaller end.  (first K and rev_comp of last K compared)
	// this is a dump of the collected compacted chains.
	std::string compacted_chain_str_filename(out_prefix);
	compacted_chain_str_filename.append("_chain.fasta");

	// compressed chain
	std::string compressed_chain_filename(out_prefix);
	compressed_chain_filename.append("_compressed_chain.debug");

	// filename for compacted chain interior kmers.  in format <K, Chain Id, pos, +/->
	// K is canonical.  + if K is on same strand as chain, - if not.
	// this is a dump of the chain map nodes.
	// chain_id may not be canonical, but is the smaller of the 5' ends of a bimolecule.
	// NOTE: to read 5' to 3', we'd never start a chain with k-mer that has in edge but not out edge,
	//   so chain id is the k-mer at 5' end.  however, since all k-mer indexing are done canonically,
	//   canonical k-mer is +.  so for a non-canonical chain id, the terminal K is at - strand.
	//   all other k-mers on this chain have strand value relative to canonical of chain id.
	std::string compacted_chain_kmers_filename(out_prefix);
	compacted_chain_kmers_filename.append("_chain.components");

	// filename for compacted chain end points, in format : < K_l, K_r, chain_id, f_l_A, f_l_C, f_l_G, f_l_T, f_r_A, f_r_C, f_r_G, f_r_T, f >
	// K_l is canonical.  as such, we need chain_id to link back to component k-mers and the compated strings.
	// K_l is the smaller of K_l and rev_comp(K_r).  (compare 5')
	// K_l and K_r are on same strand.
	// f is frequency - minimum or average of all intervening nodes?  - a reduction similar to the one used to construct the compacted chain string is needed.
	// f_l_X and f_r_X are in and out edges of K_l and K_r, respectively.
	std::string compacted_chain_ends_filename(out_prefix);
	compacted_chain_ends_filename.append("_chain.edges");

	// filename for junctions in format <K, f_l_A, f_l_C, f_l_G, f_l_T, f_r_A, f_r_C, f_r_G, f_r_T, f>
	// K is canonical.
	// f is frequency
	// f_l_X and f_r_X are in and out edges of K, respectively.
	// this is a dump of the dbg junctional nodes (filtered) to disk.
	std::string branch_filename(out_prefix);
	branch_filename.append("_branch.edges");

  std::string branch_fasta_filename(out_prefix);
  branch_fasta_filename.append("_branch.fasta");

	// std::string graph_filename(out_prefix);
	// graph_filename.append("_graph.edges.debug");
	std::string chain_biedge_filename(out_prefix);
	chain_biedge_filename.append("_chainmap.debug");
	// std::string graph_filename(out_prefix);
	// graph_filename.append("_graph.edges.debug");
	std::string chain_biedge_filename2(out_prefix);
	chain_biedge_filename2.append("_chainmap_uncomp.debug");


  BL_BENCH_INIT(app);

  // ================  read and get file
  BL_BENCH_START(app);
  ::std::vector<::bliss::io::file_data> file_data = open_files(filenames, comm, mpiio);
  BL_BENCH_COLLECTIVE_END(app, "read_file", file_data.size(), comm);
  // == DONE == reading

	::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> > selected_edges;
	ChainGraphType chainmap(comm);
  DBGType idx(comm);

	// =================  make compacted simple DBG, so that we can get chain and branch kmers.
	{
		BL_BENCH_START(app);

		// first read the input



#if defined(MIN_MEM)
		if (thresholding) {
//#if (pPARSER == FASTQ)
			build_index_thresholded_incremental(file_data, idx, lower, upper, comm, benchmark).swap(selected_edges);
//#else
//			// TODO: VERIFY FASTA is working. THIS IS WORKING.  FASTA does not support thresholded index build because fasta reader is not yet handling overlaps correctly when dealing with k+2-mers.
//			if (comm.rank() == 0) printf("ERROR: FASTA files does not yet support pre-filter by frequency.\n");
//#endif			

		} else {
			build_index_incremental(file_data, idx, comm);
		}
#else
		if (thresholding) {
			build_index_thresholded(file_data, idx, lower, upper, comm, benchmark).swap(selected_edges);
		} else {
			build_index(file_data, idx, comm);
		}
#endif

		if (benchmark) {
			file_data.clear();
		}


		BL_BENCH_COLLECTIVE_END(app, "construct", idx.local_size(), comm);


		if (!benchmark) {
			BL_BENCH_START(app);
			print_edge_histogram(idx, comm);
      		check_index(idx, comm);
      		if (comm.rank() == 0) printf("rank 0 finished checking index\n");
			BL_BENCH_COLLECTIVE_END(app, "histo", idx.local_size(), comm);
			// == DONE == make compacted simple DBG

		}
		// TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet
		//   this is done via read filtering/splitting.

//   BL_BENCH_START(app);
//   print_graph_edge_frequencies(graph_filename, idx, comm);
//   BL_BENCH_COLLECTIVE_END(app, "print graph", idx.local_size(), comm);


		// == PRINT == prep branch for printing - here ONLY BECAUSE WE ARE DISCARDING IDX AFTER MAKING CHAINMAPS
		if (!benchmark) {
			// TODO: try using CountDBGType as base instead of DBGType, then this data is captured to begin with.
			CountDBGType idx2(comm);
			{
				BL_BENCH_START(app);
				std::vector<KmerType> kmers = idx.get_branch_node_kmers();
//				get_branch_kmers(idx, kmers, comm);
				BL_BENCH_COLLECTIVE_END(app, "branch_kmers", kmers.size(), comm);

				BL_BENCH_START(app);
				// same distribution (using same hashmap params) as idx, so is_local can be set to true.
				idx2.insert(kmers, true);
				assert(idx2.local_size() == kmers.size());  // should be 1 to 1.
				BL_BENCH_COLLECTIVE_END(app, "make_branch_counter", kmers.size(), comm);
			}

			{
				// get the edges counts for these kmers.
				BL_BENCH_START(app);
				count_edges(file_data, selected_edges, thresholding, idx2, comm);
				BL_BENCH_COLLECTIVE_END(app, "branch_freq", idx2.local_size(), comm);
			}  // enforce delete kmers vec

			// == print branches.
			BL_BENCH_START(app);
			print_branch_edge_frequencies(branch_filename, idx2, comm);
			BL_BENCH_COLLECTIVE_END(app, "print branch", idx2.local_size(), comm);
		
			BL_BENCH_START(app);
			print_branch_fasta(branch_fasta_filename, idx2, comm);
			BL_BENCH_COLLECTIVE_END(app, "print branch fasta", idx2.local_size(), comm);
		}  // enforce delete idx2.

		// ==== make chain map
		BL_BENCH_START(app);
		chainmap.extract_chains(idx);
		//make_chain_map(idx, chainmap, comm);
		if (benchmark) {
			idx.get_map().clear();
			idx.get_map().reserve(0);
		}

		BL_BENCH_COLLECTIVE_END(app, "chainmap", chainmap.local_size(), comm);
		// == DONE == make chain map
	} // enforce delete idx.

#ifndef NDEBUG
		BL_BENCH_START(app);
		print_chain_biedges(chain_biedge_filename2, chainmap, comm);
		BL_BENCH_COLLECTIVE_END(app, "print_chain_biedge_uncomp", chainmap.local_size(), comm);
#endif

	// ===== parallel list ranking for chain compaction
	BL_BENCH_START(app);
	size_t iterations = 0;
	if (LRoptimized)
		iterations = chainmap.list_rank_min_update();
	else
		iterations = chainmap.list_rank();
	//auto cycle_node_kmers = list_rank(chainmap, comm);
	BL_BENCH_COLLECTIVE_END(app, "list_rank", iterations, comm);
	// == DONE == parallel list ranking for chain compaction

	if (benchmark) {

		// now print chain string - order is destroyed via psort.
		if (!compress) {
			// prepare
			BL_BENCH_START(app);
			ListRankedChainNodeVecType compacted_chain = chainmap.to_ranked_chain_nodes();
			BL_BENCH_COLLECTIVE_END(app, "compact_chain", compacted_chain.size(), comm);

			if (benchmark) {
				chainmap.clear();
			}

			BL_BENCH_START(app);
			print_chain_string(compacted_chain_str_filename, compacted_chain, comm);
			BL_BENCH_COLLECTIVE_END(app, "chain_str", compacted_chain.size(), comm);
		} else {
			// === and compress.
			BL_BENCH_START(app);
			::std::vector<::std::string> compressed_chain = chainmap.to_compressed_chains();
			BL_BENCH_COLLECTIVE_END(app, "compress_chains", compressed_chain.size(), comm);

			if (benchmark) {
				chainmap.clear();
			}

			BL_BENCH_START(app);
			size_t out_size = print_compressed_chains(compressed_chain_filename, compressed_chain, comm);
			BL_BENCH_COLLECTIVE_END(app, "chain_str", out_size, comm);
		}
	}

	// =============================================================
	// below is for printing.


	if (!benchmark) {
#ifndef NDEBUG
		BL_BENCH_START(app);
		print_chain_biedges(chain_biedge_filename, chainmap, comm);
		BL_BENCH_COLLECTIVE_END(app, "print_chain_biedge", chainmap.local_size(), comm);
#endif
		// =========== remove cycles and isolated
		BL_BENCH_START(app);
		auto cycle_kmers = chainmap.get_cycle_node_kmers();
		idx.erase(cycle_kmers);
		idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
		BL_BENCH_COLLECTIVE_END(app, "remove cycles/isolated/etc", idx.local_size(), comm);


		//  auto lidx = idx.get_map().get_local_container();
		//  for (auto it = lidx.begin(); it != lidx.end(); ++it) {
		//    std::cout << "valid kmer in index: " << (*it).first << std::endl;
		//  }

		// == PRINT == valid k-mers files
		if (thresholding) {

#if (pPARSER == FASTA)
			if (comm.rank() == 0) printf("WARNING: outputting first/last valid kmer position for each read is supported for FASTQ format only.\n");
#elif (pPARSER == FASTQ)

			BL_BENCH_START(app);
			for (size_t i = 0; i < filenames.size(); ++i) {
				std::string fn(out_prefix);

				fn.append(".");
				fn.append(std::to_string(i));
				fn.append(".valid");

				print_valid_kmer_pos_in_reads(fn, file_data[i], idx, comm);
			}
			BL_BENCH_COLLECTIVE_END(app, "print_valid_kmer_pos", filenames.size(), comm);
#endif
		}



		// == print ===  prepare for printing compacted chain and frequencies

		// search in chainmap to find canonical termini.
		BL_BENCH_START(app);
		ChainVecType chain_rep = chainmap.find_if(::bliss::debruijn::filter::chain::IsCanonicalTerminusOrIsolated());
		BL_BENCH_COLLECTIVE_END(app, "chain rep", chain_rep.size(), comm);

		BL_BENCH_START(app);
		//ChainVecType termini = chainmap.find_if(::bliss::debruijn::filter::chain::IsTerminusOrIsolated());
		std::vector<KmerType> termini = chainmap.get_terminal_node_kmers();
		BL_BENCH_COLLECTIVE_END(app, "chain termini", termini.size(), comm);


		FreqMapType freq_map(comm);
		{
			// prepare
			BL_BENCH_START(app);
			ListRankedChainNodeVecType compacted_chain = chainmap.to_ranked_chain_nodes();
			BL_BENCH_COLLECTIVE_END(app, "compacted_chain", compacted_chain.size(), comm);


			if (compress) {
				// === and compress.
				BL_BENCH_START(app);
				::std::vector<::std::string> compressed_chain = chainmap.to_compressed_chains();
				BL_BENCH_COLLECTIVE_END(app, "compress_chains", compressed_chain.size(), comm);

				BL_BENCH_START(app);
				size_t out_size = print_compressed_chains(compressed_chain_filename, compressed_chain, comm);
				BL_BENCH_COLLECTIVE_END(app, "chain_str", out_size, comm);
			}

			// do this first because we need the order of compacted chain to be same as hashed distribution.
			{
				// compute count index
				BL_BENCH_START(app);
				CountIndexType count_idx(comm);
				count_kmers(file_data, selected_edges, thresholding, count_idx, comm);
				BL_BENCH_COLLECTIVE_END(app, "count_kmers", count_idx.local_size(), comm);

				// compute freq map
				BL_BENCH_START(app);
#if defined(MIN_MEM)
				compute_freq_map_incremental(compacted_chain, count_idx, freq_map, comm);
#else
				compute_freq_map(compacted_chain, count_idx, freq_map, comm);
#endif			
				BL_BENCH_COLLECTIVE_END(app, "chain_freqs", freq_map.local_size(), comm);
			} // ensure delet count_index



			// now print chain string - order is destroyed via psort.
			BL_BENCH_START(app);
			print_chain_string(compacted_chain_str_filename, compacted_chain, comm);
			BL_BENCH_COLLECTIVE_END(app, "chain_str", compacted_chain.size(), comm);

			BL_BENCH_START(app);
			print_chain_nodes(compacted_chain_kmers_filename, compacted_chain, comm);
			BL_BENCH_COLLECTIVE_END(app, "chain_node", compacted_chain.size(), comm);

		} // ensure release compacted chain



		// get terminal k-mers
		CountDBGType idx2(comm);
		{

			BL_BENCH_START(app);
			// same distribution (using same hashmap params) as idx, so is_local can be set to true.
			idx2.insert(termini, true);
			assert(idx2.local_size() == termini.size());  // should be 1 to 1.
			BL_BENCH_COLLECTIVE_END(app, "make_terminal_counter", termini.size(), comm);


			// get the edges counts for these kmers.
			BL_BENCH_START(app);
			count_edges(file_data, selected_edges, thresholding, idx2, comm);
			BL_BENCH_COLLECTIVE_END(app, "terminal_edge_freq", idx2.local_size(), comm);
		} // ensure delete kmers.


		BL_BENCH_START(app);
		print_chain_frequencies(compacted_chain_ends_filename, chain_rep, idx2, freq_map, comm);
		BL_BENCH_COLLECTIVE_END(app, "print_chain_freq", chain_rep.size(), comm);

					// release chainmap
			BL_BENCH_START(app);
			chainmap.get_map().get_local_container().reset();
			BL_BENCH_COLLECTIVE_END(app, "chainmap_reset", chainmap.local_size(), comm);

	}

	BL_BENCH_REPORT_MPI_NAMED(app, "app", comm);


	// mpi cleanup is automatic

	return 0;

}

