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
 * compact_dbg_countIndex_construct.cpp
 *
 * goal: 	this version constructs a k+2-mer count index first, then uses that to filter (construct k+1-mer counter locally)
 *
 *  Created on: Aug 17, 2015
 *      Author: yongchao
 *
 *  Rewrote: June 10, 2016
 *      Author: tony pan
 *
 */

#include "bliss-config.hpp"

#include <unistd.h> // get hostname

#include <functional>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <chrono>
#include <iostream> // for system("pause");
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
#include "debruijn/debruijn_topo_operations.hpp"

#include "debruijn/debruijn_stats.hpp"
#include "utils/minimizer_hash.hpp"

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
//using CountType = uint16_t;  // set to uint16_t by default.  this is to ensure that most frequencies are captured.

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

		#if (((bDNA * (pK + 1)) & 0x3F) == 0)
			using K1merWordType = uint64_t;
		#elif (((bDNA * (pK + 1)) & 0x1F) == 0)
			using K1merWordType = uint32_t;
		#elif (((bDNA * (pK + 1)) & 0x0F) == 0)
			using K1merWordType = uint16_t;
		#else
			using K1merWordType = uint8_t;
		#endif	
		using K1merType = bliss::common::Kmer<(pK + 1), Alphabet, K1merWordType>;

		#if (((bDNA * (pK + 2)) & 0x3F) == 0)
			using K2merWordType = uint64_t;
		#elif (((bDNA * (pK + 2)) & 0x1F) == 0)
			using K2merWordType = uint32_t;
		#elif (((bDNA * (pK + 2)) & 0x0F) == 0)
			using K2merWordType = uint16_t;
		#else
			using K2merWordType = uint8_t;
		#endif	
		using K2merType = bliss::common::Kmer<(pK + 2), Alphabet, K2merWordType>;

	#else // pK not defined.  assume 31.  compute data type.
		#if (pDNA == 16)
			using KmerWordType = uint64_t;
			using K1merWordType = uint64_t;
			using K2merWordType = uint8_t;
		#elif (pDNA == 5)
			using KmerWordType = uint32_t;
			using K1merWordType = uint32_t;
			using K2merWordType = uint8_t;
		#elif (pDNA == 4)
			using KmerWordType = uint64_t;
			using K1merWordType = uint64_t;
			using K2merWordType = uint8_t;
		#endif
		using KmerType = bliss::common::Kmer<31, Alphabet, KmerWordType>;
		using K1merType = bliss::common::Kmer<32, Alphabet, K1merWordType>;
		using K2merType = bliss::common::Kmer<33, Alphabet, K2merWordType>;
	#endif

		using CountType = uint8_t; // set to uint16_t by default.  this is to ensure that most frequencies are captured.
#else					   // MIN_MEM not defined.  use uint64_t.  this is same as compact_debruijn_graph_refactor.cpp

	using KmerWordType = uint64_t;  // matches system architecture.
	using K1merWordType = uint64_t; // matches system architecture.
	using K2merWordType = uint64_t; // matches system architecture.

	#if defined(pK)
		using KmerType = bliss::common::Kmer<pK, Alphabet, KmerWordType>;
		using K1merType = bliss::common::Kmer<(pK + 1), Alphabet, K1merWordType>;
		using K2merType = bliss::common::Kmer<(pK + 2), Alphabet, K2merWordType>;
	#else
		using KmerType = bliss::common::Kmer<31, Alphabet, KmerWordType>;
		using K1merType = bliss::common::Kmer<32, Alphabet, K1merWordType>;
		using K2merType = bliss::common::Kmer<33, Alphabet, K2merWordType>;
	#endif
		using CountType = uint16_t; // set to uint16_t by default.  this is to ensure that most frequencies are captured.
#endif

#define FASTA 1
#define FASTQ 0

//============== index input file format
#if (pPARSER == FASTA)
#define FileParser ::bliss::io::FASTAParser
#elif (pPARSER == FASTQ)
#define FileParser ::bliss::io::FASTQParser
#endif

#if defined(MINIMIZER)
template <typename KmerType>
using KmerDistHash = ::bliss::kmer::hash::minimizer::murmur<KmerType, true>;
#else
template <typename KmerType>
using KmerDistHash = ::bliss::kmer::hash::murmur<KmerType, true>;
#endif

using EdgeEncoding = Alphabet;

// sequence iterator for use in constructing the dbg.  here we throw away reads containing N.  technically, the splitting sequence iterator works here too,
// but the frequency calculated would include kmers in reads containing N
template <typename Iterator, template <typename> class SeqParser>
using SplitSeqIterType = bliss::io::NFilterSequencesIterator<Iterator, SeqParser>;

// sequence iterator for use in identifying the first and last valid kmers in a read.  we want 1 sequence per read, so just use non-filtering and non-splitting kmer reader.
template <typename Iterator, template <typename> class SeqParser>
using SeqIterType = bliss::io::SequencesIterator<Iterator, SeqParser>;

using FileReaderType = ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, FileParser>;

// using DBGMapType = ::bliss::debruijn::graph::simple_hash_debruijn_graph_map<KmerType>;
// used by benchmark?
using DBGType = ::bliss::debruijn::graph::debruijn_graph<KmerType, bool, KmerDistHash>;

//using CountDBGMapType = ::bliss::debruijn::graph::count_hash_debruijn_graph_map<KmerType, CountType>;
using CountDBGType = ::bliss::debruijn::graph::debruijn_graph<KmerType, CountType, KmerDistHash>;

using ChainNodeType = ::bliss::debruijn::simple_biedge<KmerType>;

//template <typename K>
//using ChainMapParams = typename DBGType::map_params_template<K>;
//
//using ChainGraphType = ::dsc::densehash_map<KmerType, ChainNodeType,
//		ChainMapParams,
//		 ::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

using ChainGraphType = ::bliss::debruijn::graph::debruijn_chain_graph<KmerType, KmerDistHash>;

template <typename Key>
using FreqMapParams = ::bliss::index::kmer::CanonicalHashMapParams<Key, KmerDistHash>;

// here for compilation purpose.  there are functions in compact_dbg_stats that need these.  this class does not.
 using CountMapType = ::dsc::counting_densehash_map<KmerType, CountType,
 		FreqMapParams,
 		::bliss::kmer::hash::sparsehash::special_keys<KmerType, true>>;
 using CountIndexType = ::bliss::index::kmer::CountIndex2<CountMapType>;

using FreqSummaryType = std::tuple<size_t, size_t, CountType, CountType>;
using FreqMapType = ::dsc::reduction_densehash_map<KmerType, FreqSummaryType,
		FreqMapParams,
		::bliss::kmer::hash::sparsehash::special_keys<KmerType, true>,
		::bliss::debruijn::operation::chain::freq_summary<CountType>>;

using ListRankedChainNodeVecType = std::vector<::bliss::debruijn::chain::listranked_chain_node<KmerType>>;

using ChainVecType = ::std::vector<std::pair<KmerType, ChainNodeType>>;

#include "../common/compact_dbg_build.hpp"

#include "../common/compact_dbg_io.hpp"

#include "../common/compact_dbg_stats.hpp"

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
size_t
build_index_thresholded(::std::vector<::bliss::io::file_data> const &file_data, Index &idx,
		std::vector<size_t> const &threshes, mxx::comm const &comm,
		std::string k2mer_filename)
		{
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) 
		printf("PARSING, FILTER, and INSERT\n");

	// need to build the k2mer counter first using all files
	BL_BENCH_START(build);
  	typename Index::map_type::LocalK2CountMapType k2_counter;
	BL_BENCH_COLLECTIVE_END(build, "init_k2counter", k2_counter.size(), comm);

	// ======= count k+2-mers.  incremental by file
	{
		::std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge>> nodes2;

		for (auto x : file_data)
		{

			BL_BENCH_START(build);
			nodes2.clear();
			// construct biedges (nodes)
			::bliss::io::KmerFileHelper::template parse_file_data<
			  ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>,
			   FileParser, SplitSeqIterType>(x, nodes2, comm);
			BL_BENCH_COLLECTIVE_END(build, "parse", nodes2.size(), comm);

			// compute the frequencies.
			BL_BENCH_START(build);
			idx.get_map().compute_biedge_freqencies(nodes2, k2_counter);
			BL_BENCH_END(build, "compute_freq", nodes2.size());
		}
	}
#ifndef NDEBUG
	print_k2mer_frequencies(k2mer_filename, k2_counter, comm);
#endif
	// then filter the k2mers and insert into dbg (no need to touch files again)
	BL_BENCH_START(build);
	std::vector<size_t> lthreshes(threshes.begin(), threshes.end());
	idx.get_map().local_insert_by_freqencies(k2_counter, lthreshes);
	BL_BENCH_END(build, "insert", idx.local_size());

	size_t total = idx.size();
	if (comm.rank() == 0) 
		printf("PARSING, FILTER, and INSERT: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "filtered_construct", comm);

	return idx.local_size();
}

#if defined(MIN_MEM) // we'd need to build a k2mer counter incrementally, and then call insert at the end with local k2mer counter data.
					  
template <typename Index>
size_t build_index_thresholded_incremental(::std::vector<::bliss::io::file_data> const &file_data, Index &idx, 
	std::vector<size_t> const &threshes, mxx::comm const &comm,
		std::string k2mer_filename) 
		{
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) 
		printf("PARSING and INSERT incrementally\n");

	// need to build the k2mer counter first using all files
	BL_BENCH_START(build);
  	typename Index::map_type::LocalK2CountMapType k2_counter;
	BL_BENCH_COLLECTIVE_END(build, "init_k2counter", k2_counter.size(), comm);

    using CharIterType = typename ::bliss::io::file_data::const_iterator;
	using SeqParserType = FileParser<CharIterType>;
	using SeqIterType = SplitSeqIterType<CharIterType, FileParser>;
	using KmerParser = ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>;
	using Iter = typename ::bliss::iterator::ContainerConcatenatingIterator<SeqIterType, KmerParser>;

	for (auto x : file_data)
	{

		// initialization
		BL_BENCH_START(build);

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
    	size_t block_size = (free_mem / (8 * sizeof(typename KmerParser::value_type))); // number of elements that can be held in freemem
    	block_size = std::min(block_size, x.getRange().size());

    	if (comm.rank() == 0)
			std::cout << "estimate num elements=" << block_size << ", value_type size=" << sizeof(typename KmerParser::value_type) << " bytes" << std::endl;

		BL_BENCH_COLLECTIVE_END(build, "parse_setup", block_size, comm);

    	//=== copy into array incrementally
		BL_BENCH_START(build);
		idx.get_map().compute_biedge_freqencies_incremental(start, endd, k2_counter, block_size);
		BL_BENCH_COLLECTIVE_END(build, "count_k2mer_incr", k2_counter.size(), comm);
	}
#ifndef NDEBUG
	print_k2mer_frequencies(k2mer_filename, k2_counter, comm);
#endif

	// then filter the k2mers and insert into dbg (no need to touch files again)
	BL_BENCH_START(build);
	std::vector<size_t> lthreshes(threshes.begin(), threshes.end());
	idx.get_map().local_insert_by_freqencies(k2_counter, lthreshes);
	BL_BENCH_END(build, "insert", idx.local_size());

	size_t total = idx.size();
	if (comm.rank() == 0) 
		printf("PARSING and INSERT incremental DONE: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "construct_incr", comm);

	return idx.local_size();
}



#endif


void do_benchmark(::std::vector<::bliss::io::file_data> const & file_data, std::string const & out_prefix,
	bool thresholding, bool benchmark, bool LRoptimized, bool compress, bool mpiio, 
	std::vector<size_t> const & threshes, 
	 mxx::comm const & comm) {
	// filename for compacted chain strings
	// string starts with smaller end.  (first K and rev_comp of last K compared)
	// this is a dump of the collected compacted chains.
	std::string compacted_chain_str_filename(out_prefix);
	compacted_chain_str_filename.append("_chain.fasta");

	BL_BENCH_INIT(benchmark);

	BL_BENCH_START(benchmark);
	ChainGraphType chainmap(comm);
	DBGType idx(comm);

	// =================  make compacted simple DBG, so that we can get chain and branch kmers.
	{
		// first read the input

#if defined(MIN_MEM)
		if (thresholding) {
			std::string k2mer_filename(out_prefix);
			k2mer_filename.append("_k2mers.debug");

			build_index_thresholded_incremental(file_data, idx, threshes, comm, k2mer_filename);
		} else {
			build_index_incremental(file_data, idx, comm);
		}
#else
		if (thresholding) {
			std::string k2mer_filename(out_prefix);
			k2mer_filename.append("_k2mers.debug");

			build_index_thresholded(file_data, idx, threshes, comm, k2mer_filename);
		} else {
			build_index(file_data, idx, comm);
		}
#endif
		BL_BENCH_COLLECTIVE_END(benchmark, "construct", idx.local_size(), comm);

		// TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet
		//   this is done via read filtering/splitting.

		// ==== make chain map
		BL_BENCH_START(benchmark);
		chainmap.extract_chains(idx);
		//make_chain_map(idx, chainmap, comm);
		BL_BENCH_COLLECTIVE_END(benchmark, "chainmap", chainmap.local_size(), comm);
		// == DONE == make chain map
	} // enforce delete idx.


	// ===== parallel list ranking for chain compaction
	BL_BENCH_START(benchmark);
	size_t iterations = 0;
	if (LRoptimized)
		iterations = chainmap.list_rank_min_update();
	else
		iterations = chainmap.list_rank();
	//auto cycle_node_kmers = list_rank(chainmap, comm);
	BL_BENCH_COLLECTIVE_END(benchmark, "list_rank", iterations, comm);
	// == DONE == parallel list ranking for chain compaction

	{

	// now print chain string - order is destroyed via psort.
		// prepare
		BL_BENCH_START(benchmark);
		ListRankedChainNodeVecType compacted_chain = chainmap.to_ranked_chain_nodes();
		BL_BENCH_COLLECTIVE_END(benchmark, "compact_chain", compacted_chain.size(), comm);

		BL_BENCH_START(benchmark);
		print_chain_string(compacted_chain_str_filename, compacted_chain, comm);
		BL_BENCH_COLLECTIVE_END(benchmark, "chain_str", compacted_chain.size(), comm);
	}


	BL_BENCH_REPORT_MPI_NAMED(benchmark, "benchmark", comm);

}


void do_work(::std::vector<::bliss::io::file_data> const &file_data, std::string const &out_prefix,
	bool thresholding, bool benchmark, bool LRoptimized, bool compress, bool mpiio, 
	std::vector<size_t> const &threshes, 
	mxx::comm const &comm) 
	{
	// filename for compacted chain strings
	// string starts with smaller end.  (first K and rev_comp of last K compared)
	// this is a dump of the collected compacted chains.
	std::string compacted_chain_str_filename(out_prefix);
	compacted_chain_str_filename.append("_chain.fasta");

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

	BL_BENCH_INIT(work);

	BL_BENCH_START(work);
	ChainGraphType chainmap(comm);
	CountDBGType idx(comm);

	// =================  make compacted simple DBG, so that we can get chain and branch kmers.

#if defined(MIN_MEM)
	if (thresholding)
	{
		std::string k2mer_filename(out_prefix);
		k2mer_filename.append("_k2mers.debug");

		build_index_thresholded_incremental(file_data, idx, threshes, comm, k2mer_filename);
	}
	else
	{
		build_index_incremental(file_data, idx, comm);
	}
#else
	if (thresholding)
	{
		std::string k2mer_filename(out_prefix);
		k2mer_filename.append("_k2mers.debug");

		build_index_thresholded(file_data, idx, threshes, comm, k2mer_filename);
	}
	else
	{
		build_index(file_data, idx, comm);
	}
#endif
	BL_BENCH_COLLECTIVE_END(work, "construct", idx.local_size(), comm);

		// TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet
		//   this is done via read filtering/splitting.

  // ================== Do stats and checks
#ifndef NDEBUG  
  // ====== print edge histogram
  BL_BENCH_START(work);
  if (comm.rank() == 0)
  	printf("rank 0 checking (thresholded) index\n");
  print_edge_histogram(idx, comm);
  check_index(idx, comm);
  BL_BENCH_COLLECTIVE_END(work, "histo", idx.local_size(), comm);
  // == DONE == make compacted simple DBG

	{
	BL_BENCH_START(work);
		std::string graph_filename(out_prefix);
		graph_filename.append("_graph.all.nodes");
	print_graph_edge_frequencies(graph_filename, idx, comm);
	BL_BENCH_COLLECTIVE_END(work, "print graph", idx.local_size(), comm);
	}
#endif

  // == PRINT == prep branch for printing - here ONLY BECAUSE WE ARE DISCARDING IDX AFTER MAKING CHAINMAPS
	if (!benchmark)
	{
		BL_BENCH_START(work);
			std::string branch_filename(out_prefix);
			branch_filename.append("_branch.edges");
		print_branch_edge_frequencies(branch_filename, idx, comm);
		BL_BENCH_COLLECTIVE_END(work, "print branch edges", idx.local_size(), comm);
	}
	{
		BL_BENCH_START(work);
			std::string branch_fasta_filename(out_prefix);
			branch_fasta_filename.append("_branch.fasta");
		print_branch_fasta(branch_fasta_filename, idx, comm);
		BL_BENCH_COLLECTIVE_END(work, "print branch fasta", idx.local_size(), comm);
	}


	// ==== make chain map
	BL_BENCH_START(work);
	chainmap.extract_chains(idx);
	//make_chain_map(idx, chainmap, comm);
	BL_BENCH_COLLECTIVE_END(work, "chainmap", chainmap.local_size(), comm);
	// == DONE == make chain map

#ifndef NDEBUG  
{
	BL_BENCH_START(work);
	std::string chain_biedge_filename2(out_prefix);
	chain_biedge_filename2.append("_chainmap_uncomp.debug");
	print_chain_biedges(chain_biedge_filename2, chainmap, comm);
	BL_BENCH_COLLECTIVE_END(work, "print_chain_biedge_uncomp", chainmap.local_size(), comm);
}
#endif

	// ===== parallel list ranking for chain compaction
	BL_BENCH_START(work);
	size_t iterations = 0;
	if (LRoptimized)
		iterations = chainmap.list_rank_min_update();
	else
		iterations = chainmap.list_rank();
	//auto cycle_node_kmers = list_rank(chainmap, comm);
	BL_BENCH_COLLECTIVE_END(work, "list_rank", iterations, comm);
	// == DONE == parallel list ranking for chain compaction

	// =============================================================
	// below is for printing.

#ifndef NDEBUG  
	{

		BL_BENCH_START(work);
		std::string chain_biedge_filename(out_prefix);
		chain_biedge_filename.append("_chain.nodes.debug");
		print_chain_biedges(chain_biedge_filename, chainmap, comm);
		BL_BENCH_COLLECTIVE_END(work, "print_chain_biedge", chainmap.local_size(), comm);
	}
	{
		// =============================================================
		// generate chain_summaries
		BL_BENCH_START(work);
		auto summaries = chainmap.to_summarized_chains();
		BL_BENCH_COLLECTIVE_END(work, "chain_summaries", summaries.size(), comm);

		BL_BENCH_START(work);
		std::string chain_summary_filename(out_prefix);
		chain_summary_filename.append("_chains.summary");
		print_chain_summaries(chain_summary_filename, summaries, comm);
		BL_BENCH_COLLECTIVE_END(work, "print_chain_summaries", summaries.size(), comm);
	}
#endif

	{
		// =========== remove cycles and isolated
		BL_BENCH_START(work);
		auto cycle_kmers = chainmap.get_cycle_node_kmers();
		idx.erase(cycle_kmers);
		idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
		BL_BENCH_COLLECTIVE_END(work, "remove cycles/isolated/etc", idx.local_size(), comm);

#ifndef NDEBUG  
		BL_BENCH_START(work);
		std::string graph_filename(out_prefix);
		graph_filename.append("_graph.no_cycle.nodes");
		print_graph_edge_frequencies(graph_filename, idx, comm);
		BL_BENCH_COLLECTIVE_END(work, "print graph", idx.local_size(), comm);

		BL_BENCH_START(work);
		if (comm.rank() == 0)
			printf("rank 0 checking cycle removed index\n");
		print_edge_histogram(idx, comm);
		check_index(idx, comm);
		BL_BENCH_COLLECTIVE_END(work, "histo", idx.local_size(), comm);
#endif
	}

	// done with first compaction





	if (!benchmark)
	{
		// == PRINT == valid k-mers files
		if (thresholding)
		{

#if (pPARSER == FASTA)
			if (comm.rank() == 0)
				printf("WARNING: outputting first/last valid kmer position for each read is supported for FASTQ format only.\n");
#elif (pPARSER == FASTQ)

			BL_BENCH_START(work);
			for (size_t i = 0; i < file_data.size(); ++i)
			{
				std::string fn(out_prefix);

				fn.append(".");
				fn.append(std::to_string(i));
				fn.append(".valid");

				print_valid_kmer_pos_in_reads(fn, file_data[i], idx, comm);
			}
			BL_BENCH_COLLECTIVE_END(work, "print_valid_kmer_pos", file_data.size(), comm);
#endif
		}
	}

	// == print ===  prepare for printing compacted chain and frequencies

	FreqMapType freq_map(comm);
	{
		// prepare
		BL_BENCH_START(work);
		ListRankedChainNodeVecType compacted_chain = chainmap.to_ranked_chain_nodes();
		BL_BENCH_COLLECTIVE_END(work, "compact_chain", compacted_chain.size(), comm);

#ifndef NDEBUG
		if (compress)
		{
			// === and compress.
			BL_BENCH_START(work);
			::std::vector<::std::string> compressed_chain = chainmap.to_compressed_chains();
			BL_BENCH_COLLECTIVE_END(work, "compress_chains", compressed_chain.size(), comm);

			BL_BENCH_START(work);
			// compressed chain
			std::string compressed_chain_filename(out_prefix);
			compressed_chain_filename.append("_compressed_chain.debug");
			size_t out_size = print_compressed_chains(compressed_chain_filename, compressed_chain, comm);
			BL_BENCH_COLLECTIVE_END(work, "chain_str", out_size, comm);
		}
#endif

		// do this first because we need the order of compacted chain to be same as hashed distribution.

		if (!benchmark)
		{
			// compute freq map
			BL_BENCH_START(work);
	#if defined(MIN_MEM)
			compute_freq_map_incremental(compacted_chain, idx, freq_map, comm);
	#else
			compute_freq_map(compacted_chain, idx, freq_map, comm);
	#endif			
			BL_BENCH_COLLECTIVE_END(work, "chain_freqs", freq_map.local_size(), comm);
		}

		// now print chain string - order is destroyed via psort.
		BL_BENCH_START(work);
		std::string compacted_chain_str_filename2(out_prefix);
		compacted_chain_str_filename2.append("_chain.fasta");
		print_chain_string(compacted_chain_str_filename2, compacted_chain, comm);
		BL_BENCH_COLLECTIVE_END(work, "chain_str", compacted_chain.size(), comm);

	// compressed chain
#ifndef NDEBUG
		BL_BENCH_START(work);
		std::string compacted_chain_kmers_filename2(out_prefix);
		compacted_chain_kmers_filename2.append("_chain.components");
		print_chain_nodes(compacted_chain_kmers_filename2, compacted_chain, comm);
		BL_BENCH_COLLECTIVE_END(work, "chain_node", compacted_chain.size(), comm);
#endif

	} // ensure release compacted chain

	if (!benchmark)
	{
		// search in chainmap to find canonical termini.
		BL_BENCH_START(work);
		ChainVecType chain_rep = chainmap.find_if(::bliss::debruijn::filter::chain::IsCanonicalTerminusOrIsolated());
		BL_BENCH_COLLECTIVE_END(work, "chain rep", chain_rep.size(), comm);

		// get terminal k-mers and frequency
		// erase everything except for terminal kmers.
		BL_BENCH_START(work);
		//ChainVecType termini = chainmap.find_if(::bliss::debruijn::filter::chain::IsTerminusOrIsolated());
		std::vector<KmerType> chain_internal_kmers = chainmap.get_internal_node_kmers();
		BL_BENCH_COLLECTIVE_END(work, "chain internal", chain_internal_kmers.size(), comm);
		
		BL_BENCH_START(work);
		//remove everything except for chain terminal
		idx.erase(chain_internal_kmers);
		idx.erase_if(bliss::debruijn::filter::graph::IsBranchPoint());
		BL_BENCH_COLLECTIVE_END(work, "erase_non_termini", idx.local_size(), comm);

		// CountDBGType idx2(comm);
		// {

		// 		BL_BENCH_START(work);
		// 		// same distribution (using same hashmap params) as idx, so is_local can be set to true.
		// 		idx2.insert(termini, true);
		// 		assert(idx2.local_size() == termini.size());  // should be 1 to 1.
		// 		BL_BENCH_COLLECTIVE_END(work, "make_terminal_counter", termini.size(), comm);

		// 		// get the edges counts for these kmers.
		// 		BL_BENCH_START(work);
		// 		count_edges(file_data, selected_edges, thresholding, idx, comm);
		// 		BL_BENCH_COLLECTIVE_END(work, "terminal_edge_freq", idx2.local_size(), comm);
		// } // ensure delete kmers.

#ifndef NDEBUG
		BL_BENCH_START(work);
		print_chain_frequencies(compacted_chain_ends_filename, chain_rep, idx, freq_map, comm);
		BL_BENCH_COLLECTIVE_END(work, "print_chain_freq", chain_rep.size(), comm);
#endif
	}

	// release chainmap
	if (!benchmark) 
	{
	BL_BENCH_START(work);
	chainmap.clear();
	BL_BENCH_COLLECTIVE_END(work, "chainmap_reset", chainmap.local_size(), comm);
	}
	BL_BENCH_REPORT_MPI_NAMED(work, "work", comm);
}

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv)
{

	//////////////// init logging
	LOG_INIT();

	//////////////// initialize MPI and openMP

	mxx::env e(argc, argv);
	mxx::comm comm;

	if (comm.rank() == 0)
		printf("EXECUTING %s\n", argv[0]);

	comm.barrier();

	//////////////// parse parameters

	//////////////// parse parameters
	std::vector<std::string> filenames;
	std::string filename;

	std::string out_prefix;
	out_prefix.assign("./output");

	  // thresholds capped at uint32_t, but uses size_t for representation.
	std::vector<size_t> threshes(6, ((sizeof(CountType) > 4) ? ::std::numeric_limits<uint32_t>::max() : ::std::numeric_limits<CountType>::max()) + 1);

	bool thresholding = false;
	bool benchmark = false;
	bool LRoptimized = false;
	bool compress = false;
  	bool mpiio = false;

	//  std::string queryname(filename);
	//  int sample_ratio = 100;

	// Wrap everything in a try block.  Do this every time,
	// because exceptions will be thrown for problems.
	try
	{

		// Define the command line object, and insert a message
		// that describes the program. The "Command description message"
		// is printed last in the help text. The second argument is the
		// delimiter (usually space) and the last one is the version number.
		// The CmdLine object parses the argv array based on the Arg objects
		// that it contains.
		TCLAP::CmdLine cmd("Parallel de bruijn graph compaction", ' ', "0.5");

		// MPI friendly commandline output.
		::bliss::utils::tclap::MPIOutput cmd_output(comm);
		cmd.setOutput(&cmd_output);

		// Define a value argument and add it to the command line.
		// A value arg defines a flag and a type of value that it expects,
		// such as "-n Bishop".
		TCLAP::ValueArg<std::string> outputArg("O", "output_prefix", "Prefix for output files, including directory", false, "", "string", cmd);

		TCLAP::SwitchArg threshArg("T", "thresholding", "on/off for thresholding", cmd, false);
		TCLAP::MultiArg<size_t> lowerThreshArg("L", "lower_thresholds",
			"Lower frequency thresholds, inclusive. Single value maps to K1mer. 3 values correspond to kmer, k1mer, k2mer.", false, 
			"uint16", cmd);
		TCLAP::MultiArg<size_t> upperThreshArg("U", "upper_thresholds",
			"Upper frequency thresholds, exclusive. Single value maps to K1mer. 3 values correspond to kmer, k1mer, k2mer.", false, 
			"uint16", cmd);

		TCLAP::SwitchArg benchmarkArg("B", "benchmark", "on/off for benchmarking (no file output)", cmd, false);

		TCLAP::SwitchArg lrOptimizeArg("R", "list_rank_opt", "on/off for list ranking optimization", cmd, false);

		TCLAP::SwitchArg compressArg("C", "compress", "on/off for in-mem compressed chains", cmd, false);

		TCLAP::SwitchArg mpiioArg("M", "mpiio", "on - use mpiio for input.  off - use posix io", cmd, false);

		//    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);
		TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names", false, "string", cmd);

		// Parse the argv array.
		cmd.parse(argc, argv);

		filenames = fileArg.getValue();
		out_prefix = outputArg.getValue();

		// ==== thresholding
		thresholding = threshArg.getValue();

		auto lower = lowerThreshArg.getValue();
		auto upper = upperThreshArg.getValue();

		if (lower.size() == 1)
		{
			threshes[2] = ::std::min(lower[0], threshes[2]);

      			threshes[0] = threshes[4] = 0;
		}
		else if (lower.size() == 3)
		{
			threshes[0] = ::std::min(lower[0], threshes[0]);
			threshes[2] = ::std::min(lower[1], threshes[2]);
			threshes[4] = ::std::min(lower[2], threshes[4]);
		}
		if (upper.size() == 1)
		{
			threshes[3] = ::std::min(upper[0], threshes[3]);
		}
		else if (upper.size() == 3)
		{
			threshes[1] = ::std::min(upper[0], threshes[1]);
			threshes[3] = ::std::min(upper[1], threshes[3]);
			threshes[5] = ::std::min(upper[2], threshes[5]);
		}

		// ====
		benchmark = benchmarkArg.getValue();
		LRoptimized = lrOptimizeArg.getValue();
		compress = compressArg.getValue();
    	mpiio = mpiioArg.getValue();
	}
	catch (TCLAP::ArgException &e) // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		exit(-1);
	}

	if (comm.rank() == 0)
		std::cout << "Parallel de bruijn graph compaction v0.5" << std::endl;

	if (thresholding && (comm.rank() == 0))
	{
	  std::cout << "THRESHOLDS: ";
	  std::cout << "k0:  " << threshes[0] << "-" << threshes[1] << std::endl;
    std::cout << "k1:  " << threshes[2] << "-" << threshes[3] << std::endl;
    std::cout << "k2:  " << threshes[4] << "-" << threshes[5] << std::endl;
	}

//#if (pPARSER == FASTA)
//    if (thresholding)
//      throw std::invalid_argument("ERROR: FASTA version of debruijn graph compaction does not currently support filtering.");
//#endif

	if (filenames.size() == 0)
	{
		filename.assign(PROJ_SRC_DIR);

#if (pPARSER == FASTA)
		filename.append("/test/data/test.debruijn.small.fasta");
#elif (pPARSER == FASTQ)
		filename.append("/test/data/test.debruijn.small.fastq");
#endif

		filenames.push_back(filename);
	}

	BL_BENCH_INIT(app);

	// ================  read and get file
	BL_BENCH_START(app);
	::std::vector<::bliss::io::file_data> file_data = open_files(filenames, comm, mpiio);
	BL_BENCH_COLLECTIVE_END(app, "read_file", file_data.size(), comm);
	// == DONE == reading

	BL_BENCH_START(app);
	// if (benchmark) {

	// 	do_benchmark(file_data, out_prefix, thresholding, benchmark, LRoptimized, compress, mpiio, threshes, comm);

	// } else {

		do_work(file_data, out_prefix, thresholding, benchmark, LRoptimized, compress, mpiio, threshes, comm);

	// }
	BL_BENCH_COLLECTIVE_END(app, "processing", file_data.size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(app, "app", comm);

	// mpi cleanup is automatic

	return 0;
}
