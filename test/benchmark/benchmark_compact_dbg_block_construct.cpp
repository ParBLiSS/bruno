// TODO: the logic in this file needs to be updated.

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

	using CountType = uint8_t;

#else   // MIN_MEM not defined.  use uint64_t.  this is same as compact_debruijn_graph_refactor.cpp

	using KmerWordType = uint64_t;  // matches system architecture.

	#if defined(pK)
		using KmerType = bliss::common::Kmer<pK, Alphabet, KmerWordType>;
	#else
		using KmerType = bliss::common::Kmer<31, Alphabet, KmerWordType>;
	#endif

	using CountType = uint16_t;
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

using CountMapType = ::dsc::counting_densehash_map<KmerType, CountType,
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

#include "../common/compact_dbg_build.hpp"

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
		CountType const & lower_thresh, CountType const & upper_thresh,  mxx::comm const & comm) {
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
			results.emplace_back(::bliss::debruijn::biedge::filter::extract_biedges<KmerType>(nodes2));
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
		CountType const & lower_thresh, CountType const & upper_thresh,  mxx::comm const & comm) {
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
		::fsc::back_emplace_iterator<std::vector<EdgeType> > emplace_iter(edges);
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
		total2 += ::bliss::debruijn::biedge::filter::freq_filter_insert_biedges<
				KmerType, Iter, CountMap1Type, Index,
				::fsc::back_emplace_iterator<std::vector<EdgeType> > >(start, endd, counter2, block_size, idx, emplace_iter, comm);
		BL_BENCH_LOOP_PAUSE(build, 2);

		BL_BENCH_LOOP_RESUME(build, 3);
		total3 += edges.size();
		// save the results.
		results.emplace_back(std::move(edges));
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

	if (comm.rank() == 0) printf("WARNING: logic in this program %s needs to be updated to latest\n", argv[0]);

	if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);

	comm.barrier();


	//////////////// parse parameters

	//////////////// parse parameters
	std::vector<std::string> filenames;
	std::string filename;


	std::string out_prefix;
	out_prefix.assign("./output");

//#if (pPARSER == FASTQ)
	CountType lower, upper;
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

		TCLAP::SwitchArg threshArg("T", "thresholding", "on/off for thresholding", cmd, false);
		TCLAP::ValueArg<int64_t> lowerThreshArg("L", "lower_thresh", "Lower Threshold for Kmer and Edge frequency", false, 0, "uint16", cmd);
		TCLAP::ValueArg<int64_t> upperThreshArg("U", "upper_thresh", "Upper Threshold for Kmer and Edge frequency", false,
				std::numeric_limits<CountType>::max(), "uint16", cmd);

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
    lower = ::std::min(static_cast<size_t>(lowerThreshArg.getValue()),
		static_cast<size_t>(::std::numeric_limits<CountType>::max()));
    upper = ::std::min(static_cast<size_t>(upperThreshArg.getValue()),
		static_cast<size_t>(::std::numeric_limits<CountType>::max()));
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
			selected_edges = build_index_thresholded_incremental(file_data, idx, lower, upper, comm);
//#else
//			// TODO: VERIFY FASTA is working. THIS IS WORKING.  FASTA does not support thresholded index build because fasta reader is not yet handling overlaps correctly when dealing with k+2-mers.
//			if (comm.rank() == 0) printf("ERROR: FASTA files does not yet support pre-filter by frequency.\n");
//#endif

			::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> >().swap(selected_edges);
			

		} else {
			build_index_incremental(file_data, idx, comm);
		}
#else
		if (thresholding) {
			selected_edges = build_index_thresholded(file_data, idx, lower, upper, comm);
			::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> >().swap(selected_edges);
		} else {
			build_index(file_data, idx, comm);
		}
#endif


		BL_BENCH_COLLECTIVE_END(app, "construct", idx.local_size(), comm);
		// TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet
		//   this is done via read filtering/splitting.

		// ==== make chain map
		BL_BENCH_START(app);
		chainmap.extract_chains(idx);
		//make_chain_map(idx, chainmap, comm);
		BL_BENCH_COLLECTIVE_END(app, "chainmap", chainmap.local_size(), comm);
		// == DONE == make chain map
	} // enforce delete idx.

	if (benchmark) {
		idx.get_map().clear();
		idx.get_map().reserve(0);
	}

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

	{

    // now print chain string - order is destroyed via psort.
    if (!compress) {
      // prepare
      BL_BENCH_START(app);
      ListRankedChainNodeVecType compacted_chain = chainmap.to_ranked_chain_nodes();
      BL_BENCH_COLLECTIVE_END(app, "compact_chain", compacted_chain.size(), comm);

      BL_BENCH_START(app);
      print_chain_string(compacted_chain_str_filename, compacted_chain, comm);
      BL_BENCH_COLLECTIVE_END(app, "chain_str", compacted_chain.size(), comm);
    } else {
      // === and compress.
      BL_BENCH_START(app);
      ::std::vector<::std::string> compressed_chain = chainmap.to_compressed_chains();
      BL_BENCH_COLLECTIVE_END(app, "compress_chains", compressed_chain.size(), comm);

      BL_BENCH_START(app);
      size_t out_size = print_compressed_chains(compressed_chain_filename, compressed_chain, comm);
      BL_BENCH_COLLECTIVE_END(app, "chain_str", out_size, comm);
    }
	}

	// =============================================================
	// below is for printing.

	BL_BENCH_REPORT_MPI_NAMED(app, "app", comm);


	// mpi cleanup is automatic

	return 0;

}

