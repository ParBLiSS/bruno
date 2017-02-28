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
 * goal: minimize memory usage during kmer counting and node filtering phase.
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
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
#elif (pDNA == 4)
using Alphabet = bliss::common::DNA;
#endif


#if defined(pK)
using KmerType = bliss::common::Kmer<pK, Alphabet, uint16_t>;
#else
using KmerType = bliss::common::Kmer<31, Alphabet, uint16_t>;
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

using CountType = uint16_t;
using CountDBGMapType = ::bliss::debruijn::graph::count_hash_debruijn_graph_map<KmerType, CountType>;
using CountDBGType = ::bliss::debruijn::graph::debruijn_graph<KmerType, CountType>;

using ChainNodeType = ::bliss::debruijn::simple_biedge<KmerType>;

//template <typename K>
//using ChainMapParams = typename DBGType::map_params_template<K>;
//
//using ChainMapType = ::dsc::densehash_map<KmerType, ChainNodeType,
//		ChainMapParams,
//		 ::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

using ChainMapType = ::bliss::debruijn::graph::debruijn_chain_graph<KmerType>;

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


std::string get_error_string(std::string const & filename, std::string const & op_name, int const & return_val, mxx::comm const & comm) {
	char error_string[BUFSIZ];
	int length_of_error_string, error_class;
	std::stringstream ss;

	MPI_Error_class(return_val, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);

	ss << "ERROR in mpiio: rank " << comm.rank() << " " << op_name << " " << filename << " error: " << error_string << std::endl;
	return ss.str();
}

std::string get_error_string(std::string const & filename, std::string const & op_name, int const & return_val, MPI_Status const & stat, mxx::comm const & comm) {
	char error_string[BUFSIZ];
	int length_of_error_string, error_class;
	std::stringstream ss;

	MPI_Error_class(return_val, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);

	ss << "ERROR in mpiio: rank " << comm.rank() << " " << op_name << " " << filename << " error: " << return_val << " [" << error_string << "]";

	//		// status.MPI_ERROR does not appear to be decodable by error_class.  google search did not find how to decode it.
	//		MPI_Error_class(stat.MPI_ERROR, &error_class);
	//		MPI_Error_string(error_class, error_string, &length_of_error_string);

	ss << " MPI_Status error: [" << stat.MPI_ERROR << "]" << std::endl;

	return ss.str();
}


void write_mpiio(std::string const & filename, const char* data, size_t len, mxx::comm const & comm ) {
	// TODO: subcommunicator to work with only nodes that have data.

	/// MPI file handle
	MPI_File fh;

	int res = MPI_File_open(comm, const_cast<char *>(filename.c_str()), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if (res != MPI_SUCCESS) {
		throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "open", res, comm));
	}

	res = MPI_File_set_size(fh, 0);
	if (res != MPI_SUCCESS) {
		throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "truncate", res, comm));
	}

	// ensure atomicity is turned off
	MPI_File_set_atomicity(fh, 1);

	// get the global offset.
	size_t global_offset = ::mxx::exscan(len, comm);

	size_t step = (0x1 << 30);
	size_t iterations = (len + step - 1) / step;

	//std::cout << "rank " << comm.rank() << " mpiio write offset is " << global_offset << " len " << len << " iterations " << iterations << ::std::endl;

	// get the maximum number of iterations
	iterations = ::mxx::allreduce(iterations, [](size_t const & x, size_t const & y){
		return (x >= y) ? x : y;
	}, comm);

	size_t remainder = len;
	size_t curr_step = step;
	MPI_Status stat;
	int count = 0;
	for (size_t i = 0; i < iterations; ++i) {
		curr_step = std::min(remainder, step);

		res = MPI_File_write_at_all( fh, global_offset, const_cast<char*>(data), curr_step, MPI_BYTE, &stat);

		if (res != MPI_SUCCESS)
			throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "write", res, stat, comm));

		res = MPI_Get_count(&stat, MPI_BYTE, &count);
		if (res != MPI_SUCCESS)
			throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "write count", res, stat, comm));

		if (static_cast<size_t>(count) != curr_step) {
			std::stringstream ss;
			ss << "ERROR in mpiio: rank " << comm.rank() << " write error. request " << curr_step << " bytes got " << count << " bytes" << std::endl;

			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
		}

		global_offset += curr_step;
		data += curr_step;
		remainder -= curr_step;
	}

	// close the file when done.
	MPI_File_close(&fh);
}

::std::vector<::bliss::io::file_data> open_files(std::vector<std::string> const & filenames, mxx::comm const & comm, bool mpiio = false) {
	::std::vector<::bliss::io::file_data> file_data;

	BL_BENCH_INIT(open);

	BL_BENCH_START(open);
	size_t total = 0;
	for (auto fn : filenames) {
		if (comm.rank() == 0) printf("READING %s via posix\n", fn.c_str());

		if (mpiio) {
		  ::bliss::io::parallel::mpiio_file<FileParser> fobj(fn, KmerType::size + 1, comm);

      file_data.push_back(fobj.read_file());
		} else {
		  FileReaderType fobj(fn, KmerType::size + 1, comm);

		  file_data.push_back(fobj.read_file());
		}

		total += file_data.back().getRange().size();
	}
	BL_BENCH_COLLECTIVE_END(open, "read", total, comm);

	total = mxx::allreduce(total, comm);
	if (comm.rank() == 0) printf("total size read is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(open, "open", comm);

	return file_data;
}

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


/**
 * @brief pases input file_data, and parse k2mers
 */
::std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge> >
parse_nodes(::bliss::io::file_data const & file_data,
            mxx::comm const & comm) {

  ::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > nodes;

  // the parser needs to collectively parse the records.
  // for kmer at the beginning and end of sequence, generated a padded version.
  ::bliss::io::KmerFileHelper::template parse_file_data<
   ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>,
   FileParser, SplitSeqIterType>(file_data, nodes, comm);

//  std::cout << "rank " << comm.rank() << " read " << nodes.size() << " nodes " << std::endl;
//  std::cout << "rank " << comm.rank() << " first " << nodes.front() << std::endl;
//  std::cout << "      last " << nodes.back() << std::endl;

  return nodes;
}


/**
 * @brief pases input file_data to make k2mers, and filter out k2mers with low frequency k1mer edges.
 * @param file_data
 * @param selected   bit vector indicating whether an edge has the matching frequency.
 */
::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> >
 parse_and_filter_nodes(::bliss::io::file_data const & file_data,
                        ::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> const & selected,
                         mxx::comm const & comm) {

	::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > temp =
			parse_nodes(file_data, comm);

//	::std::vector<KmerType> temp;
//  // the parser needs to collectively parse the records.
//  // for kmer at the beginning and end of sequence, generated a padded version.
//  ::bliss::io::KmerFileHelper::template parse_file_data<::bliss::index::kmer::KmerParser<KmerType>,
//   FileParser, SplitSeqIterType>(file_data, temp, comm);

  // reconstruct kmer biedge pairs, and filter out isolated.
  return ::bliss::debruijn::biedge::filter::reconstruct_filter_nodes(temp, selected);

}


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
		uint8_t const & lower_thresh, uint8_t const & upper_thresh,  mxx::comm const & comm) {
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) printf("PARSING, FILTER, and INSERT\n");

	// k+1-mer count map.  note that it should use the same hash function as Index.
	using K1merType	= ::bliss::common::Kmer<KmerType::size + 1, typename KmerType::KmerAlphabet, typename KmerType::KmerWordType>;
	using CountMap1Type = ::dsc::saturating_counting_densehash_map<K1merType, uint8_t,
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

template <typename Index>
::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> >
build_index_thresholded_incremental(::std::vector<::bliss::io::file_data> const & file_data, Index & idx,
		uint8_t const & lower_thresh, uint8_t const & upper_thresh,  mxx::comm const & comm) {
	BL_BENCH_INIT(build);

	if (comm.rank() == 0) printf("PARSING, FILTER, and INSERT\n");


	// ========  count the k+1-mers.
	using K1merType	= ::bliss::common::Kmer<KmerType::size + 1, typename KmerType::KmerAlphabet, typename KmerType::KmerWordType>;
	// k+1-mer count map.  note that it should use the same hash function as Index.
	using CountMap1Type = ::dsc::saturating_counting_densehash_map<K1merType, uint8_t,
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


//#endif


template <typename Index>
void check_index(Index const & idx, mxx::comm const & comm) {

	// ============== testing to ensure that all the possible edges are present.
	std::vector<KmerType> query;

	std::vector<KmerType> neighbors;
	neighbors.reserve(5);

	auto cc = idx.get_map().get_local_container();
	for (auto it = cc.begin(); it != cc.end(); ++it) {

		//std::cout << "kmer: " << ::bliss::utils::KmerUtils::toASCIIString(it->first) << " edge " << it->second << std::endl;

		neighbors.clear();
		it->second.get_out_neighbors(it->first, neighbors);
		query.insert(query.end(), neighbors.begin(), neighbors.end());

		neighbors.clear();
		it->second.get_in_neighbors(it->first, neighbors);
		query.insert(query.end(), neighbors.begin(), neighbors.end());
	}

	// =============== check to see if index is superset of query.  (count should have every result entry showing 1.)
	{
	  auto lquery = query;
	  auto counts = idx.count(lquery);

	  auto absent_end = std::partition(counts.begin(), counts.end(), [](std::pair<KmerType, size_t> const & x){
		  return x.second == 0;
	  });
	  printf(" total query = %lu, unique query = %lu, unique absent = %lu\n", query.size(), counts.size(), std::distance(counts.begin(), absent_end));

	  for (auto it = counts.begin(); it != absent_end; ++it) {
		  std::cout << "absent k-mer " << ::bliss::utils::KmerUtils::toASCIIString(it->first) << std::endl;
	  }
	  assert( std::distance(counts.begin(), absent_end) == 0);
	}

	comm.barrier();

	// =============== check to see if query is superset of index.  (erase should result in empty)
	{
		auto lquery = query;
		Index idx_copy(comm);
		idx_copy.get_map().get_local_container().insert(idx.get_map().get_local_container().begin(), idx.get_map().get_local_container().end());

		size_t erased = idx_copy.get_map().erase(lquery);

		printf(" total query = %lu, erased = %lu, remaining = %lu\n", query.size(), erased, idx.local_size());
		assert(idx_copy.size() == 0);
	}


}

template <typename Index>
void print_edge_histogram(Index const & idx, mxx::comm const & comm) {
	// get histogram for the edge types in debruijn graph
	if (comm.rank() == 0) printf("HISTOGRAM\n");

	BL_BENCH_INIT(histo);

	BL_BENCH_START(histo);
	// then compute histogram
	::bliss::debruijn::graph::print_compact_multi_biedge_histogram(idx.get_map().get_local_container().begin(),
			idx.get_map().get_local_container().end(), comm);
	BL_BENCH_COLLECTIVE_END(histo, "histogram", idx.local_size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(histo, "insert", comm);
}


/**
 * @brief counting just the kmers that have been selected, accounting for selected edge (based on frequency, e.g.)
 * @details:  insert the selected kmers first into count dbg with 0 freqs.
 * 			  then read the file for kmers, attaching the selected edges, and insert into count to get the actual counts.
 * 			  return the count dbg.
 */
void count_edges(::std::vector<::bliss::io::file_data> const & file_data,
		 ::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> > const & selected_edges,
		  bool thresholding,
		 CountDBGType & idx2, mxx::comm const & comm) {

//	BL_BENCH_INIT(count_edge);
//
//	// parser the file data
//	BL_BENCH_START(count_edge);
	::bliss::debruijn::operation::graph::compact_multi_biedge_update<KmerType> updater;

	::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > temp;
	size_t count = 0;
	for (size_t i = 0; i < file_data.size(); ++i) {
//#if (pPARSER == FASTQ)
    if (thresholding)
      parse_and_filter_nodes(file_data[i], selected_edges[i], comm).swap(temp);

    else
//#endif
      parse_nodes(file_data[i], comm).swap(temp);

	  count += idx2.update(temp, false, updater);
	}
//	BL_BENCH_COLLECTIVE_END(count_edge, "update_count", count, comm);
//
//	BL_BENCH_REPORT_MPI_NAMED(count_edge, "edge counts", comm);
}

void count_edges_old(std::vector<KmerType> const & selected,
		::std::vector<::bliss::io::file_data> const & file_data,
		 ::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> > const & selected_edges,
		  bool thresholding,
		 CountDBGType & idx2, mxx::comm const & comm) {

	BL_BENCH_INIT(count_edge);

	using mapped_type = typename CountDBGMapType::mapped_type;

	// pre populate the count index.  SAME DISTRIBUTION AS input (input is from hashmap with same hashed distribution
	BL_BENCH_START(count_edge);
	idx2.get_map().get_local_container().resize(selected.size());   // reserve size

	// insert empty entry first.
	for (auto x : selected) {
		idx2.get_map().get_local_container().insert(std::make_pair(x, mapped_type()));
	}
	BL_BENCH_COLLECTIVE_END(count_edge, "alloc_count", idx2.local_size(), comm);

	// parser the file data
	BL_BENCH_START(count_edge);
	::bliss::debruijn::operation::graph::compact_multi_biedge_update<KmerType> updater;

	::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > temp;
	size_t count = 0;
	for (size_t i = 0; i < file_data.size(); ++i) {
//#if (pPARSER == FASTQ)
    if (thresholding)
      parse_and_filter_nodes(file_data[i], selected_edges[i], comm).swap(temp);

    else
//#endif
      parse_nodes(file_data[i], comm).swap(temp);

		count += idx2.get_map().update(temp, false, updater);
	}
	BL_BENCH_COLLECTIVE_END(count_edge, "update_count", count, comm);

	assert(selected.size() == idx2.local_size());

	BL_BENCH_REPORT_MPI_NAMED(count_edge, "edge counts", comm);
}



void print_branch_edge_frequencies(
		std::string const & filename,
		CountDBGType const & idx2,
		mxx::comm const & comm) {

	if (comm.rank() == 0) printf("PRINT BRANCHES\n");
	BL_BENCH_INIT(branch_print);

	// then find branches.
	BL_BENCH_START(branch_print);
	std::vector<typename CountDBGType::mutable_value_type> branch_pts =
			idx2.find_if(::bliss::debruijn::filter::graph::IsBranchPoint());
	BL_BENCH_COLLECTIVE_END(branch_print, "get_branch_counts", branch_pts.size(), comm);

	// sort the branches
	int has_data = (branch_pts.size() == 0) ? 0 : 1;
	int all_has_data = mxx::allreduce(has_data, comm);
	if (all_has_data > 0) {
		// global sort
		BL_BENCH_START(branch_print);
		mxx::comm subcomm = (all_has_data == comm.size()) ? comm.copy() : comm.split(has_data);
		if (has_data == 1) {
			mxx::sort(branch_pts.begin(), branch_pts.end(), [](typename CountDBGType::mutable_value_type const & x,
					typename CountDBGType::mutable_value_type const & y){
				return x.first < y.first;
			}, subcomm);
		}
		BL_BENCH_COLLECTIVE_END(branch_print, "psort branches", branch_pts.size(), comm);   // this is for ordered output.
	}

	// and print.
	BL_BENCH_START(branch_print);

	std::stringstream ss;
	ss.clear();
	std::for_each(branch_pts.begin(), branch_pts.end(),
			::bliss::debruijn::operation::graph::print_graph_node<KmerType>(ss));
	write_mpiio(filename, ss.str().c_str(), ss.str().length(), comm);

	BL_BENCH_COLLECTIVE_END(branch_print, "print branches (4)", branch_pts.size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(branch_print, "branch_print", comm);

}

void print_branch_fasta(
    std::string const & filename,
    CountDBGType const & idx2,
    mxx::comm const & comm) {

  if (comm.rank() == 0) printf("PRINT BRANCH KMERS\n");
  BL_BENCH_INIT(branch_print);

  // then find branches.
  BL_BENCH_START(branch_print);
  std::vector<typename CountDBGType::mutable_value_type> branch_pts =
      idx2.find_if(::bliss::debruijn::filter::graph::IsBranchPoint());
  BL_BENCH_COLLECTIVE_END(branch_print, "get_branch_counts", branch_pts.size(), comm);

  // sort the branches
  int has_data = (branch_pts.size() == 0) ? 0 : 1;
  int all_has_data = mxx::allreduce(has_data, comm);
  if (all_has_data > 0) {
    // global sort
    BL_BENCH_START(branch_print);
    mxx::comm subcomm = (all_has_data == comm.size()) ? comm.copy() : comm.split(has_data);
    if (has_data == 1) {
      mxx::sort(branch_pts.begin(), branch_pts.end(),
    		  [](typename CountDBGType::mutable_value_type const & x,
          typename CountDBGType::mutable_value_type const & y){
        return x.first < y.first;
      }, subcomm);
    }
    BL_BENCH_COLLECTIVE_END(branch_print, "psort branches", branch_pts.size(), comm);   // this is for ordered output.
  }

  // and print.
  BL_BENCH_START(branch_print);

  std::stringstream ss;
  ss.clear();
  std::for_each(branch_pts.begin(), branch_pts.end(),
      ::bliss::debruijn::operation::graph::print_graph_node_fasta<KmerType>(ss));
  write_mpiio(filename, ss.str().c_str(), ss.str().length(), comm);

  BL_BENCH_COLLECTIVE_END(branch_print, "print branches (6)", branch_pts.size(), comm);

  BL_BENCH_REPORT_MPI_NAMED(branch_print, "branch_print", comm);

}

/// printsthe first and last valid k-mer positions in each read.
/// this is for rahul's distance constraints, no no
template <typename FP = FileParser<typename ::bliss::io::file_data::container::const_iterator>, typename Index,
		typename std::enable_if<std::is_same<FP,
				::bliss::io::FASTQParser<typename ::bliss::io::file_data::container::const_iterator> >::value, int>::type = 1>
void print_valid_kmer_pos_in_reads(std::string const & filename,
		::bliss::io::file_data const & fdata, Index & idx,
		mxx::comm const & comm) {
	// parse through the reads
	BL_BENCH_INIT(valid_print);

	BL_BENCH_START(valid_print);
	std::vector<KmerType> kmers;
  ::bliss::debruijn::lex_less<KmerType> canonical;


	// also get the read's starting positions in the output kmer vector.  first is the length of the sequence, second is whether it contains N or not.
  // true indicate that this is a valid read (without N).
	std::vector<std::pair<size_t, bool> >  local_offsets;

	using LocalCountMapType = ::dsc::counting_densehash_map<KmerType, size_t,
			FreqMapParams,
			::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

	typename LocalCountMapType::local_container_type counter;
	typename LocalCountMapType::local_container_type::const_iterator it;

	std::stringstream ss;

	kmers.clear();
	local_offsets.clear();

	ss.str(std::string());
	BL_BENCH_COLLECTIVE_END(valid_print, "init", fdata.getRange().size(), comm);

	// note that we are not using a split sequence iterator here, or a filtering sequence iterator, sicne we NEED to identify the actual first entry.

	BL_BENCH_START(valid_print);
	// get the kmers
  ::bliss::io::KmerFileHelper::template parse_file_data<::bliss::index::kmer::KmerParser<KmerType>, FileParser, SeqIterType>(fdata, kmers, comm);
	std::transform(kmers.begin(), kmers.end(), kmers.begin(), canonical);
	BL_BENCH_COLLECTIVE_END(valid_print, "parse kmers", kmers.size(), comm);


	BL_BENCH_START(valid_print);
	// get the read lengths
  ::bliss::io::KmerFileHelper::template parse_file_data_old<::bliss::debruijn::ReadLengthParser<KmerType>, FileParser, SeqIterType>(fdata, local_offsets, comm);
	// prefix scan to get the offsets
	for (size_t i = 1; i < local_offsets.size(); ++i) {
		local_offsets[i].first += local_offsets[i-1].first;
	}
	// global prefix scan.  only for procs that have data.
//	::mxx::comm subcomm = comm.split(local_offsets.size() > 0);
//	size_t global_offset = 0;
//	if (local_offsets.size() > 0) {
//		global_offset = ::mxx::exscan(local_offsets.back(), subcomm);
//	}
	BL_BENCH_COLLECTIVE_END(valid_print, "read size", local_offsets.size(), comm);


	BL_BENCH_START(valid_print);
	// for all k_mers, check existence.  use the node's kmers.  put into a local count map
	{
		counter.clear();
		auto results = idx.count(kmers); // idx.find(kmers);
//		std::cout << "rank " << comm.rank() << " result size " << results.size() << " reserved " << results.capacity() << std::endl;
		counter.resize(results.size());
//    std::cout << "rank " << comm.rank() << " before insert counter size " << counter.size() << " bucket count " << counter.bucket_count() <<  std::endl;
		for (auto it = results.begin(); it != results.end(); ++it) {
		  if ((*it).second != 0) counter.insert(*it);
		  //counter.insert(std::make_pair((*it).first, 1));
		}

//		std::cout << "rank " << comm.rank() << " after insert counter size " << counter.size() << " bucket count " << counter.bucket_count() << 	std::endl;
	}
	BL_BENCH_COLLECTIVE_END(valid_print, "local count", counter.size(), comm);


	BL_BENCH_START(valid_print);
	// get the kmers again - earlier kmer vector is scrambled by idx.count.  doing this instead of saving another copy because of space constraints.
	kmers.clear();
  ::bliss::io::KmerFileHelper::template parse_file_data<::bliss::index::kmer::KmerParser<KmerType>, FileParser, SeqIterType>(fdata, kmers, comm);
	std::transform(kmers.begin(), kmers.end(), kmers.begin(), canonical);
	BL_BENCH_COLLECTIVE_END(valid_print, "reparse", kmers.size(), comm);


	BL_BENCH_START(valid_print);
	// linear scan to output the valid positions
	int64_t rstart = 0;
	int64_t rend, vstart, vend;
	for (size_t i = 0; i < local_offsets.size(); ++i) {
		rend = static_cast<int64_t>(local_offsets[i].first);

		vstart = -1;
		vend = -1;

		if (local_offsets[i].second == true) {  // only compute if there is no N.

		//std::cout << "global offset " << global_offset << " local offsets : " << rstart << " - " << rend << std::endl;

      for (int64_t j = rstart; j < rend; ++j) {


        it = counter.find(kmers[j]);

        if (it != counter.end()) {
  //				std::cout << "rank " << comm.rank() << " pos " << j << " rstart " << rstart <<
  //						" query " << bliss::utils::KmerUtils::toASCIIString(kmers[j]) <<
  //						" result " << bliss::utils::KmerUtils::toASCIIString((*it).first) <<
  //						" start count " << (*it).second << std::endl;
          // kmer exists in the debruijn graph
          if ((*it).second > 0) {
            //print the first pos.
            vstart = j - rstart;
            break;
          }
        }
      }

      if (vstart > -1) {
        for (int64_t j = rend - 1; j >= rstart; --j) {
          it = counter.find(kmers[j]);

          if (it != counter.end()) {
  //					std::cout << "rank " << comm.rank() << " pos " << j << " rstart " << rstart <<
  //							" query " << bliss::utils::KmerUtils::toASCIIString(kmers[j]) <<
  //							" result " << bliss::utils::KmerUtils::toASCIIString((*it).first) <<
  //							" end count " << (*it).second << std::endl;
            if ((*it).second > 0) {
              vend = j - rstart;
              break;
            }
          }
        }
      }
		}  // only compute start and end if there is no N.

//		// DEBUG
//		for (int64_t j = rstart; j < rend; ++j) {
//		  it = counter.find(kmers[j]);
//		  std::cout << (it != counter.end() ? "1 " : "0 ");
//		}
//		std::cout << std::endl;
//
		// get ready for next read.
		rstart = rend;

		// print start and end, and
		ss << vstart << "\t" << vend << std::endl;
	}
	BL_BENCH_COLLECTIVE_END(valid_print, "find start-end", local_offsets.size(), comm);


	// write out to file
	BL_BENCH_START(valid_print);
	write_mpiio(filename, ss.str().c_str(), ss.str().length(), comm);

	BL_BENCH_COLLECTIVE_END(valid_print, "print", ss.str().length(), comm);


	BL_BENCH_REPORT_MPI_NAMED(valid_print, "read_pos_print", comm);

}

// DEPRECATED
//ListRankedChainNodeVecType to_compacted_chain(ChainMapType const & chainmap,
//		mxx::comm const & comm) {
//
//	BL_BENCH_INIT(chain_convert);
//
//	BL_BENCH_START(chain_convert);
//	ListRankedChainNodeVecType compacted_chain;
//	compacted_chain.reserve(chainmap.size());
//	::fsc::back_emplace_iterator<ListRankedChainNodeVecType > back_emplacer(compacted_chain);
//
//	//== first transform nodes so that we are pointing to canonical terminus k-mers.
//	std::transform(chainmap.get_local_container().cbegin(), chainmap.get_local_container().cend(), back_emplacer,
//			::bliss::debruijn::operation::chain::to_listranked_chain_node<KmerType>());
//	BL_BENCH_COLLECTIVE_END(chain_convert, "transform chain", chainmap.local_size(), comm);
//
//	BL_BENCH_REPORT_MPI_NAMED(chain_convert, "convert_chain", comm);
//
//	return compacted_chain;
//}

void print_chain_string(std::string const & filename,
		ListRankedChainNodeVecType & compacted_chain,
		mxx::comm const & comm) {
	// ========== construct new graph with compacted chains and junction nodes.
	BL_BENCH_INIT(print_chain_string);

	if (comm.rank() == 0) printf("PRINT CHAIN String\n");

	int has_data = (compacted_chain.size() == 0) ? 0 : 1;
	int all_has_data = mxx::allreduce(has_data, comm);
	if (all_has_data > 0) {
		// global sort
		BL_BENCH_START(print_chain_string);
		mxx::comm subcomm = (all_has_data == comm.size()) ? comm.copy() : comm.split(has_data);
		if (has_data == 1) {
			mxx::sort(compacted_chain.begin(), compacted_chain.end(), ::bliss::debruijn::operation::chain::chain_rep_less<KmerType>(), subcomm);
		}
		BL_BENCH_COLLECTIVE_END(print_chain_string, "psort lmer", compacted_chain.size(), comm);   // this is for constructing the chains


		// print out.
		BL_BENCH_START(print_chain_string);
		if (has_data == 1) {
				std::cout << "rank " << comm.rank() << " printing " << std::endl << std::flush;

			std::stringstream ss;
			std::for_each(compacted_chain.begin(), compacted_chain.end(), ::bliss::debruijn::operation::chain::print_chain_as_fasta<KmerType>(ss));
			// above will produce an extra newline character at the beginning of the first.  below special cases it to not print that character
			if (subcomm.rank() == 0) {
				write_mpiio(filename, ss.str().c_str() + 1, ss.str().length() - 1, subcomm);
			} else {
				write_mpiio(filename, ss.str().c_str(), ss.str().length(), subcomm);
			}
		}
		BL_BENCH_COLLECTIVE_END(print_chain_string, "print chains (3)", compacted_chain.size(), comm);

	}


	BL_BENCH_REPORT_MPI_NAMED(print_chain_string, "convert_chain", comm);
}

//void compress_chains(ListRankedChainNodeVecType & compacted_chains, std::vector<::std::string> & compressed_chains,
//    mxx::comm const & comm) {
//  // ========== construct new graph with compacted chains and junction nodes.
//  BL_BENCH_INIT(compress_chain);
//
//  if (comm.rank() == 0) printf("Compress Chain\n");
//
//  int has_data = (compacted_chains.size() == 0) ? 0 : 1;
//  int all_has_data = mxx::allreduce(has_data, comm);
//  if (all_has_data > 0) {
//    // global sort
//    BL_BENCH_START(compress_chain);
//
//    // distribute the data
//    bool sorted = false;
//    std::vector<size_t> recv_counts =
//        ::dsc::distribute(compacted_chains,
//                          ::bliss::debruijn::operation::chain::chain_node_to_proc<::bliss::debruijn::CanonicalDeBruijnHashMapParams, KmerType>(comm.size()),
//                           sorted, comm);
//    BL_BENCH_COLLECTIVE_END(compress_chain, "distribute nodes", compacted_chains.size(), comm);  // this is for output ordering.
//  } else
//    return;
//
//  // next local sort the data by terminus kmer and position
//  BL_BENCH_START(compress_chain);
//  ::std::sort(compacted_chains.begin(), compacted_chains.end(), ::bliss::debruijn::operation::chain::chain_rep_less<KmerType>());
//  BL_BENCH_COLLECTIVE_END(compress_chain, "sort lmer nodes", compacted_chains.size(), comm);   // this is for constructing the chains
//
//  // compressing.
//  std::stringstream ss;
//  //=== approach is to scan for the end of a chain
//  auto curr = compacted_chains.begin();
//  auto next = curr;
//  BL_BENCH_START(compress_chain);
//  compressed_chains.clear();
//  while (curr != compacted_chains.end()) {
//    // scan forward for start of next chain, using adjacent find.
//    next = ::std::adjacent_find(curr, compacted_chains.end(), [](
//      ::bliss::debruijn::chain::listranked_chain_node<KmerType> const & x,
//       ::bliss::debruijn::chain::listranked_chain_node<KmerType> const & y
//    ){
//      // adjacent_find returns first occurrence of 2 consecutive elements that the predicate evaluates to true.
//      return std::get<1>(x) != std::get<1>(y);
//    });
//
//    if (next != compacted_chains.end()) ++next;  // get the second of the pair, == start of next.
//
//    // now do something with the <curr, next> range.
//    ss.clear();
//    ss.str(std::string());
//    std::for_each(curr, next, ::bliss::debruijn::operation::chain::print_chain_as_fasta<KmerType>(ss, true));  // chain_only
//    compressed_chains.push_back(ss.str());
//
//    // get ready for next segment.
//    curr = next;
//
//  }
//  BL_BENCH_COLLECTIVE_END(compress_chain, "toASCII", compressed_chains.size(), comm);
//
//  BL_BENCH_REPORT_MPI_NAMED(compress_chain, "compress_chain", comm);
//
//}

size_t print_compressed_chains(std::string const & filename,
                             ::std::vector<::std::string> const & compressed_chain,
                              mxx::comm const & comm) {


  // aggregate into a single stringstream object, then print.

  BL_BENCH_INIT(print_compressed_chain);

  if (comm.rank() == 0) printf("PRINT COMPRESSED CHAIN String\n");
  size_t pernode = 0;

  int has_data = (compressed_chain.size() == 0) ? 0 : 1;
  int all_has_data = mxx::allreduce(has_data, comm);
  if (all_has_data > 0) {
    // global sort
    BL_BENCH_START(print_compressed_chain);
    mxx::comm subcomm = (all_has_data == comm.size()) ? comm.copy() : comm.split(has_data);
    if (has_data == 1) {
      std::stringstream ss;
      std::for_each(compressed_chain.begin(), compressed_chain.end(), [&ss](::std::string const & x){
        ss << x << std::endl;
      });

      pernode = ss.str().length();
      write_mpiio(filename, ss.str().c_str(), ss.str().length(), subcomm);
    }
    BL_BENCH_COLLECTIVE_END(print_compressed_chain, "print", pernode, comm);   // this is for constructing the chains
  }

  BL_BENCH_REPORT_MPI_NAMED(print_compressed_chain, "print compressed", comm);

  return pernode;
}

void print_chain_nodes(std::string const & filename,
		ListRankedChainNodeVecType & compacted_chain,
		mxx::comm const & comm) {

	BL_BENCH_INIT(print_chain_nodes);

	//===  print chain nodes (1)
	if (comm.rank() == 0) printf("PRINT CHAIN Nodes\n");

	int has_data = (compacted_chain.size() == 0) ? 0 : 1;
	int all_has_data = mxx::allreduce(has_data, comm);
	if (all_has_data > 0) {
		// global sort
		BL_BENCH_START(print_chain_nodes);
		mxx::comm subcomm = (all_has_data == comm.size()) ? comm.copy() : comm.split(has_data);
		if (has_data == 1) {
		  // sort by node's kmer (not the chain rep and position, since we need to have sorted k-mers.
			mxx::sort(compacted_chain.begin(), compacted_chain.end(), ::bliss::debruijn::operation::chain::chain_node_less<KmerType>(), subcomm);
		}
		BL_BENCH_COLLECTIVE_END(print_chain_nodes, "psort kmer", compacted_chain.size(), comm);  // this is for output ordering.
	}

	// print out.
	BL_BENCH_START(print_chain_nodes);
	{
		std::stringstream ss2;
		std::for_each(compacted_chain.begin(), compacted_chain.end(), ::bliss::debruijn::operation::chain::print_chain_node<KmerType>(ss2));
		write_mpiio(filename, ss2.str().c_str(), ss2.str().length(), comm);

	}
	//      std::cout << ss.str() << std::endl;
	BL_BENCH_COLLECTIVE_END(print_chain_nodes, "print chains (1)", compacted_chain.size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(print_chain_nodes, "convert_chain", comm);

}


/// count kmers in chains.  compacted chain should have same distribution as would be for count index.
void count_kmers(::std::vector<::bliss::io::file_data> const & file_data,
                 ::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> > const & selected_edges,
                  bool thresholding,
		 CountIndexType & count_idx,
		 mxx::comm const & comm) {

	BL_BENCH_INIT(count_kmer);

	// parser the file data and insert in the count index

	BL_BENCH_START(count_kmer);
	::std::vector<::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> > temp1;
	::std::vector< KmerType > temp;
	::fsc::back_emplace_iterator<std::vector<KmerType> > back_emplacer(temp);


	for (size_t i = 0; i < file_data.size(); ++i) {

		// TODO: this part can be simplified? so that memory usage is further reduced?

	  // use the same mechanism as the one for building the graph, so we can count.
//#if (pPARSER == FASTQ)
	  if (thresholding)
		  parse_and_filter_nodes(file_data[i], selected_edges[i], comm).swap(temp1);
	  else
//#endif
	    parse_nodes(file_data[i], comm).swap(temp1);

		// copy out the kmer only.  overlap is k+1 plus any newline chars.  because of the newline chars, not practical to truncate x.
		// this is safer.
		temp.clear();
		size_t counter = 0;
		std::transform(temp1.begin(), temp1.end(), back_emplacer,
				[&comm, &counter](::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> const & y){
//			std::cout << "kmer counting: rank " << comm.rank() << " counter " << counter << " "<< y.first << std::endl;
			++counter;
			return y.first;
		});

		count_idx.insert(temp);   // distributed insert.
	}
	BL_BENCH_COLLECTIVE_END(count_kmer, "insert", count_idx.local_size(), comm);


	BL_BENCH_REPORT_MPI_NAMED(count_kmer, "kmer counts", comm);
}

/// compute frequency of chains.  compacted chain should have same distribution as would be for count index.
void compute_freq_map(ListRankedChainNodeVecType const & compacted_chain,
		CountIndexType const & count_idx,
		FreqMapType & chain_freq_map,
		mxx::comm const & comm) {

//	BL_BENCH_INIT(chain_freq);
//
//	// ==  first compute frequency summary, and store into a reduction map
//	// allocate input
//	std::vector< std::pair<KmerType, FreqSummaryType > > freqs;
//
//	BL_BENCH_START(chain_freq);
//	// extract frequencies.
//	for (auto x : compacted_chain) {
//		// compute the chain rep
//		CountType c = count_idx.get_map().get_local_container().find(std::get<0>(x))->second;
//
//		// new key is the chain rep.
//		freqs.emplace_back(std::get<1>(x), FreqSummaryType(1, c, c, c));
//	}
//	BL_BENCH_COLLECTIVE_END(chain_freq, "get_node_freqs", freqs.size(), comm);
//
//	// create a reduction map
//	BL_BENCH_START(chain_freq);
//	chain_freq_map.insert(freqs);   // collective comm.
//	BL_BENCH_COLLECTIVE_END(chain_freq, "reduce_freq", chain_freq_map.local_size(), comm);
//
//	BL_BENCH_REPORT_MPI_NAMED(chain_freq, "chain counts", comm);


	BL_BENCH_INIT(chain_freq);

	// ==  first compute frequency summary, and store into a reduction map
	// allocate input
	std::vector< std::pair<KmerType, FreqSummaryType > > freqs;

	// estimate the largest amount of memory to use.
	unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);

	// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
	size_t step = (free_mem / (8 * sizeof(std::pair<KmerType, FreqSummaryType >)));  // number of elements that can be held in freemem
  step = std::min(step, compacted_chain.size());

	if (comm.rank() == 0) std::cout << "estimate num elements=" << step << ", value_type size=" <<
			sizeof(std::pair<KmerType, FreqSummaryType >) << " bytes" << std::endl;

	freqs.reserve(step);   // do in steps of 1000000
	size_t nsteps = (compacted_chain.size() + step - 1) / step;
	nsteps = mxx::allreduce(nsteps, [](size_t const & x, size_t const & y){
		return std::max(x, y);
	}, comm);

	BL_BENCH_START(chain_freq);
	auto iter = compacted_chain.begin();
	auto end = compacted_chain.end();
	CountType c = 0;

	::bliss::debruijn::lex_less<KmerType> canonical;

	for (size_t s = 0; s < nsteps; ++s) {

		freqs.clear();

		// extract frequencies.
		for (size_t i = 0; (i < step) && (iter != end); ++i, ++iter) {
			// compute the chain rep
			c = count_idx.get_map().get_local_container().find(canonical(std::get<0>(*iter)))->second;

			// new key is the chain rep.
			freqs.emplace_back(std::get<1>(*iter), FreqSummaryType(1, c, c, c));
			//std::cout << "rank " << comm.rank() << " kmer " << std::get<0>(*iter) << " freq " << c << std::endl;
		}

		// create a reduction map
		chain_freq_map.insert(freqs);   // collective comm.
	}
	BL_BENCH_COLLECTIVE_END(chain_freq, "compute_freq", chain_freq_map.local_size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(chain_freq, "chain freqs", comm);
}


/// compute frequency of chains.  chain_rep should have same distribution as would be for count index.
void print_chain_frequencies(std::string const & filename,
		ChainVecType const & chain_reps,
		CountDBGType const & idx2,
		FreqMapType const & freq_map,
		mxx::comm const & comm ) {

	BL_BENCH_INIT(print_chain_freq);

	// ===========  print compacted chain with left and right
	if (comm.rank() == 0) printf("COMPUTE CHAIN FREQ SUMMARY\n");

	// freq_map, idx2, and chainmap now all have same distribution.
	KmerType L, R, cL, cR;
	using edge_freq_type = std::tuple<KmerType, KmerType,
			std::tuple<CountType, CountType, CountType, CountType>,
			std::tuple<CountType, CountType, CountType, CountType>,
			std::tuple<CountType, CountType, CountType> >;
	std::vector<edge_freq_type> edge_freqs;


	//========= get the R frequencies (from remote) and insert into a local map
	if (comm.rank() == 0) printf("GATHER NON_REP_END EDGE FREQUENCY\n");

	// get query vector
	BL_BENCH_START(print_chain_freq);
	typename CountDBGMapType::local_container_type R_freq_map;
	{
		std::vector<KmerType> R_query;
		for (auto x : chain_reps) {
			if ((::std::get<2>(x.second) == 0) && (::std::get<3>(x.second) == 0)) {
				R_query.emplace_back(x.first);
			} else if (::std::get<2>(x.second) == 0) {
				R_query.emplace_back(std::get<1>(x.second));
			} else if (::std::get<3>(x.second) == 0) {
				R_query.emplace_back(std::get<0>(x.second));
			} else {
				std::cout << "rank " << comm.rank() << " canonical chain terminal not really terminal" << std::endl;
				continue;
			}
			//				  std::cout << "rank " << comm.rank() << " R query " << bliss::utils::KmerUtils::toASCIIString(R_query.back()) << std::endl;
		}
		// do query and insert results into local map.
		auto R_results = idx2.find(R_query);
		R_freq_map.insert(R_results.begin(), R_results.end());
		//			  for (auto x : R_freq_map) {
		//				  std::cout << "rank " << comm.rank() << " R result " << bliss::utils::KmerUtils::toASCIIString(x.first) << std::endl;
		//			  }
	}
	BL_BENCH_COLLECTIVE_END(print_chain_freq, "local_R_freqs", R_freq_map.size(), comm);

	// convert to tuple.
	if (comm.rank() == 0) printf("CREATE CHAIN EDGE FREQUENCIES\n");

	BL_BENCH_START(print_chain_freq);
	bliss::debruijn::lex_less<KmerType> canonical;
	auto compact_edgeL = idx2.get_map().get_local_container().find(KmerType());
  typename CountDBGMapType::local_container_type::iterator compact_edgeR;
  typename FreqMapType::local_container_type::const_iterator fre;

	for (auto x : chain_reps) {
		// first get the kmer strings

		if ((::std::get<2>(x.second) == 0) && (::std::get<3>(x.second) == 0)) {
			L = x.first;
			R = x.first;
		} else if (::std::get<2>(x.second) == 0) {
			L = x.first;
			//        		  R = std::get<1>(x.second).reverse_complement();  // opposite strand
			R = std::get<1>(x.second); // we now want same strand as L.

		} else if (::std::get<3>(x.second) == 0) {
			L = x.first.reverse_complement();
			//        		  R = std::get<0>(x.second);  // opposite strand
			R = std::get<0>(x.second).reverse_complement(); // we now want same strand as L.
		} else {
			std::cout << "ERROR" << std::endl;
			continue;
		}

		edge_freq_type ef;
		::std::get<0>(ef) = L;
		::std::get<1>(ef) = R;


		// next print the left and right edges.
		cL = canonical(L);
		compact_edgeL = idx2.get_map().get_local_container().find(cL);
		assert(compact_edgeL != idx2.get_map().get_local_container().end());
		if (cL == L) {  // already canonical.  can use in edge directly.
			// get in edges of L
			std::get<0>(std::get<2>(ef)) = compact_edgeL->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
			std::get<1>(std::get<2>(ef)) = compact_edgeL->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
			std::get<2>(std::get<2>(ef)) = compact_edgeL->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
			std::get<3>(std::get<2>(ef)) = compact_edgeL->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
		} else {   // not canonical
			// get out edges of L, then complement each.  (reverse order)
			std::get<0>(std::get<2>(ef)) = compact_edgeL->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
			std::get<1>(std::get<2>(ef)) = compact_edgeL->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
			std::get<2>(std::get<2>(ef)) = compact_edgeL->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
			std::get<3>(std::get<2>(ef)) = compact_edgeL->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
		}

		cR = canonical(R);
		compact_edgeR = R_freq_map.find(cR);  // previously retrieved from remote.
		assert(compact_edgeR != R_freq_map.end());
		// we now assume R is on same strand as L
		if (cR == R) {  // already canonical
			// get out edges of R
			std::get<0>(std::get<3>(ef)) = compact_edgeR->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
			std::get<1>(std::get<3>(ef)) = compact_edgeR->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
			std::get<2>(std::get<3>(ef)) = compact_edgeR->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
			std::get<3>(std::get<3>(ef)) = compact_edgeR->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
		} else {   // not canonical
			// get in edges of R, then complement each.  (reverse order)
			std::get<0>(std::get<3>(ef)) = compact_edgeR->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
			std::get<1>(std::get<3>(ef)) = compact_edgeR->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
			std::get<2>(std::get<3>(ef)) = compact_edgeR->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
			std::get<3>(std::get<3>(ef)) = compact_edgeR->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
		}

		fre = freq_map.get_local_container().find(cL);
		assert(fre != freq_map.get_local_container().end());
		std::get<0>(std::get<4>(ef)) = (static_cast<float>(std::get<1>(fre->second)) /  static_cast<float>(std::get<0>(fre->second)));
		std::get<1>(std::get<4>(ef)) = std::get<2>(fre->second);
		std::get<2>(std::get<4>(ef)) = std::get<3>(fre->second);

		edge_freqs.emplace_back(ef);
	}
	BL_BENCH_COLLECTIVE_END(print_chain_freq, "gather_edge_freqs", edge_freqs.size(), comm);

	// psort
	int has_data = (edge_freqs.size() == 0) ? 0 : 1;
	int all_has_data = mxx::allreduce(has_data, comm);
	if (all_has_data > 0) {
		// global sort
		BL_BENCH_START(print_chain_freq);
		mxx::comm subcomm = (all_has_data == comm.size()) ? comm.copy() : comm.split(has_data);
		if (has_data == 1) {
			mxx::sort(edge_freqs.begin(), edge_freqs.end(), [](edge_freq_type const & x, edge_freq_type const & y){
				return std::get<0>(x) < std::get<0>(y);
			}, subcomm);
		}
		BL_BENCH_COLLECTIVE_END(print_chain_freq, "psort edge_freqs", edge_freqs.size(), comm);   // this is for constructing the chains
	}

	if (comm.rank() == 0) printf("PRINT CHAIN EDGE FREQS\n");

	// print
	BL_BENCH_START(print_chain_freq);
	{
		std::stringstream ss;
		ss.clear();
		for (auto x : edge_freqs) {
			ss << bliss::utils::KmerUtils::toASCIIString(std::get<0>(x)) << "\t" <<
					bliss::utils::KmerUtils::toASCIIString(std::get<1>(x)) << "\t" <<
					std::get<0>(std::get<2>(x)) << "\t" <<
					std::get<1>(std::get<2>(x)) << "\t" <<
					std::get<2>(std::get<2>(x)) << "\t" <<
					std::get<3>(std::get<2>(x)) << "\t" <<
					std::get<0>(std::get<3>(x)) << "\t" <<
					std::get<1>(std::get<3>(x)) << "\t" <<
					std::get<2>(std::get<3>(x)) << "\t" <<
					std::get<3>(std::get<3>(x)) << "\t" <<
					std::get<0>(std::get<4>(x)) << "\t" <<
					std::get<1>(std::get<4>(x)) << "\t" <<
					std::get<2>(std::get<4>(x)) << "\t" << std::endl;

		}
		write_mpiio(filename, ss.str().c_str(), ss.str().length(), comm);

	}
	BL_BENCH_COLLECTIVE_END(print_chain_freq, "print_edge_freqs", edge_freqs.size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(print_chain_freq, "print_chain_freqs", comm);

}




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
		TCLAP::ValueArg<CountType> lowerThreshArg("L", "lower_thresh", "Lower Threshold for Kmer and Edge frequency", false, 0, "uint16", cmd);
		TCLAP::ValueArg<CountType> upperThreshArg("U", "upper_thresh", "Upper Threshold for Kmer and Edge frequency", false,
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
    lower = lowerThreshArg.getValue();
    upper = upperThreshArg.getValue();
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

  BL_BENCH_INIT(app);

  // ================  read and get file
  BL_BENCH_START(app);
  ::std::vector<::bliss::io::file_data> file_data = open_files(filenames, comm, mpiio);
  BL_BENCH_COLLECTIVE_END(app, "read_file", file_data.size(), comm);
  // == DONE == reading

	::std::vector<::std::vector<::bliss::debruijn::biedge::compact_simple_biedge> > selected_edges;
	ChainMapType chainmap(comm);
  DBGType idx(comm);

	// =================  make compacted simple DBG, so that we can get chain and branch kmers.
	{
		BL_BENCH_START(app);

		// first read the input



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



		BL_BENCH_COLLECTIVE_END(app, "construct", idx.local_size(), comm);

		if (!benchmark) {
			BL_BENCH_START(app);
			print_edge_histogram(idx, comm);
			BL_BENCH_COLLECTIVE_END(app, "histo", idx.local_size(), comm);
			// == DONE == make compacted simple DBG

			check_index(idx, comm);
			printf("rank %d finished checking index\n", comm.rank());
		}
		// TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet
		//   this is done via read filtering/splitting.

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

	if (!benchmark) {
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


			// release chainmap
			BL_BENCH_START(app);
			chainmap.get_map().get_local_container().reset();
			BL_BENCH_COLLECTIVE_END(app, "chainmap_reset", chainmap.local_size(), comm);


			// do this first because we need the order of compacted chain to be same as hashed distribution.
			{
				// compute count index
				BL_BENCH_START(app);
				CountIndexType count_idx(comm);
				count_kmers(file_data, selected_edges, thresholding, count_idx, comm);
				BL_BENCH_COLLECTIVE_END(app, "count_kmers", count_idx.local_size(), comm);

				// compute freq map
				BL_BENCH_START(app);
				compute_freq_map(compacted_chain, count_idx, freq_map, comm);
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
	}

	BL_BENCH_REPORT_MPI_NAMED(app, "app", comm);


	// mpi cleanup is automatic

	return 0;

}

