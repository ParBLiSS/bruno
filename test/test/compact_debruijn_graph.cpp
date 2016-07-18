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
 * test_de_brujin_graph_construction.cpp
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

#include "utils/logging.h"

#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"

#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"

#include "iterators/transform_iterator.hpp"

#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "index/quality_score_iterator.hpp"

#include "index/kmer_index.hpp"
#include "index/kmer_hash.hpp"

#include "debruijn/debruijn_common.hpp"
#include "debruijn/edge_iterator.hpp"

#include "debruijn/debruijn_graph_node.hpp"
#include "debruijn/debruijn_graph_map.hpp"
#include "debruijn/debruijn_graph_loader.hpp"
#include "debruijn/debruijn_graph_filters.hpp"
#include "debruijn/debruijn_graph_operations.hpp"

#include "debruijn/debruijn_chain_filters.hpp"
#include "debruijn/debruijn_chain_node.hpp"
#include "debruijn/debruijn_chain_operations.hpp"

#include "debruijn/debruijn_stats.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"
#include "utils/tclap_utils.hpp"

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
using KmerType = bliss::common::Kmer<pK, Alphabet, WordType>;
#else
using KmerType = bliss::common::Kmer<31, Alphabet, WordType>;
#endif

//============== index input file format
#if (pPARSER == FASTA)
#define FileParser ::bliss::io::FASTAParser
#elif (pPARSER == FASTQ)
#define FileParser ::bliss::io::FASTQParser
#endif


using KmerType = bliss::common::Kmer<31, Alphabet, WordType>;
using EdgeEncoding = Alphabet;

using FileReaderType = ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, FileParser >;

using DBGNodeParser = bliss::debruijn::debruijn_graph_parser<KmerType>;

using DBGMapType = ::bliss::debruijn::graph::simple_hash_compact_debruijn_graph_map<KmerType>;
using DBGType = ::bliss::index::kmer::Index<DBGMapType, DBGNodeParser>;

using CountType = uint32_t;
using CountDBGMapType = ::bliss::debruijn::graph::count_hash_compact_debruijn_graph_map<KmerType, CountType>;
using CountDBGType = ::bliss::index::kmer::Index<CountDBGMapType, DBGNodeParser>;

using ChainNodeType = ::bliss::debruijn::simple_biedge<KmerType>;
//template <typename K>
//using ChainMapParams = ::bliss::index::kmer::CanonicalDebuijnHashMapParams<K>;
using ChainMapType = ::dsc::densehash_map<KmerType, ChainNodeType,
    ::bliss::debruijn::CanonicalDeBruijnHashMapParams,
     ::bliss::kmer::hash::sparsehash::special_keys<KmerType> >;

template <typename Key>
using FreqMapParams = ::bliss::index::kmer::CanonicalHashMapParams<Key>;


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

	std::cout << "rank " << comm.rank() << " mpiio write offset is " << global_offset << " len " << len << " iterations " << iterations << ::std::endl;

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

		res = MPI_File_write_at_all( fh, global_offset, data, curr_step, MPI_BYTE, &stat);

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
    TCLAP::CmdLine cmd("Parallel de bruijn graph compaction", ' ', "0.1");

    // MPI friendly commandline output.
    ::bliss::utils::tclap::MPIOutput cmd_output(comm);
    cmd.setOutput(&cmd_output);

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<std::string> outputArg("O", "output_prefix", "Prefix for output files, including directory", false, "", "string", cmd);

//    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);
    TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names", false, "string", cmd);

    // Parse the argv array.
    cmd.parse( argc, argv );

    filenames = fileArg.getValue();
    out_prefix = outputArg.getValue();

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(-1);
  }

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
//  compacted_chain_str_filename.append("_chain.");
//  compacted_chain_str_filename.append(std::to_string(comm.rank()));
//  compacted_chain_str_filename.append(".fasta");
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
//  compacted_chain_kmers_filename.append("_chain.");
//  compacted_chain_kmers_filename.append(std::to_string(comm.rank()));
//  compacted_chain_kmers_filename.append(".components");
  compacted_chain_kmers_filename.append("_chain.components");

  // filename for compacted chain end points, in format : < K_l, K_r, chain_id, f_l_A, f_l_C, f_l_G, f_l_T, f_r_A, f_r_C, f_r_G, f_r_T, f >
  // K_l and K_r are canonical.  as such, we need chain_id to link back to component k-mers and the compated strings.
  // K_l is the smaller of K_l and rev_comp(original K_r).
  // K_l and original K_r are on same strand.
  // f is frequency - minimum or average of all intervening nodes?  - a reduction similar to the one used to construct the compacted chain string is needed.
  // f_l_X and f_r_X are in and out edges of K_l and K_r, respectively.
  std::string compacted_chain_ends_filename(out_prefix);
//  compacted_chain_ends_filename.append("_chain.");
//  compacted_chain_ends_filename.append(std::to_string(comm.rank()));
//  compacted_chain_ends_filename.append(".edges");
  compacted_chain_ends_filename.append("_chain.edges");

  // filename for junctions in format <K, f_l_A, f_l_C, f_l_G, f_l_T, f_r_A, f_r_C, f_r_G, f_r_T, f>
  // K is canonical.
  // f is frequency
  // f_l_X and f_r_X are in and out edges of K_l and K_r, respectively.
  // this is a dump of the dbg junctional nodes (filtered) to disk.
  std::string branch_filename(out_prefix);
//  branch_filename.append("_branch.");
//  branch_filename.append(std::to_string(comm.rank()));
//  branch_filename.append(".edges");
  branch_filename.append("_branch.edges");


  // ================  read and get file

  DBGType idx(comm);
  ChainMapType chainmap(comm);



  std::vector<KmerType> query;
  std::vector<KmerType> neighbors;
  neighbors.reserve(4);

//	::std::vector<typename DBGNodeParser::value_type> temp;
  ::std::vector<typename ::bliss::io::file_data> file_data;

//  KmerType testKmer(std::string("CAAGATGGGTGGAATGGCCAGTTAACCACTG"));

  // TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet
  BL_BENCH_INIT(test);

  {
	BL_BENCH_RESET(test);

    BL_BENCH_START(test);
    size_t total = 0;
    for (auto fn : filenames) {
		if (comm.rank() == 0) printf("reading %s via posix\n", fn.c_str());

		file_data.push_back(idx.open_file<FileReaderType>(fn, comm));
		total += file_data.back().getRange().size();
//		idx.read_file_posix<FileParser, DBGNodeParser>(fn, temp1, comm);
    }
    BL_BENCH_COLLECTIVE_END(test, "read", total, comm);

    //		  for (auto t : temp) {
    //		    std::cout << "input kmer " << bliss::utils::KmerUtils::toASCIIString(t.first) << " edges " << t.second << std::endl;
    //		  }
    //

    // all possible k-mers from input should already be present, so no need to test.

    // ============== insert

    total = mxx::allreduce(total, comm);
    if (comm.rank() == 0) printf("total size read is %lu\n", total);

    BL_BENCH_REPORT_MPI_NAMED(test, "open", comm);
  }



  {
	  BL_BENCH_RESET(test);

	  BL_BENCH_START(test);
	  {
		::std::vector<typename DBGNodeParser::value_type> temp;
		for (auto x : file_data) {
			temp.clear();
			idx.parse_file_data<FileParser, DBGNodeParser>(x, temp, comm);
			idx.insert(temp);
		}
	  }
	  BL_BENCH_COLLECTIVE_END(test, "parse and insert", idx.local_size(), comm);

    size_t total = idx.size();
    if (comm.rank() == 0) printf("total size after insert/rehash is %lu\n", total);


    // get all possible edges.
    //============== testing to ensure that all the possible edges are present.
    //    	std::cout << "COMPRESSED MAP" << std::endl;
    //	    auto cc = idx.get_map().get_local_container();
    //	    for (auto it = cc.begin(); it != cc.end(); ++it) {
    //
    //	    	std::cout << "kmer: " << it->first << " edge " << it->second << std::endl;
    //
    //	    	neighbors.clear();
    //	    	it->second.get_out_neighbors(it->first, neighbors);
    //	        query.insert(query.end(), neighbors.begin(), neighbors.end());
    //
    //	        neighbors.clear();
    //	        it->second.get_in_neighbors(it->first, neighbors);
    //	        query.insert(query.end(), neighbors.begin(), neighbors.end());
    //	    }

    //	    // =============== check to see if index is superset of query.  (count should have every result entry showing 1.)
    //	    {
    //	      auto lquery = query;
    //	      BL_BENCH_START(test);
    //	      auto counts = idx.count(lquery);
    //	      BL_BENCH_COLLECTIVE_END(test, "count", counts.size(), comm);
    //
    //	      auto absent_end = std::partition(counts.begin(), counts.end(), [](std::pair<KmerType, size_t> const & x){
    //	    	  return x.second == 0;
    //	      });
    //	      printf(" total query = %lu, unique query = %lu, unique absent = %lu\n", query.size(), counts.size(), std::distance(counts.begin(), absent_end));
    //
    //	      for (auto it = counts.begin(); it != absent_end; ++it) {
    //	    	  std::cout << "  " << it->first << std::endl;
    // 	      }
    //	    }
    //
    //	    comm.barrier();
    //
    //	    // =============== check to see if query is superset of index.  (erase should result in empty)
    //	    {
    //	    	auto lquery = query;
    //		      BL_BENCH_START(test);
    //		      size_t erased = idx.get_map().erase(lquery);
    //		      BL_BENCH_COLLECTIVE_END(test, "erase", erased, comm);
    //
    //		      printf(" total query = %lu, erased = %lu, remaining = %lu\n", query.size(), erased, idx.local_size());
    //
    //	    }
    //
    //	      comm.barrier();
    //	      idx.get_map().erase(::fsc::TruePredicate());
    //	      comm.barrier();
    //		  BL_BENCH_START(test);
    //		  idx.insert(temp);
    //		  BL_BENCH_COLLECTIVE_END(test, "reinsert", idx.local_size(), comm);


    {
      // get histogram for the edge types in debruijn graph

      ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > nodes;

      BL_BENCH_START(test);
      // get the content as an array
      idx.get_map().to_vector(nodes);
      BL_BENCH_COLLECTIVE_END(test, "to_vector", nodes.size(), comm);

      BL_BENCH_START(test);
      // then compute histogram
      ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(nodes, comm);
      BL_BENCH_COLLECTIVE_END(test, "histogram", nodes.size(), comm);
    }


    BL_BENCH_REPORT_MPI_NAMED(test, "insert", comm);

  }

  {
	BL_BENCH_RESET(test);

    //	      {
    //	        auto lquery = query;
    //
    //	      BL_BENCH_START(test);
    //	      auto found = idx.find_overlap(lquery);
    //	      BL_BENCH_COLLECTIVE_END(test, "find_overlap", found.size(), comm);
    //	      }
    //	    // separate test because of it being potentially very slow depending on imbalance.



    {
      ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > chain_nodes;

      BL_BENCH_START(test);
      // find chain nodes
      chain_nodes = idx.find_if(::bliss::debruijn::filter::graph::IsChainNode());  // not isolated and not branching
      BL_BENCH_COLLECTIVE_END(test, "get_chains", chain_nodes.size(), comm);

//      BL_BENCH_START(test);
//      // then compute histogram
//      ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(chain_nodes, comm);
//      BL_BENCH_COLLECTIVE_END(test, "chain_histogram", chain_nodes.size(), comm);


      // insert into local container inside chainmap.
      BL_BENCH_START(test);

      // initialize the chain map.  note that key kmers are canonical, same as in DBG.
      // note also that edge k-mers have same orientation as key kmers, but may not be canonical.
      for (auto t : chain_nodes) {
        ChainNodeType node(KmerType(), KmerType(), 0, 0);   // default node

        // get the in neighbor
        neighbors.clear();
        t.second.get_in_neighbors(t.first, neighbors);
        assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
        if (neighbors.size() == 1) {
          std::get<0>(node) = neighbors[0];
          std::get<2>(node) = 1;
        }

        // get the out neighbor
        neighbors.clear();
        t.second.get_out_neighbors(t.first, neighbors);
        assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
        if (neighbors.size() == 1) {
          std::get<1>(node) = neighbors[0];
          std::get<3>(node) = 1;
        }

        // chainmap uses the same distribution hash function and transform, so can insert locally.
        chainmap.get_local_container().insert(::std::make_pair(::std::move(t.first), ::std::move(node)));
        //            std::cout << "BIEDGE\tin dist\t" << std::get<2>(node) << " kmer: " << std::get<0>(node) << std::endl;
        //            std::cout << "\tKmer bi edge:\t" << t.first << std::endl;
        //            std::cout << "\tout dist\t" << std::get<3>(node) << " kmer: " << std::get<1>(node) << std::endl;

      }
      BL_BENCH_COLLECTIVE_END(test, "insert in chainmap.", chainmap.local_size(), comm);

      //========= report.
//      auto result = chainmap.find(::bliss::debruijn::filter::chain::IsTerminus());
//      auto result2 = chainmap.find(::bliss::debruijn::filter::chain::IsIsolated());
//      printf("chain map contains %lu chained termini and  %lu isolated\n", result.size(), result2.size());

    }


    {
      //=== find branching nodes. local computation.
      BL_BENCH_START(test);
      auto nodes = idx.find_if(::bliss::debruijn::filter::graph::IsBranchPoint());
      BL_BENCH_COLLECTIVE_END(test, "get_branches", nodes.size(), comm);

//      BL_BENCH_START(test);
//      // then compute histogram
//      ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(nodes, comm);
//      BL_BENCH_COLLECTIVE_END(test, "branch_histogram", nodes.size(), comm);

      //          //=== get the neighbors of the branch points.  for information only.
      //          BL_BENCH_START(test);
      //          std::vector<KmerType> all_neighbors2;
      //          all_neighbors2.reserve(nodes.size() * 4);
      //
      //          for (auto t : nodes) {
      //            neighbors.clear();
      //            t.second.get_out_neighbors(t.first, neighbors);
      //            all_neighbors2.insert(all_neighbors2.end(), neighbors.begin(), neighbors.end());
      //
      //            neighbors.clear();
      //            t.second.get_in_neighbors(t.first, neighbors);
      //            all_neighbors2.insert(all_neighbors2.end(), neighbors.begin(), neighbors.end());
      //          }
      //          BL_BENCH_COLLECTIVE_END(test, "branch_neighbors", all_neighbors2.size(), comm);
      //
      //          // now check to see which are chain nodes.  these are chain nodes adjacent to branch points.
      //          // include chain termini that are adjacent to branch points, so we can mark them in the chainmap.
      //          BL_BENCH_START(test);
      //          auto found = idx.find_if_overlap(all_neighbors2, ::bliss::debruijn::filter::graph::IsChainNode());
      //          BL_BENCH_COLLECTIVE_END(test, "terminal_neighbors", found.size(), comm);
      //
      //          BL_BENCH_START(test);
      //          // then compute histogram
      //          ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(found, comm);
      //          BL_BENCH_COLLECTIVE_END(test, "termini_histogram", found.size(), comm);
      //
      //
      //          // now check to see which are chain nodes.  these are chain nodes adjacent to branch points.
      //          // include chain termini that are adjacent to branch points, so we can mark them in the chainmap.
      //          BL_BENCH_START(test);
      //          auto found2 = idx.find_if_overlap(all_neighbors2, ::bliss::debruijn::filter::graph::IsTerminus());
      //          BL_BENCH_COLLECTIVE_END(test, "stump_termini", found2.size(), comm);
      //
      //          BL_BENCH_START(test);
      //          // then compute histogram
      //          ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(found2, comm);
      //          BL_BENCH_COLLECTIVE_END(test, "stump_histogram", found2.size(), comm);


      //=== mark neighbors of branch points.
      BL_BENCH_START(test);
      std::vector<std::pair<KmerType, bliss::debruijn::operation::chain::terminus_update_md<KmerType> > > all_neighbors;
      all_neighbors.reserve(nodes.size() * 4);

      for (auto t : nodes) {
        neighbors.clear();
        t.second.get_out_neighbors(t.first, neighbors);
        for (auto n : neighbors) {
          // insert as is.  let lex_less handle flipping it.
          all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::debruijn::operation::IN));
        }

        neighbors.clear();
        t.second.get_in_neighbors(t.first, neighbors);
        for (auto n : neighbors) {
          // insert as is.  let lex_less handle flipping it.
          all_neighbors.emplace_back(n, bliss::debruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::debruijn::operation::OUT));
        }
      }
      BL_BENCH_COLLECTIVE_END(test, "branch_neighbors_2", all_neighbors.size(), comm);

      // now check to see which are chain nodes.  these are chain nodes adjacent to branch points.
      // include chain termini that are adjacent to branch points, so we can mark them in the chainmap.
      BL_BENCH_START(test);
      ::bliss::debruijn::operation::chain::terminus_update<KmerType> updater;
      size_t count = chainmap.update(all_neighbors, false, updater );
      BL_BENCH_COLLECTIVE_END(test, "update_termini", count, comm);

      //          for (auto it = chainmap.get_local_container().begin(); it != chainmap.get_local_container().end(); ++it) {
      //              std::cout << "TERMED\tin dist\t" << std::get<2>(it->second) << " kmer: " << std::get<0>(it->second) << std::endl;
      //              std::cout << "\tKmer bi edge:\t" << it->first << std::endl;
      //              std::cout << "\tout dist\t" << std::get<3>(it->second) << " kmer: " << std::get<1>(it->second) << std::endl;
      //          }

      //          //========= split singleton entries from chainmap.
      //          {
      //            BL_BENCH_START(test);
      //            auto single_chains = chainmap.find(::bliss::debruijn::filter::chain::IsIsolated());
      //            BL_BENCH_COLLECTIVE_END(test, "singletons", single_chains.size(), comm);
      //            printf("found %lu singletons in %lu chainmap\n", single_chains.size(), chainmap.local_size());
      //          }
      //
      //          {
      //            BL_BENCH_START(test);
      //            auto single_chains = chainmap.find(::bliss::debruijn::filter::chain::IsTerminus());
      //            BL_BENCH_COLLECTIVE_END(test, "termini", single_chains.size(), comm);
      //            printf("found %lu termini in %lu chainmap\n", single_chains.size(), chainmap.local_size());
      //          }
    }
    BL_BENCH_REPORT_MPI_NAMED(test, "chain_map", comm);

    BL_BENCH_RESET(test);

    //        {
    //          // NOW: do the list ranking
    //
    //          // search unfinished
    //          BL_BENCH_START(test);
    //          auto unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
    //          BL_BENCH_COLLECTIVE_END(test, "unfinished", unfinished.size(), comm);
    //
    //          // get global unfinished count
    //          bool all_compacted = (unfinished.size() == 0);
    //          all_compacted = ::mxx::all_of(all_compacted, comm);
    //
    //
    //          // while have unfinished,  run
    //          BL_BENCH_START(test);
    //          size_t iterations = 0;
    //
    //          std::vector<KmerType> qq;
    //          qq.reserve(unfinished.size() * 2);
    //
    //
    //          while (!all_compacted) {
    //
    //        	  qq.clear();
    //            // get left and right edges, generate updates
    //
    ////            std::cout << "iteration " << iterations << " kmer: " << unfinished[0].first <<
    ////                " in: " << std::get<0>(unfinished[0].second) << " out: " << std::get<1>(unfinished[0].second) <<
    ////                " in dist " << std::get<2>(unfinished[0].second) << " out dist " << std::get<3>(unfinished[0].second) << std::endl;
    //
    //        	  // get the end points.
    //            for (auto t : unfinished) {
    //                md = t.second;
    //            	std::cout << "it " << iterations;
    //            	std::cout << "\tin dist " << std::get<2>(md) << " kmer: " << std::get<0>(md) << std::endl;
    //            	std::cout << "\tkmer: " << t.first << std::endl;
    //            	std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << std::get<1>(md) << std::endl;
    //
    //            	// add the in edge target if it needs to be updated.
    //            	if (std::get<2>(md) > 0) {
    //            		qq.emplace_back(std::get<0>(md));
    //            	}
    //            	if (std::get<3>(md) > 0) {
    //            		qq.emplace_back(std::get<1>(md));
    //            	}
    //
    //            }
    //            printf("iter %ld query size = %ld\n", iterations, qq.size());
    //
    //            comm.barrier();
    //
    //            // now query.
    //            auto results = chainmap.find(qq);
    //            printf("iter %ld query results size = %ld\n", iterations, results.size());
    //
    //
    //            // now do local updates.
    //            // go throught the unfinished list again
    //            // need info about where the query came from (which k-mer?)
    //            // else use a local map to index the results.
    //
    //
    //
    //
    //            // search unfinished.
    //            unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
    //
    //
    //            // get global unfinished count
    //            all_compacted = (unfinished.size() == 0);
    //            all_compacted = ::mxx::all_of(all_compacted, comm);
    //
    //            ++iterations;
    //
    //          }
    //          BL_BENCH_COLLECTIVE_END(test, "compact", iterations, comm);
    //
    //
    //          // ========= check compacted chains.
    //
    //
    //
    //
    //        }



    size_t iterations = 0;
    size_t cycle_nodes = 0;

    {
      // NOW: do the list ranking

      // search unfinished
      BL_BENCH_START(test);
      auto unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
      BL_BENCH_COLLECTIVE_END(test, "unfinished", unfinished.size(), comm);

      // get global unfinished count
      bool all_compacted = (unfinished.size() == 0);
      all_compacted = ::mxx::all_of(all_compacted, comm);


      // while not same, run
      BL_BENCH_START(test);

      std::vector<std::pair<KmerType, bliss::debruijn::operation::chain::chain_update_md<KmerType> > > updates;
      updates.reserve(unfinished.size() * 2);

      int dist = 0;
      KmerType ll, rr;
      bliss::debruijn::simple_biedge<KmerType> md;

      //          size_t last_updated = 0;

      while (!all_compacted) {

        updates.clear();

        // get left and right edges, generate updates

        //            std::cout << "iteration " << iterations << " kmer: " << unfinished[0].first <<
        //                " in: " << std::get<0>(unfinished[0].second) << " out: " << std::get<1>(unfinished[0].second) <<
        //                " in dist " << std::get<2>(unfinished[0].second) << " out dist " << std::get<3>(unfinished[0].second) << std::endl;


        for (auto t : unfinished) {
          md = t.second;



          // each is a pair with kmer, <in kmer, out kmer, in dist, out dist>
          // constructing 2 edges <in, out> and <out, in>  distance is sum of the 2.
          // indication of whether edge destination is a terminus depend only on the sign of distance to that node.

          dist = abs(std::get<2>(md)) + abs(std::get<3>(md));   // this double the distance...

          // and below pointer jumps.

          // construct forward edge, from ll to rr, only if current node is not a terminus for the "in" side
          if (std::get<2>(md) != 0)  {
            // send rr to ll.  also let ll know if rr is a terminus.  orientation is OUT
            updates.emplace_back(std::get<0>(md),
                                 bliss::debruijn::operation::chain::chain_update_md<KmerType>((std::get<3>(md) == 0) ? t.first : std::get<1>(md),
                                     ((std::get<3>(md) > 0) ? dist : -dist),
                                     bliss::debruijn::operation::OUT));
            // if out dist is 0, then this node is a terminal node.  sent self as target.  else use right kmer.
            // if out dist is negative, then out kmer (rr) points to a terminus, including self (dist = 0), set update distance to negative to indicate so.

          }  // else case is same as below

          // construct backward edge, from out to in, only if current node is not a terminus for the "out" side
          if (std::get<3>(md) != 0) {
            // send ll to rr.  also let rr know if ll is a terminus.  orientation is IN
            updates.emplace_back(std::get<1>(md),
                                 bliss::debruijn::operation::chain::chain_update_md<KmerType>((std::get<2>(md) == 0) ? t.first : std::get<0>(md),
                                     ((std::get<2>(md) > 0) ? dist : -dist),
                                     bliss::debruijn::operation::IN));
            // if target is a terminus, then set self as target.  else use left kmer
            // if target points to a terminus, including self (dist = 0), set update distance to negative to indicate so.
          }  // else case is same as above



          // if ((std::get<2>(md) == 0) && (std::get<3>(md) == 0)) continue;  // singleton.   next.
        }

        ////            if (last_updated == updates.size()) {
        //                for (auto t : unfinished) {
        //                	md = t.second;
        //
        //                	KmerType tt = testKmer.reverse_complement();
        //                	if ( ((t.first == testKmer) || (t.first == tt) ) ||
        //                		((std::get<0>(md) == testKmer) || (std::get<0>(md) == tt)) ||
        //                			((std::get<1>(md) == testKmer) || (std::get<1>(md) == tt))) {
        //
        //						std::cout << "it " << iterations;
        //						std::cout << "\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) <<
        //								" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) << std::endl;
        //						std::cout << "\tkmer: " << bliss::utils::KmerUtils::toASCIIString(t.first) << " rc: " << bliss::utils::KmerUtils::toASCIIString(t.first.reverse_complement()) << std::endl;
        //						std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) <<
        //								" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
        //                	}
        //                }
        ////            }


        comm.barrier();

        // now perform update
        ::bliss::debruijn::operation::chain::chain_update<KmerType> chain_updater;
        size_t count = chainmap.update(updates, false, chain_updater );

        //            last_updated = count;

        // search unfinished.
        unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());

        //            {
        //              auto t = unfinished[0];
        //              auto md = t.second;
        //  			std::cout << "it " << iterations;
        //  			std::cout << "\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) <<
        //  					" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) << std::endl;
        //  			std::cout << "\tkmer: " << bliss::utils::KmerUtils::toASCIIString(t.first) << " rc: " << bliss::utils::KmerUtils::toASCIIString(t.first.reverse_complement()) << std::endl;
        //  			std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) <<
        //  					" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
        //            }




        // at this point, the new distances in lists are 2^(iterations + 1)
        ++iterations;

        cycle_nodes = std::count_if(unfinished.begin(), unfinished.end(),
                                    ::bliss::debruijn::filter::chain::IsCycleNode(iterations));

        // going over 30 makes the max_dist in IsCycleNode go to -1, then it is no longer valid as distances are int.  stop at 30
        if (iterations >= 30) {

        	// locally find the one that are not cycle nodes and print them.
        	::bliss::debruijn::filter::chain::IsCycleNode check_cycle(iterations);
        	for (auto t : unfinished) {
        		if (check_cycle(t)) continue;

        		auto md = t.second;

        		std::cout << "rank " << comm.rank() << " max iter " << iterations <<
        		"\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) <<
							"\tkmer: " << bliss::utils::KmerUtils::toASCIIString(t.first) << " rc: " << bliss::utils::KmerUtils::toASCIIString(t.first.reverse_complement()) <<
							"\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) <<
									" rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
        	}

//        	printf("rank %d max iter %lu updated %lu, unfinished %lu cycle nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_nodes);
        	all_compacted = true;
        	continue;
        }

        // get global unfinished count

        all_compacted = (count == 0) || (cycle_nodes == unfinished.size());
        if (!all_compacted) printf("rank %d iter %lu updated %lu, unfinished %lu cycle nodes %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_nodes);
		all_compacted = ::mxx::all_of(all_compacted, comm);

      }
      BL_BENCH_COLLECTIVE_END(test, "compact", cycle_nodes, comm);

    }

    BL_BENCH_REPORT_MPI_NAMED(test, "list_rank", comm);

    BL_BENCH_RESET(test);

    // ==== finally update all nodes that are not yet set to pointing to termini.
    // should just be the ones that are pointing to termini but values are not yet marked as negative.
    // previously has been sending out remote updates based on local information.
    // due to pointer doubling, and local updating to the farthest remote jump,  some remote
    // sources of updates may not be updated again, near the ends, even when they are already pointing to terminal
    // these are resolved using a single query.  note that previous code terminates when nothing is updated or when only cycles remain.
    {

      // search unfinished
      BL_BENCH_START(test);
      auto unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
      BL_BENCH_COLLECTIVE_END(test, "unfinished2", unfinished.size(), comm);

      BL_BENCH_START(test);

      // get global unfinished count
      bool all_compacted = (unfinished.size() == 0);
      all_compacted = ::mxx::all_of(all_compacted, comm);


      // while have unfinished,  run.  qq contains kmers not necessarily canonical

      std::vector<KmerType> qq;
      qq.reserve(unfinished.size() * 2);

      if (!all_compacted) {

        qq.clear();
        // get left and right edges, generate updates

        //            std::cout << "iteration " << iterations << " kmer: " << unfinished[0].first <<
        //                " in: " << std::get<0>(unfinished[0].second) << " out: " << std::get<1>(unfinished[0].second) <<
        //                " in dist " << std::get<2>(unfinished[0].second) << " out dist " << std::get<3>(unfinished[0].second) << std::endl;

        // get the end points.
        bliss::debruijn::simple_biedge<KmerType> md;

        for (auto t : unfinished) {
          md = t.second;
          //				std::cout << "\tin dist " << std::get<2>(md) << " kmer: " << std::get<0>(md) << std::endl;
          //				std::cout << "\tkmer: " << t.first << std::endl;
          //				std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << std::get<1>(md) << std::endl;

          // add the in edge target if it needs to be updated.
          if (std::get<2>(md) > 0) {
            qq.emplace_back(std::get<0>(md));
          }
          if (std::get<3>(md) > 0) {
            qq.emplace_back(std::get<1>(md));
          }

        }

        comm.barrier();

        // now query.
        auto results = chainmap.find(qq);


        // put query results in a map.  key is CANONICAL
        typename ChainMapType::local_container_type res_map(results.begin(), results.end());

        // now update.  for each unfinished
        KmerType kk;
        ::bliss::kmer::transform::lex_less<KmerType> lexless;
        for (auto t : unfinished) {
          // get left (in) kmer.
          kk = std::get<0>(t.second);

          // lookup left by canonical
          auto it = res_map.find(lexless(kk));

          // if left is at terminus, and we haven't updated it,
          if (((std::get<2>(it->second) == 0) || (std::get<3>(it->second) == 0)) &&
              (std::get<2>(t.second) > 0)) {
            // update the current left
            std::get<2>((*(chainmap.get_local_container().find(t.first))).second) = -(std::get<2>(t.second));
          }

          // get right (out) kmer
          kk = std::get<1>(t.second);

          // lookup left by canonical
          it = res_map.find(lexless(kk));

          // if left is at terminus, and we haven't updated it,
          if (((std::get<2>(it->second) == 0) || (std::get<3>(it->second) == 0)) &&
              (std::get<3>(t.second) > 0)) {
            // update the current left
            std::get<3>((*(chainmap.get_local_container().find(t.first))).second) = -(std::get<3>(t.second));
          }

        }

      }
      BL_BENCH_COLLECTIVE_END(test, "cleanup", unfinished.size(), comm);


      // search unfinished.count
      BL_BENCH_START(test);

      unfinished = chainmap.find(::bliss::debruijn::filter::chain::IsUncompactedNode(iterations));


      // get global unfinished count
      all_compacted = (unfinished.size() == 0);
      all_compacted = ::mxx::all_of(all_compacted, comm);

      assert(all_compacted);
      BL_BENCH_COLLECTIVE_END(test, "unfinished", unfinished.size(), comm);

    }
    BL_BENCH_REPORT_MPI_NAMED(test, "finalize", comm);


    BL_BENCH_RESET(test);
    {
		// remove cycle nodes.  has to do after query, to ensure that exactly middle is being treated not as a cycle node.
		BL_BENCH_START(test);
		size_t cycle_node_count = chainmap.erase(::bliss::debruijn::filter::chain::IsCycleNode(iterations));
		BL_BENCH_COLLECTIVE_END(test, "erase cycle", cycle_node_count, comm);

		BL_BENCH_START(test);
		auto unfinished = chainmap.find(::bliss::debruijn::filter::chain::IsUncompactedNode(iterations));

		bool all_compacted = (unfinished.size() == 0);
		all_compacted = ::mxx::all_of(all_compacted, comm);

		assert(all_compacted);
		BL_BENCH_COLLECTIVE_END(test, "unfinished", unfinished.size(), comm);

	  }

  BL_BENCH_REPORT_MPI_NAMED(test, "rem_cycle", comm);


  //============= PRINTING....

    BL_BENCH_RESET(test);
    {
      // ========== construct edge count index
    	CountDBGType idx2(comm);
  	  BL_BENCH_START(test);
  	  {
		::std::vector<typename DBGNodeParser::value_type> temp;
		for (auto x : file_data) {
			temp.clear();
			idx2.parse_file_data<FileParser, DBGNodeParser>(x, temp, comm);
			idx2.insert(temp);
		}
		// idx2.insert(temp);
  	  }
      BL_BENCH_COLLECTIVE_END(test, "count insert", idx2.local_size(), comm);

      {
		  // then find branches.
		  BL_BENCH_START(test);
		  std::vector<typename CountDBGType::TupleType> branch_pts =
				  idx2.find_if(::bliss::debruijn::filter::graph::IsBranchPoint());
		  BL_BENCH_COLLECTIVE_END(test, "get_branches", branch_pts.size(), comm);

		  bool no_branches = mxx::all_of((branch_pts.size() == 0), comm);
		  if (!no_branches) {
			  // global sort
			  BL_BENCH_START(test);
			  mxx::sort(branch_pts.begin(), branch_pts.end(), [](typename CountDBGType::TupleType const & x,
					  typename CountDBGType::TupleType const & y){
				  return x.first < y.first;
			  }, comm);
			  BL_BENCH_COLLECTIVE_END(test, "psort branches", branch_pts.size(), comm);   // this is for ordered output.
		  }

		  // and print.
		  BL_BENCH_START(test);
		  {
			  std::stringstream ss;
			  ss.clear();
			  std::for_each(branch_pts.begin(), branch_pts.end(),
					  ::bliss::debruijn::operation::graph::print_graph_node<KmerType>(ss));
			  write_mpiio(branch_filename, ss.str().c_str(), ss.str().length(), comm);

//			  std::ofstream ofs_branch_nodes(branch_filename);
//			  ofs_branch_nodes << ss.str();
//			  ofs_branch_nodes.close();
		  }
		  BL_BENCH_COLLECTIVE_END(test, "print branches (4)", branch_pts.size(), comm);
      }


      // ========== construct new graph with compacted chains and junction nodes.
      {
          BL_BENCH_START(test);
          std::vector<::bliss::debruijn::chain::compacted_chain_node<KmerType> > result;
          result.reserve(chainmap.size());
          ::fsc::back_emplace_iterator<std::vector<::bliss::debruijn::chain::compacted_chain_node<KmerType> > > back_emplacer(result);

          //== first transform nodes so that we are pointing to canonical terminus k-mers.
          std::transform(chainmap.get_local_container().begin(), chainmap.get_local_container().end(), back_emplacer,
        		  ::bliss::debruijn::operation::chain::to_compacted_chain_node<KmerType>());
          BL_BENCH_COLLECTIVE_END(test, "transform chain", chainmap.local_size(), comm);

		  bool no_result = mxx::all_of((result.size() == 0), comm);
		  if (!no_result) {

			  // sort
			  BL_BENCH_START(test);
			  // first transform nodes so that we are pointing to canonical terminus k-mers.
			  mxx::sort(result.begin(), result.end(), ::bliss::debruijn::operation::chain::chain_rep_less<KmerType>());
			  // global sort?

			  BL_BENCH_COLLECTIVE_END(test, "psort lmer", result.size(), comm);   // this is for constructing the chains
		  }
		  // print out.
		  BL_BENCH_START(test);
		  {
			  std::stringstream ss;

			  std::for_each(result.begin(), result.end(), ::bliss::debruijn::operation::chain::print_chain<KmerType>(ss));
			  write_mpiio(compacted_chain_str_filename, ss.str().c_str(), ss.str().length(), comm);

//			  std::ofstream ofs_chain_str(compacted_chain_str_filename);
//			  ofs_chain_str << ss.str();
//			  ofs_chain_str.close();

		  }
	//      std::cout << ss.str() << std::endl;
		  BL_BENCH_COLLECTIVE_END(test, "print chains (3)", result.size(), comm);


		//===  print chain nodes (1)

		  no_result = mxx::all_of((result.size() == 0), comm);
		  if (!no_result) {

			  // sort
			  BL_BENCH_START(test);
			  // first transform nodes so that we are pointing to canonical terminus k-mers.
			  mxx::sort(result.begin(), result.end(), ::bliss::debruijn::operation::chain::chain_node_less<KmerType>(), comm);
			  // global sort?
			  BL_BENCH_COLLECTIVE_END(test, "psort kmer", result.size(), comm);  // this is for output ordering.
		  }

		  // print out.
		  BL_BENCH_START(test);
		  {
			  std::stringstream ss2;
			  std::for_each(result.begin(), result.end(), ::bliss::debruijn::operation::chain::print_chain_node<KmerType>(ss2));
			  write_mpiio(compacted_chain_kmers_filename, ss2.str().c_str(), ss2.str().length(), comm);
//
//			  std::ofstream ofs_chain_nodes(compacted_chain_kmers_filename);
//			  ofs_chain_nodes << ss2.str();
//			  ofs_chain_nodes.close();

		  }
	//      std::cout << ss.str() << std::endl;
		  BL_BENCH_COLLECTIVE_END(test, "print chains (1)", result.size(), comm);
      }


      // ===========  print compacted chain with left and right
      {
    	  // ==  first compute frequency summary, and store into a reduction map
    	  // allocate input
    	  using freq_type = std::pair<KmerType, std::tuple<CountType, size_t, CountType, CountType> >;
    	  std::vector< freq_type > freqs;

		  BL_BENCH_START(test);
    	  ::bliss::debruijn::operation::chain::to_compacted_chain_node<KmerType> get_chain_rep;

    	  // extract frequencies.
    	  for (auto x : chainmap.get_local_container()) {
    		  // compute the chain rep
    		  CountType c = idx2.get_map().get_local_container().find(x.first)->second.get_self_frequency();

    		  freqs.emplace_back(std::get<1>(get_chain_rep(x)), ::std::tuple<CountType, size_t, CountType, CountType>(1, c, c, c));
    	  }
		  BL_BENCH_COLLECTIVE_END(test, "get_freqs", freqs.size(), comm);

    	  // create a reduction map
		  BL_BENCH_START(test);
    	  using FreqMapType = ::dsc::reduction_densehash_map<KmerType, ::std::tuple<CountType, size_t, CountType, CountType>,
    			  FreqMapParams,
				   ::bliss::kmer::hash::sparsehash::special_keys<KmerType>,
					::bliss::debruijn::operation::chain::freq_summary<CountType> >;

    	  FreqMapType freq_map(comm);
    	  freq_map.insert(freqs);   // collective comm.
		  BL_BENCH_COLLECTIVE_END(test, "reduce_freq", freq_map.local_size(), comm);


    	  // freq_map, idx2, and chainmap now all have same distribution.
    	  // search in chainmap to find canonical termini.
          BL_BENCH_START(test);
          auto chain_rep = chainmap.find(::bliss::debruijn::filter::chain::IsCanonicalTerminus());
          BL_BENCH_COLLECTIVE_END(test, "chain rep", chain_rep.size(), comm);


          KmerType L, R, cL, cR;
          //========= get the R frequencies (remote) and insert into a local map
          using edge_freq_type = std::tuple<KmerType, KmerType, std::tuple<CountType, CountType, CountType, CountType>,
		  	  std::tuple<CountType, CountType, CountType, CountType>,
			  std::tuple<CountType, CountType, CountType> >;
          std::vector<edge_freq_type> edge_freqs;

          // get query vector
          BL_BENCH_START(test);
          typename CountDBGMapType::local_container_type R_freq_map;
          {
			  std::vector<KmerType> R_query;
			  for (auto x : chain_rep) {
				  if (::std::get<2>(x.second) == 0) {
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
			  auto R_results = idx2.find_overlap(R_query);
			  R_freq_map.insert(R_results.begin(), R_results.end());
//			  for (auto x : R_freq_map) {
//				  std::cout << "rank " << comm.rank() << " R result " << bliss::utils::KmerUtils::toASCIIString(x.first) << std::endl;
//			  }
          }
          BL_BENCH_COLLECTIVE_END(test, "local_R_freqs", R_freq_map.size(), comm);

          // convert to tuple.
          BL_BENCH_START(test);
          bliss::debruijn::lex_less<KmerType> canonical;
          for (auto x : chain_rep) {
        	  // first get the kmer strings

        	  if (::std::get<2>(x.second) == 0) {
        		  L = x.first;
        		  R = std::get<1>(x.second).reverse_complement();
        		  // left
//        		  ofs_chain_ends << bliss::utils::KmerUtils::toASCIIString(L) << "\t" <<
//        				  bliss::utils::KmerUtils::toASCIIString(R) << "\t";

        	  } else if (::std::get<3>(x.second) == 0) {
        		  L = x.first.reverse_complement();
        		  R = std::get<0>(x.second);
//        		  ofs_chain_ends << bliss::utils::KmerUtils::toASCIIString(L) << "\t" <<
//        				  bliss::utils::KmerUtils::toASCIIString(R) << "\t";
        	  } else {
        		  std::cout << "ERROR" << std::endl;
        		  continue;
        	  }

        	  edge_freq_type ef;
        	  ::std::get<0>(ef) = L;
        	  ::std::get<1>(ef) = R;


        	  // next print the left and right edges.
        	  cL = canonical(L);
        	  auto compact_edge = idx2.get_map().get_local_container().find(cL);
			  if (cL == L) {  // already canonical.  can use in edge directly.
				  // get in edges of L
				  std::get<0>(std::get<2>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
				  std::get<1>(std::get<2>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
				  std::get<2>(std::get<2>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
				  std::get<3>(std::get<2>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);


			  } else {   // not canonical
				  // get out edges of L, then complement each.  (reverse order)
				  std::get<0>(std::get<2>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
				  std::get<1>(std::get<2>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
				  std::get<2>(std::get<2>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
				  std::get<3>(std::get<2>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
			  }

        	  cR = canonical(R);
        	  compact_edge = R_freq_map.find(cR);  // previously retrieved from remote.
			  if (cR == R) {  // already canonical
				  // get in edges of R
				  std::get<0>(std::get<3>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
				  std::get<1>(std::get<3>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
				  std::get<2>(std::get<3>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
				  std::get<3>(std::get<3>(ef)) = compact_edge->second.get_in_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
			  } else {   // not canonical
				  // get out edges of L, then complement each.  (reverse order)
				  std::get<0>(std::get<3>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['T']);
				  std::get<1>(std::get<3>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['G']);
				  std::get<2>(std::get<3>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['C']);
				  std::get<3>(std::get<3>(ef)) = compact_edge->second.get_out_edge_frequency(KmerType::KmerAlphabet::FROM_ASCII['A']);
			  }

        	  auto fre = freq_map.get_local_container().find(cL);
        	  std::get<0>(std::get<4>(ef)) = (static_cast<float>(std::get<1>(fre->second)) /  static_cast<float>(std::get<0>(fre->second)));
        	  std::get<1>(std::get<4>(ef)) = std::get<2>(fre->second);
        	  std::get<2>(std::get<4>(ef)) = std::get<3>(fre->second);


        	  edge_freqs.emplace_back(ef);
          }
          BL_BENCH_COLLECTIVE_END(test, "gather_edge_freqs", edge_freqs.size(), comm);


		  bool no_result = mxx::all_of((edge_freqs.size() == 0), comm);
		  if (!no_result) {

			  // sort
			  BL_BENCH_START(test);
			  mxx::sort(edge_freqs.begin(), edge_freqs.end(), [](edge_freq_type const & x, edge_freq_type const & y){
				  return std::get<0>(x) < std::get<0>(y);
			  }, comm);
			  BL_BENCH_COLLECTIVE_END(test, "psort_edge_freqs", edge_freqs.size(), comm);
		  }
          // print
          BL_BENCH_START(test);
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
			  write_mpiio(compacted_chain_ends_filename, ss.str().c_str(), ss.str().length(), comm);

//			  std::ofstream ofs_chain_ends(compacted_chain_ends_filename);
//			  ofs_chain_ends << ss.str();
//			  ofs_chain_ends.close();
          }
          BL_BENCH_COLLECTIVE_END(test, "print_edge_freqs", edge_freqs.size(), comm);

//        std::cout << "COMPACTED CHAIN END POINTS" << std::endl;
//        for (auto t : chain_rep) {
//          auto md = t.second;
//          std::cout << "terminus\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) << std::endl;
//          std::cout << "\tkmer: " << bliss::utils::KmerUtils::toASCIIString(t.first) << std::endl;
//          std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) << std::endl;
//        }
      }
    }


    BL_BENCH_REPORT_MPI_NAMED(test, "output", comm);

  }




  // mpi cleanup is automatic
  comm.barrier();

  return 0;

}

