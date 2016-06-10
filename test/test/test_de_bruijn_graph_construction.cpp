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
#include "debruijn/edge_iterator.hpp"
#include "debruijn/de_bruijn_construct_engine.hpp"
#include "debruijn/de_bruijn_nodes_distributed.hpp"
#include "debruijn/de_bruijn_stats.hpp"
#include "debruijn/de_bruijn_filter.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"

#include "mxx/env.hpp"
#include "mxx/comm.hpp"



using Alphabet = bliss::common::DNA;
using KmerType = bliss::common::Kmer<31, Alphabet, WordType>;
using EdgeEncoding = bliss::common::DNA16;

#define FileParser ::bliss::io::FASTQParser
using Parser = bliss::de_bruijn::de_bruijn_parser<KmerType, EdgeEncoding>;

using MapType = bliss::de_bruijn::simple_hash_de_bruijn_map<KmerType>;

using IndexType = ::bliss::index::kmer::Index<MapType, Parser>;


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_posix(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  // default to including quality score iterators.
  IndexType::template read_file_posix<FileParser, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}

template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed, mxx::comm const & comm) {
  std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));

  size_t n_p = (n / comm.size());
  std::vector<size_t> send_counts(comm.size(), n_p);

  if (n < static_cast<size_t>(comm.size())) {
    n_p = 1;

    for (size_t i = 0; i < n; ++i) {
      send_counts[(i + comm.rank()) % comm.size()] = 1;
    }
    for (int i = n; i < comm.size(); ++i) {
      send_counts[(i + comm.rank()) % comm.size()] = 0;
    }
  }

  std::vector<KmerType> out = ::mxx::all2allv(query, send_counts, comm);
  query.swap(out);
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

	  std::string filename;
	  filename.assign(PROJ_SRC_DIR);
    filename.append("/test/data/test.debruijn.small.fastq");

	  std::string queryname(filename);

	  int sample_ratio = 100;

	  // Wrap everything in a try block.  Do this every time,
	  // because exceptions will be thrown for problems.
	  try {

	    // Define the command line object, and insert a message
	    // that describes the program. The "Command description message"
	    // is printed last in the help text. The second argument is the
	    // delimiter (usually space) and the last one is the version number.
	    // The CmdLine object parses the argv array based on the Arg objects
	    // that it contains.
	    TCLAP::CmdLine cmd("Benchmark parallel de bruijn graph building", ' ', "0.1");

	    // Define a value argument and add it to the command line.
	    // A value arg defines a flag and a type of value that it expects,
	    // such as "-n Bishop".
	    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);
	    TCLAP::ValueArg<std::string> queryArg("Q", "query", "FASTQ file path for query. default to same file as index file", false, "", "string", cmd);


	    TCLAP::ValueArg<int> sampleArg("S",
	                                 "query-sample", "sampling ratio for the query kmers. default=100",
	                                 false, sample_ratio, "int", cmd);


	    // Parse the argv array.
	    cmd.parse( argc, argv );

	    // Get the value parsed by each arg.
	    queryname = queryArg.getValue();   // get this first
	    if (queryname.empty()) // at default  set to same as input.
	      queryname = fileArg.getValue();

	    filename = fileArg.getValue();
	    sample_ratio = sampleArg.getValue();

	    // set the default for query to filename, and reparse



	    // Do what you intend.

	  } catch (TCLAP::ArgException &e)  // catch any exceptions
	  {
	    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	    exit(-1);
	  }





	  // ================  read and get file

	  IndexType idx(comm);

    if (true) {
      KmerType empty_key = ::bliss::kmer::hash::sparsehash::empty_key<KmerType>::generate();
      KmerType deleted_key = ::bliss::kmer::hash::sparsehash::deleted_key<KmerType>::generate();

	  	idx.get_map().reserve_keys(empty_key, deleted_key);

//	  	// upper key is negation of lower keys
//	  	KmerType upper_empty_key = empty_key;
//	  	KmerType upper_deleted_key = deleted_key;
//	  	for (size_t i = 0; i < KmerType::nWords; ++i) {
//	  		upper_empty_key.getDataRef()[i] = ~(upper_empty_key.getDataRef()[i]);
//	  		upper_deleted_key.getDataRef()[i] = ~(upper_deleted_key.getDataRef()[i]);
//	  	}
//
//	  	idx.get_map().reserve_upper_keys(upper_empty_key, upper_deleted_key);
	  }

	  BL_BENCH_INIT(test);

	  if (comm.rank() == 0) printf("reading query %s via posix\n", queryname.c_str());
	  BL_BENCH_START(test);
	  auto query = readForQuery_posix<IndexType>(queryname, comm);
	  BL_BENCH_COLLECTIVE_END(test, "read_query", query.size(), comm);

	  BL_BENCH_START(test);
	  sample(query, query.size() / sample_ratio, comm.rank(), comm);
	  BL_BENCH_COLLECTIVE_END(test, "sample", query.size(), comm);


	  {
		  ::std::vector<typename Parser::value_type> temp;

		  BL_BENCH_START(test);
			if (comm.rank() == 0) printf("reading %s via posix\n", filename.c_str());
			idx.read_file_posix<FileParser, Parser>(filename, temp, comm);
		  BL_BENCH_COLLECTIVE_END(test, "read", temp.size(), comm);

		  size_t total = mxx::allreduce(temp.size(), comm);
		  if (comm.rank() == 0) printf("total size is %lu\n", total);

		  BL_BENCH_START(test);
		  idx.insert(temp);
		  BL_BENCH_COLLECTIVE_END(test, "insert", idx.local_size(), comm);

	    total = idx.size();
	    if (comm.rank() == 0) printf("total size after insert/rehash is %lu\n", total);
	  }

	  {

	    {
	      auto lquery = query;
	      BL_BENCH_START(test);
	      auto counts = idx.count(lquery);
	      BL_BENCH_COLLECTIVE_END(test, "count", counts.size(), comm);
	    }
	      {
	        auto lquery = query;

	      BL_BENCH_START(test);
	      auto found = idx.find_overlap(lquery);
	      BL_BENCH_COLLECTIVE_END(test, "find_overlap", found.size(), comm);
	      }
	    // separate test because of it being potentially very slow depending on imbalance.

	      {
	        BL_BENCH_START(test);
	        ::std::vector<std::pair<typename MapType::key_type, typename MapType::mapped_type> > nodes;
	        // get the content as an array
	        idx.get_map().to_vector(nodes);
	        BL_BENCH_COLLECTIVE_END(test, "to_vector", nodes.size(), comm);

          BL_BENCH_START(test);
	        // then compute histogram
	        ::bliss::de_bruijn::print_dbg_edge_histogram(nodes, comm);
	        BL_BENCH_COLLECTIVE_END(test, "histogram", nodes.size(), comm);
	      }


	      {
          BL_BENCH_START(test);
	        // find chain nodes
	        auto nodes = idx.find_if(::bliss::de_bruijn::filter::IsChainNode());
          BL_BENCH_COLLECTIVE_END(test, "get_chains", nodes.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(nodes, comm);
          BL_BENCH_COLLECTIVE_END(test, "chain_histogram", nodes.size(), comm);

	      }

        {
          BL_BENCH_START(test);
          // find chain nodes
          auto nodes = idx.find_if(::bliss::de_bruijn::filter::IsBranchPoint());
          BL_BENCH_COLLECTIVE_END(test, "get_branches", nodes.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(nodes, comm);
          BL_BENCH_COLLECTIVE_END(test, "branch_histogram", nodes.size(), comm);

          BL_BENCH_START(test);
          std::vector<KmerType> all_neighbors;
          all_neighbors.reserve(nodes.size() * 4);

          std::vector<KmerType> neighbors;
          neighbors.reserve(4);
          for (auto t : nodes) {
            bliss::de_bruijn::node::node_utils<KmerType, typename MapType::mapped_type>::get_out_neighbors(t.first, t.second, neighbors);
            all_neighbors.insert(all_neighbors.end(), neighbors.begin(), neighbors.end());

            bliss::de_bruijn::node::node_utils<KmerType, typename MapType::mapped_type>::get_in_neighbors(t.first, t.second, neighbors);
            all_neighbors.insert(all_neighbors.end(), neighbors.begin(), neighbors.end());
          }
          BL_BENCH_COLLECTIVE_END(test, "get_branch_neighbors", all_neighbors.size(), comm);

          BL_BENCH_START(test);
          auto found = idx.find_if_overlap(all_neighbors, ::bliss::de_bruijn::filter::IsChainNode());
          BL_BENCH_COLLECTIVE_END(test, "find_chain_termini", found.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(found, comm);
          BL_BENCH_COLLECTIVE_END(test, "termini_histogram", found.size(), comm);


        }

	  }


	  BL_BENCH_REPORT_MPI_NAMED(test, "app", comm);


	  // mpi cleanup is automatic
	  comm.barrier();

	  return 0;


//
//
//
//	::std::cerr<<"Using DNA16 to present each edge" << ::std::endl;
//	testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<CountNodeMapType>, bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, count."));
//
//
//  ::std::cerr<<"Using DNA16 to represent each edge" << ::std::endl;
//  testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<ExistNodeMapType>,  bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, exist."));


}

