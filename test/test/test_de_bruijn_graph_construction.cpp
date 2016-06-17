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
#include "index/kmer_hash.hpp"
#include "debruijn/edge_iterator.hpp"
#include "debruijn/de_bruijn_construct_engine.hpp"
#include "debruijn/de_bruijn_nodes_distributed.hpp"
#include "debruijn/de_bruijn_stats.hpp"
#include "debruijn/de_bruijn_filter.hpp"
#include "debruijn/de_bruijn_operations.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"

#include "mxx/env.hpp"
#include "mxx/comm.hpp"



using Alphabet = bliss::common::DNA;
using KmerType = bliss::common::Kmer<31, Alphabet, WordType>;
using EdgeEncoding = bliss::common::DNA16;

#define FileParser ::bliss::io::FASTQParser

using DBGNodeParser = bliss::de_bruijn::de_bruijn_parser<KmerType, EdgeEncoding>;

using DBGMapType = ::bliss::de_bruijn::simple_hash_de_bruijn_map<KmerType>;
using DBGType = ::bliss::index::kmer::Index<DBGMapType, DBGNodeParser>;

using ChainNodeType = ::bliss::de_bruijn::operation::chain::compaction_metadata<KmerType>;
//template <typename K>
//using ChainMapParams = ::bliss::index::kmer::CanonicalHashMapParams<K>;
using ChainMapType = ::dsc::densehash_map<KmerType, ChainNodeType,
		::bliss::de_bruijn::operation::chain::CanonicalDeBruijnChainMapParams,
    ::bliss::kmer::hash::sparsehash::special_keys<KmerType> >;

using ChainVecType = ::std::vector<std::pair<KmerType, ChainNodeType> >;


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

	    // set the default for query to filename, and/home/DATA/1000genome/HG00096/SRR077487_1.filt.0_03125noN.fastq reparse



	    // Do what you intend.

	  } catch (TCLAP::ArgException &e)  // catch any exceptions
	  {
	    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	    exit(-1);
	  }





	  // ================  read and get file

	  DBGType idx(comm);
	  ChainMapType chainmap(comm);


	  BL_BENCH_INIT(test);

	  std::vector<KmerType> query;
    std::vector<KmerType> neighbors;
    neighbors.reserve(4);


    // TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet

	  {
		  ::std::vector<typename DBGNodeParser::value_type> temp;

		  BL_BENCH_START(test);
			if (comm.rank() == 0) printf("reading %s via posix\n", filename.c_str());
			idx.read_file_posix<FileParser, DBGNodeParser>(filename, temp, comm);
		  BL_BENCH_COLLECTIVE_END(test, "read", temp.size(), comm);


		  // all possible k-mers from input should already be present, so no need to test.

      // ============== insert

		  size_t total = mxx::allreduce(temp.size(), comm);
		  if (comm.rank() == 0) printf("total size is %lu\n", total);

		  BL_BENCH_START(test);
		  idx.insert(temp);
		  BL_BENCH_COLLECTIVE_END(test, "insert", idx.local_size(), comm);

	    total = idx.size();
	    if (comm.rank() == 0) printf("total size after insert/rehash is %lu\n", total);


	    // get all possible edges.
	      //============== testing to ensure that all the possible edges are present.
	    auto cc = idx.get_map().get_local_container();
	    for (auto it = cc.begin(); it != cc.end(); ++it) {
	    	neighbors.clear();
	    	bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_out_neighbors(it->first, it->second, neighbors);
	        query.insert(query.end(), neighbors.begin(), neighbors.end());

	        neighbors.clear();
	        bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_in_neighbors(it->first, it->second, neighbors);
	        query.insert(query.end(), neighbors.begin(), neighbors.end());
	    }

	    // =============== check to see if index is superset of query.  (count should have every result entry showing 1.)
	    {
	      auto lquery = query;
	      BL_BENCH_START(test);
	      auto counts = idx.count(lquery);
	      BL_BENCH_COLLECTIVE_END(test, "count", counts.size(), comm);

	      auto absent_end = std::partition(counts.begin(), counts.end(), [](std::pair<KmerType, size_t> const & x){
	    	  return x.second == 0;
	      });
	      printf(" total query = %lu, unique query = %lu, unique absent = %lu\n", query.size(), counts.size(), std::distance(counts.begin(), absent_end));

	      for (auto it = counts.begin(); it != absent_end; ++it) {
	    	  std::cout << "  " << it->first << std::endl;
 	      }
	    }

	    comm.barrier();

	    // =============== check to see if query is superset of index.  (erase should result in empty)
	    {
	    	auto lquery = query;
		      BL_BENCH_START(test);
		      size_t erased = idx.get_map().erase(lquery);
		      BL_BENCH_COLLECTIVE_END(test, "erase", erased, comm);

		      printf(" total query = %lu, erased = %lu, remaining = %lu\n", query.size(), erased, idx.local_size());

	    }

	      comm.barrier();
	      idx.get_map().erase(::fsc::TruePredicate());
	      comm.barrier();
		  BL_BENCH_START(test);
		  idx.insert(temp);
		  BL_BENCH_COLLECTIVE_END(test, "reinsert", idx.local_size(), comm);
	  }

	  {

	      {
	        auto lquery = query;

	      BL_BENCH_START(test);
	      auto found = idx.find_overlap(lquery);
	      BL_BENCH_COLLECTIVE_END(test, "find_overlap", found.size(), comm);
	      }
	    // separate test because of it being potentially very slow depending on imbalance.

	      {
	        // get histogram for the edge types in debruijn graph

	        ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > nodes;

	        BL_BENCH_START(test);
	        // get the content as an array
	        idx.get_map().to_vector(nodes);
	        BL_BENCH_COLLECTIVE_END(test, "to_vector", nodes.size(), comm);

          BL_BENCH_START(test);
	        // then compute histogram
	        ::bliss::de_bruijn::print_dbg_edge_histogram(nodes, comm);
	        BL_BENCH_COLLECTIVE_END(test, "histogram", nodes.size(), comm);
	      }


	      {
          ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > chain_nodes;

          BL_BENCH_START(test);
	        // find chain nodes
	        chain_nodes = idx.find_if(::bliss::de_bruijn::filter::IsChainNode());
          BL_BENCH_COLLECTIVE_END(test, "get_chains", chain_nodes.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(chain_nodes, comm);
          BL_BENCH_COLLECTIVE_END(test, "chain_histogram", chain_nodes.size(), comm);


          // insert into local container inside chainmap.
          BL_BENCH_START(test);

          // initialize the chain map.  note that key kmers are canonical, same as in DBG.
          // note also that edge k-mers have same orientation as key kmers, but may not be canonical.
          for (auto t : chain_nodes) {
            ChainNodeType node(KmerType(), KmerType(), 0, 0);   // default node

            // get the in neighbor
            neighbors.clear();
            bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_in_neighbors(t.first, t.second, neighbors);
            assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
            if (neighbors.size() == 1) {
              std::get<0>(node) = neighbors[0];
              std::get<2>(node) = 1;
            }

            // get the out neighbor
            neighbors.clear();
            bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_out_neighbors(t.first, t.second, neighbors);
            assert(neighbors.size() < 2);   // should not have more than 1 neighbors.
            if (neighbors.size() == 1) {
              std::get<1>(node) = neighbors[0];
              std::get<3>(node) = 1;
            }

            chainmap.get_local_container().insert(::std::make_pair(::std::move(t.first), ::std::move(node)));
          }
          BL_BENCH_COLLECTIVE_END(test, "insert in chainmap.", chainmap.local_size(), comm);

          //========= report.
          auto result = chainmap.find(::bliss::de_bruijn::filter::chain::IsTerminus());
          auto result2 = chainmap.find(::bliss::de_bruijn::filter::chain::IsIsolated());
          printf("chain map contains %lu chained termini and  %lu isolated\n", result.size(), result2.size());

	      }


	      ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > chain_termini;
        {
          //=== find branching nodes. local computation.
          BL_BENCH_START(test);
          auto nodes = idx.find_if(::bliss::de_bruijn::filter::IsBranchPoint());
          BL_BENCH_COLLECTIVE_END(test, "get_branches", nodes.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(nodes, comm);
          BL_BENCH_COLLECTIVE_END(test, "branch_histogram", nodes.size(), comm);

          //=== get the neighbors of the branch points.  for information only.
          BL_BENCH_START(test);
          std::vector<KmerType> all_neighbors2;
          all_neighbors2.reserve(nodes.size() * 4);

          for (auto t : nodes) {
            neighbors.clear();
            bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_out_neighbors(t.first, t.second, neighbors);
            all_neighbors2.insert(all_neighbors2.end(), neighbors.begin(), neighbors.end());

            neighbors.clear();
            bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_in_neighbors(t.first, t.second, neighbors);
            all_neighbors2.insert(all_neighbors2.end(), neighbors.begin(), neighbors.end());
          }
          BL_BENCH_COLLECTIVE_END(test, "branch_neighbors", all_neighbors2.size(), comm);

          // now check to see which are chain nodes.  these are chain nodes adjacent to branch points.
          // include chain termini that are adjacent to branch points, so we can mark them in the chainmap.
          BL_BENCH_START(test);
          auto found = idx.find_if_overlap(all_neighbors2, ::bliss::de_bruijn::filter::IsChainNode());
          BL_BENCH_COLLECTIVE_END(test, "terminal_neighbors", found.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(found, comm);
          BL_BENCH_COLLECTIVE_END(test, "termini_histogram", found.size(), comm);


          // now check to see which are chain nodes.  these are chain nodes adjacent to branch points.
          // include chain termini that are adjacent to branch points, so we can mark them in the chainmap.
          BL_BENCH_START(test);
          auto found2 = idx.find_if_overlap(all_neighbors2, ::bliss::de_bruijn::filter::IsTerminus());
          BL_BENCH_COLLECTIVE_END(test, "stump_termini", found2.size(), comm);

          BL_BENCH_START(test);
          // then compute histogram
          ::bliss::de_bruijn::print_dbg_edge_histogram(found2, comm);
          BL_BENCH_COLLECTIVE_END(test, "stump_histogram", found2.size(), comm);


          //=== mark neighbors of branch points.
          BL_BENCH_START(test);
          std::vector<std::pair<KmerType, bliss::de_bruijn::operation::chain::terminus_update_md<KmerType> > > all_neighbors;
          all_neighbors.reserve(nodes.size() * 4);

          for (auto t : nodes) {
            neighbors.clear();
            bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_out_neighbors(t.first, t.second, neighbors);
            for (auto n : neighbors) {
//              if (n <= n.reverse_complement())  {  // already canonical.  insert source kmer as is, as in neighbor to n.
                all_neighbors.emplace_back(n, bliss::de_bruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::de_bruijn::operation::IN));
//              } else {  // not canonical.  insert revcomp of source kmer, as out neighbor to n.
//                all_neighbors.emplace_back(n, bliss::de_bruijn::operation::chain::terminus_update_md<KmerType>(t.first.reverse_complement(), bliss::de_bruijn::operation::OUT));
//              }
            }

            neighbors.clear();
            bliss::de_bruijn::node::node_utils<KmerType, typename DBGMapType::mapped_type>::get_in_neighbors(t.first, t.second, neighbors);
            for (auto n : neighbors) {
//              if (n <= n.reverse_complement())  {  // already canonical.  insert source kmer as is, as out neighbor to n.
                all_neighbors.emplace_back(n, bliss::de_bruijn::operation::chain::terminus_update_md<KmerType>(t.first, bliss::de_bruijn::operation::OUT));
//              } else {  // not canonical.  insert source kmer revcomp, as in neighbor to n.
//                all_neighbors.emplace_back(n, bliss::de_bruijn::operation::chain::terminus_update_md<KmerType>(t.first.reverse_complement(), bliss::de_bruijn::operation::IN));
//              }
            }
          }
          BL_BENCH_COLLECTIVE_END(test, "branch_neighbors_2", all_neighbors.size(), comm);

          printf("chainmap size = %ld\n", chainmap.size());

          // now check to see which are chain nodes.  these are chain nodes adjacent to branch points.
          // include chain termini that are adjacent to branch points, so we can mark them in the chainmap.
          BL_BENCH_START(test);
          ::bliss::de_bruijn::operation::chain::terminus_update<KmerType> updater;
          size_t count = chainmap.update(all_neighbors, false, updater );
          BL_BENCH_COLLECTIVE_END(test, "update_termini", count, comm);

          //========= split singleton entries from chainmap.
          BL_BENCH_START(test);
          count = chainmap.erase(::bliss::de_bruijn::filter::chain::IsIsolated());
          BL_BENCH_COLLECTIVE_END(test, "clear_singleton2", count, comm);

          auto result = chainmap.find(::bliss::de_bruijn::filter::chain::IsTerminus());
          printf("chain map now contains %lu chained termini, after removing %lu isolated\n", result.size(), count);
        }


        {
          // NOW: do the list ranking












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

