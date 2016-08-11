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
#include "debruijn/debruijn_graph_loader.hpp"
#include "debruijn/debruijn_stats.hpp"
#include "debruijn/debruijn_graph_filters.hpp"
#include "debruijn/debruijn_chain_filters.hpp"
#include "debruijn/debruijn_chain_node.hpp"
#include "debruijn/debruijn_chain_operations.hpp"

#include "debruijn/debruijn_common.hpp"
#include "debruijn/debruijn_graph_map.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"

#include "mxx/env.hpp"
#include "mxx/comm.hpp"

#if (pDNA == 16)
using Alphabet = bliss::common::DNA16;
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
//#elif (pDNA == 4)
//using Alphabet = bliss::common::DNA;
#endif


using KmerType = bliss::common::Kmer<31, Alphabet, WordType>;
using EdgeEncoding = Alphabet;

#define FileParser ::bliss::io::FASTQParser

using DBGNodeParser = bliss::debruijn::debruijn_graph_parser<KmerType>;

using DBGMapType = ::bliss::debruijn::graph::simple_hash_compact_debruijn_graph_map<KmerType>;
using DBGType = ::bliss::index::kmer::Index<DBGMapType, DBGNodeParser>;

using ChainNodeType = ::bliss::debruijn::simple_biedge<KmerType>;
//template <typename K>
//using ChainMapParams = ::bliss::index::kmer::CanonicalHashMapParams<K>;
using ChainMapType = ::dsc::densehash_map<KmerType, ChainNodeType,
    ::bliss::debruijn::CanonicalDeBruijnHashMapParams,
     ::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

using ChainVecType = ::std::vector<std::pair<KmerType, ChainNodeType> >;



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
    TCLAP::CmdLine cmd("Parallel de bruijn graph building", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);

    TCLAP::ValueArg<std::string> queryArg("O", "output_prefix", "Output file prefix", false, "", "string", cmd);

//    TCLAP::ValueArg<std::string> queryArg("Q", "query", "FASTQ file path for query. default to same file as index file", false, "", "string", cmd);
//    TCLAP::ValueArg<int> sampleArg("S",
//                                   "query-sample", "sampling ratio for the query kmers. default=100",
//                                   false, sample_ratio, "int", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    filename = fileArg.getValue();

//    // Get the value parsed by each arg.
// set the default for query to filename, and/home/DATA/1000genome/HG00096/SRR077487_1.filt.0_03125noN.fastq reparse
//    queryname = queryArg.getValue();   // get this first
//    if (queryname.empty()) // at default  set to same as input.
//      queryname = filename;
//    sample_ratio = sampleArg.getValue();


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


  KmerType testKmer(std::string("CAAGATGGGTGGAATGGCCAGTTAACCACTG"));

  // TODO: filter out, or do something, about "N".  May have to add back support for ASCII edge encoding so that we can use DNA5 alphabet

  {
    ::std::vector<typename DBGNodeParser::value_type> temp;

    BL_BENCH_START(test);
    if (comm.rank() == 0) printf("reading %s via posix\n", filename.c_str());
    idx.read_file_posix<FileParser, DBGNodeParser>(filename, temp, comm);
    BL_BENCH_COLLECTIVE_END(test, "read", temp.size(), comm);

    //		  for (auto t : temp) {
    //		    std::cout << "input kmer " << bliss::utils::KmerUtils::toASCIIString(t.first) << " edges " << t.second << std::endl;
    //		  }
    //

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
  }

  {

    //	      {
    //	        auto lquery = query;
    //
    //	      BL_BENCH_START(test);
    //	      auto found = idx.find_overlap(lquery);
    //	      BL_BENCH_COLLECTIVE_END(test, "find_overlap", found.size(), comm);
    //	      }
    //	    // separate test because of it being potentially very slow depending on imbalance.

    {
      // get histogram for the edge types in debruijn graph

      ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > nodes;

      BL_BENCH_START(test);
      // get the content as an array
      idx.get_map().to_vector(nodes);
      BL_BENCH_COLLECTIVE_END(test, "to_vector", nodes.size(), comm);

      BL_BENCH_START(test);
      // then compute histogram
      ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(nodes.begin(), nodes.end(), comm);
      BL_BENCH_COLLECTIVE_END(test, "histogram", nodes.size(), comm);
    }


    {
      ::std::vector<std::pair<typename DBGMapType::key_type, typename DBGMapType::mapped_type> > chain_nodes;

      BL_BENCH_START(test);
      // find chain nodes
      chain_nodes = idx.find_if(::bliss::debruijn::filter::graph::IsChainNode());
      BL_BENCH_COLLECTIVE_END(test, "get_chains", chain_nodes.size(), comm);

      BL_BENCH_START(test);
      // then compute histogram
      ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(chain_nodes.begin(), chain_nodes.end(), comm);
      BL_BENCH_COLLECTIVE_END(test, "chain_histogram", chain_nodes.size(), comm);


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


        chainmap.get_local_container().insert(::std::make_pair(::std::move(t.first), ::std::move(node)));
        //            std::cout << "BIEDGE\tin dist\t" << std::get<2>(node) << " kmer: " << std::get<0>(node) << std::endl;
        //            std::cout << "\tKmer bi edge:\t" << t.first << std::endl;
        //            std::cout << "\tout dist\t" << std::get<3>(node) << " kmer: " << std::get<1>(node) << std::endl;

      }
      BL_BENCH_COLLECTIVE_END(test, "insert in chainmap.", chainmap.local_size(), comm);

      //========= report.
      auto result = chainmap.find(::bliss::debruijn::filter::chain::IsTerminus());
      auto result2 = chainmap.find(::bliss::debruijn::filter::chain::IsIsolated());
      printf("chain map contains %lu chained termini and  %lu isolated\n", result.size(), result2.size());

    }


    {
      //=== find branching nodes. local computation.
      BL_BENCH_START(test);
      auto nodes = idx.find_if(::bliss::debruijn::filter::graph::IsBranchPoint());
      BL_BENCH_COLLECTIVE_END(test, "get_branches", nodes.size(), comm);

      BL_BENCH_START(test);
      // then compute histogram
      ::bliss::debruijn::graph::print_compact_multi_biedge_histogram(nodes.begin(), nodes.end(), comm);
      BL_BENCH_COLLECTIVE_END(test, "branch_histogram", nodes.size(), comm);

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

        // get global unfinished count
        cycle_nodes = std::count_if(unfinished.begin(), unfinished.end(),
                                    ::bliss::debruijn::filter::chain::IsCycleNode(iterations));
        if (comm.rank() < 4) printf("rank %d iter %lu updated %lu, unfinished %lu noncycle %lu\n", comm.rank(), iterations, count, unfinished.size(), cycle_nodes);
        all_compacted = (count == 0) || (cycle_nodes == unfinished.size());
        all_compacted = ::mxx::all_of(all_compacted, comm);

      }
      BL_BENCH_COLLECTIVE_END(test, "compact", cycle_nodes, comm);
    }

    // ==== final query to see if at terminus
    {
      // NOW: do the list ranking

      // search unfinished
      BL_BENCH_START(test);
      auto unfinished = chainmap.find(::bliss::debruijn::filter::chain::PointsToInternalNode());
      BL_BENCH_COLLECTIVE_END(test, "unfinished2", unfinished.size(), comm);

      // get global unfinished count
      bool all_compacted = (unfinished.size() == 0);
      all_compacted = ::mxx::all_of(all_compacted, comm);


      // while have unfinished,  run.  qq contains kmers not necessarily canonical
      BL_BENCH_START(test);

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

          // lookup left by canonicaliter %ld
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
      unfinished = chainmap.find(::bliss::debruijn::filter::chain::IsUncompactedNode(iterations));


      // get global unfinished count
      all_compacted = (unfinished.size() == 0);
      all_compacted = ::mxx::all_of(all_compacted, comm);

      assert(all_compacted);
    }


    {
      // ========= get compacted chains.

      // remove cycle nodes
      BL_BENCH_START(test);
      size_t cycle_node_count = chainmap.erase(::bliss::debruijn::filter::chain::IsCycleNode(iterations));
      BL_BENCH_COLLECTIVE_END(test, "erase cycle", cycle_node_count, comm);

      BL_BENCH_START(test);
      auto termini = chainmap.find(::bliss::debruijn::filter::chain::IsTerminus());
      BL_BENCH_COLLECTIVE_END(test, "chain termini", termini.size(), comm);

      std::cout << "COMPACTED CHAIN END POINTS" << std::endl;
      for (auto t : termini) {
        auto md = t.second;
        std::cout << "terminus\tin dist " << std::get<2>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md)) << " rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(md).reverse_complement()) << std::endl;
        std::cout << "\tkmer: " << bliss::utils::KmerUtils::toASCIIString(t.first) << " rc: " << bliss::utils::KmerUtils::toASCIIString(t.first.reverse_complement()) << std::endl;
        std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md)) << " rc: " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(md).reverse_complement()) << std::endl;
      }

      // ========== construct new graph with compacted chains and junction nodes.

      BL_BENCH_START(test);
      std::vector<::bliss::debruijn::chain::compacted_chain_node<KmerType> > result;
      result.reserve(chainmap.size());
      ::fsc::back_emplace_iterator<std::vector<::bliss::debruijn::chain::compacted_chain_node<KmerType> > > back_emplacer(result);

      // first transform nodes so that we are pointing to canonical terminus k-mers.
      std::transform(chainmap.get_local_container().begin(), chainmap.get_local_container().end(), back_emplacer,
    		  ::bliss::debruijn::operation::chain::to_compacted_chain_node<KmerType>());
      BL_BENCH_COLLECTIVE_END(test, "transform chain", chainmap.size(), comm);

//      for (auto r : result) {
//    	  std::cout << "k-mer: " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(r)) <<
//    			  " l-mer id " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(r)) << " dist " << std::get<2>(r) << std::endl;
 //     }

      // sort
      BL_BENCH_START(test);
      // first transform nodes so that we are pointing to canonical terminus k-mers.
      std::sort(result.begin(), result.end(), ::bliss::debruijn::operation::chain::chain_rep_less<KmerType>());
      BL_BENCH_COLLECTIVE_END(test, "sort chain node", result.size(), comm);

      // print out.
      BL_BENCH_START(test);
      std::for_each(result.begin(), result.end(), ::bliss::debruijn::operation::chain::print_chain<KmerType>(std::cout));
      BL_BENCH_COLLECTIVE_END(test, "print chain", result.size(), comm);

    }


  }


  BL_BENCH_REPORT_MPI_NAMED(test, "app", comm);


  // mpi cleanup is automatic
  comm.barrier();

  return 0;

}

