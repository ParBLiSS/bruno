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
#ifndef COMPACT_DBG_STATS_CPP
#define COMPACT_DBG_STATS_CPP


// ==========  choices:
//  1. no filtering:  parse simple nodes directly, insert
//  2. de novo filtering:  parse simple nodes, parse edge k1mers, transform edges, copy edges, filter and erase, insert
//  3. reconstructive filtering:  parser kmers, add edges, filter and copy, insert.

// FILTERING ONLY WORKS WITH FASTQ files right now.  TODO: STILL NEEDED?
//#if (pPARSER == FASTQ)



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

		printf("check query is superset of content:  total query = %lu, erased = %lu, remaining = %lu\n", query.size(), erased, idx_copy.local_size());
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
	std::vector<std::pair<size_t, bool> >  read_ends;

	using LocalCountMapType = ::dsc::counting_densehash_map<KmerType, size_t,
			FreqMapParams,
			::bliss::kmer::hash::sparsehash::special_keys<KmerType, true> >;

	typename LocalCountMapType::local_container_type counter;
	typename LocalCountMapType::local_container_type::const_iterator it;

	std::stringstream ss;
  std::stringstream ss2;

	kmers.clear();
	read_ends.clear();

	ss.str(std::string());
  ss2.str(std::string());
	BL_BENCH_COLLECTIVE_END(valid_print, "init", fdata.getRange().size(), comm);

	// note that we are not using a split sequence iterator here, or a filtering sequence iterator, sicne we NEED to identify the actual first entry.

	BL_BENCH_START(valid_print);
	// get the kmers
  ::bliss::io::KmerFileHelper::template parse_file_data<::bliss::index::kmer::KmerParser<KmerType>, FileParser, SeqIterType>(fdata, kmers, comm);
	std::transform(kmers.begin(), kmers.end(), kmers.begin(), canonical);
	BL_BENCH_COLLECTIVE_END(valid_print, "parse kmers", kmers.size(), comm);


	BL_BENCH_START(valid_print);
	// get the read lengths
  ::bliss::io::KmerFileHelper::template parse_file_data_old<::bliss::debruijn::ReadLengthParser<KmerType>, FileParser, SeqIterType>(fdata, read_ends, comm);
	// inclusive prefix scan to get the offset to the end of the reads.
	for (size_t i = 1; i < read_ends.size(); ++i) {
		read_ends[i].first += read_ends[i-1].first;
	}
	// global prefix scan.  only for procs that have data.
//	::mxx::comm subcomm = comm.split(read_ends.size() > 0);
//	size_t global_offset = 0;
//	if (read_ends.size() > 0) {
//		global_offset = ::mxx::exscan(read_ends.back(), subcomm);
//	}
	BL_BENCH_COLLECTIVE_END(valid_print, "read size", read_ends.size(), comm);


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
	//std::transform(kmers.begin(), kmers.end(), kmers.begin(), canonical);
	BL_BENCH_COLLECTIVE_END(valid_print, "reparse", kmers.size(), comm);


	BL_BENCH_START(valid_print);
	// linear scan to output the valid positions
	int64_t rstart = 0;
	int64_t rend, vstart, vend;
	size_t len = 0;

	// exclusive prefix scan to get correct read id.
	size_t prev_i = read_ends.size();
//	if (comm.rank() == comm.size() - 1) {
//	  std::cerr << "rank " << comm.rank() << " size " << prev_i << " -> ";
//	}
	prev_i = mxx::exscan(prev_i, comm);
//  if (comm.rank() == comm.size() - 1) {
//    std::cerr << " exscan " << prev_i << std::endl;
//  }

	for (size_t i = 0; i < read_ends.size(); ++i) {
		rend = static_cast<int64_t>(read_ends[i].first);

		vstart = -1;
		vend = -1;
		len = (i == 0 ? read_ends[i].first : read_ends[i].first - read_ends[i-1].first);
		ss2 << "read " << (i+1+prev_i) << " local endpos " << read_ends[i].first << " len " << len << std::flush;


		//std::cout << "global offset " << global_offset << " local offsets : " << rstart << " - " << rend << std::endl;
    if ((read_ends[i].second == true) && (len > 0)) {  // only compute if there is no N.

      for (int64_t j = rstart; j < rend; ++j) {


        it = counter.find(canonical(kmers[j]));

        if (it != counter.end()) {
  //				std::cout << "rank " << comm.rank() << " pos " << j << " rstart " << rstart <<
  //						" query " << bliss::utils::KmerUtils::toASCIIString(kmers[j]) <<
  //						" result " << bliss::utils::KmerUtils::toASCIIString((*it).first) <<
  //						" start count " << (*it).second << std::endl;
          // kmer exists in the debruijn graph
          if ((*it).second > 0) {
            //print the first pos.
            vstart = j - rstart;
            ss2 << " local rstart: " << rstart << " pos " << vstart <<
                " " << bliss::utils::KmerUtils::toASCIIString(kmers[j]) <<
                " " << bliss::utils::KmerUtils::toASCIIString((*it).first) << std::flush;
            break;
          }
        }
      }
      ss2 << " <=> " << std::flush;

      if (vstart > -1) {
        for (int64_t j = rend - 1; j >= rstart; --j) {
          it = counter.find(canonical(kmers[j]));

          if (it != counter.end()) {
  //					std::cout << "rank " << comm.rank() << " pos " << j << " rstart " << rstart <<
  //							" query " << bliss::utils::KmerUtils::toASCIIString(kmers[j]) <<
  //							" result " << bliss::utils::KmerUtils::toASCIIString((*it).first) <<
  //							" end count " << (*it).second << std::endl;
            if ((*it).second > 0) {
              vend = j - rstart;
              ss2 << " local rstart: " << rstart << " pos " << vend <<
                  " " << bliss::utils::KmerUtils::toASCIIString(kmers[j]) <<
                  " " << bliss::utils::KmerUtils::toASCIIString((*it).first) << std::flush;
              break;
            }
          }
        }
      }

    }  // only compute start and end if there is no N.
    ss2 << std::endl;

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
	BL_BENCH_COLLECTIVE_END(valid_print, "find start-end", read_ends.size(), comm);


	// write out to file
	BL_BENCH_START(valid_print);
	write_mpiio(filename, ss.str().c_str(), ss.str().length(), comm);

	BL_BENCH_COLLECTIVE_END(valid_print, "print", ss.str().length(), comm);

  BL_BENCH_START(valid_print);
  write_mpiio(filename + ".debug", ss2.str().c_str(), ss2.str().length(), comm);

  BL_BENCH_COLLECTIVE_END(valid_print, "debug print", ss2.str().length(), comm);

	BL_BENCH_REPORT_MPI_NAMED(valid_print, "read_pos_print", comm);

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
				[&counter](::bliss::debruijn::biedge::compact_simple_biedge_kmer_node<KmerType> const & y){
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

	BL_BENCH_INIT(chain_freq);

	// ==  first compute frequency summary, and store into a reduction map
	// allocate input
	std::vector< std::pair<KmerType, FreqSummaryType > > freqs;

  ::bliss::debruijn::lex_less<KmerType> canonical;   // count_idx is distributed same way as the DBG, which is canonical here.  so to query via local container, need to canoncialize the compacted_chain content.

	BL_BENCH_START(chain_freq);
	// extract frequencies.
	for (auto x : compacted_chain) {
		// compute the chain rep
		CountType c = count_idx.get_map().get_local_container().find(canonical(std::get<0>(x)))->second;

		// new key is the chain rep.
		freqs.emplace_back(std::get<1>(x), FreqSummaryType(1, c, c, c));
	}
	BL_BENCH_COLLECTIVE_END(chain_freq, "get_node_freqs", freqs.size(), comm);

	// create a reduction map
	BL_BENCH_START(chain_freq);
	chain_freq_map.insert(freqs);   // collective comm.
	BL_BENCH_COLLECTIVE_END(chain_freq, "reduce_freq", chain_freq_map.local_size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(chain_freq, "chain counts", comm);
}

void compute_freq_map_incremental(ListRankedChainNodeVecType const & compacted_chain,
CountIndexType const & count_idx,
FreqMapType & chain_freq_map,
mxx::comm const & comm) {


	BL_BENCH_INIT(chain_freq);

	// ==  first compute frequency summary, and store into a reduction map
	// allocate input
	std::vector< std::pair<KmerType, FreqSummaryType > > freqs;

	// estimate the largest amount of memory to use.
	unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);

	// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
	size_t step = (free_mem / (8 * sizeof(std::pair<KmerType, FreqSummaryType >)));  // number of elements that can be held in freemem
  step = std::max(std::min(step, compacted_chain.size()), static_cast<size_t>(1));

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

/// compute frequency of chains.  compacted chain should have same distribution as would be for count index.
void compute_freq_map(ListRankedChainNodeVecType const & compacted_chain,
		CountDBGType const & count_idx,
		FreqMapType & chain_freq_map,
		mxx::comm const & comm) {

	BL_BENCH_INIT(chain_freq);

	// ==  first compute frequency summary, and store into a reduction map
	// allocate input
	std::vector< std::pair<KmerType, FreqSummaryType > > freqs;

  ::bliss::debruijn::lex_less<KmerType> canonical;   // count_idx is distributed same way as the DBG, which is canonical here.  so to query via local container, need to canoncialize the compacted_chain content.


	BL_BENCH_START(chain_freq);
	typename CountDBGType::map_type::mapped_type e;
	CountType c = 0;
	// extract frequencies.
	for (auto x : compacted_chain) {
		// compute the chain rep
		c = count_idx.get_map().get_local_container().find(canonical(std::get<0>(x)))->second.get_self_frequency();

		// new key is the chain rep.
		freqs.emplace_back(std::get<1>(x), FreqSummaryType(1, c, c, c));
	}
	BL_BENCH_COLLECTIVE_END(chain_freq, "get_node_freqs", freqs.size(), comm);

	// create a reduction map
	BL_BENCH_START(chain_freq);
	chain_freq_map.insert(freqs);   // collective comm.
	BL_BENCH_COLLECTIVE_END(chain_freq, "reduce_freq", chain_freq_map.local_size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(chain_freq, "chain counts", comm);
}

void compute_freq_map_incremental(ListRankedChainNodeVecType const & compacted_chain,
CountDBGType const & count_idx,
FreqMapType & chain_freq_map,
mxx::comm const & comm) {


	BL_BENCH_INIT(chain_freq);

	// ==  first compute frequency summary, and store into a reduction map
	// allocate input
	std::vector< std::pair<KmerType, FreqSummaryType > > freqs;

	// estimate the largest amount of memory to use.
	unsigned long free_mem = ::utils::get_free_mem_per_proc(comm);

	// use 1/8 of space, local 1x, remote 1x, insert 1x, rest is just to be conservative.  this is assuming input is evenly distributed.
	size_t step = (free_mem / (8 * sizeof(std::pair<KmerType, FreqSummaryType >)));  // number of elements that can be held in freemem
	step = std::max(std::min(step, compacted_chain.size()), static_cast<size_t>(1));

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
			c = count_idx.get_map().get_local_container().find(canonical(std::get<0>(*iter)))->second.get_self_frequency();

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
					static_cast<size_t>(std::get<0>(std::get<2>(x))) << "\t" <<
					static_cast<size_t>(std::get<1>(std::get<2>(x))) << "\t" <<
					static_cast<size_t>(std::get<2>(std::get<2>(x))) << "\t" <<
					static_cast<size_t>(std::get<3>(std::get<2>(x))) << "\t" <<
					static_cast<size_t>(std::get<0>(std::get<3>(x))) << "\t" <<
					static_cast<size_t>(std::get<1>(std::get<3>(x))) << "\t" <<
					static_cast<size_t>(std::get<2>(std::get<3>(x))) << "\t" <<
					static_cast<size_t>(std::get<3>(std::get<3>(x))) << "\t" <<
					static_cast<size_t>(std::get<0>(std::get<4>(x))) << "\t" <<
					static_cast<size_t>(std::get<1>(std::get<4>(x))) << "\t" <<
					static_cast<size_t>(std::get<2>(std::get<4>(x))) << "\t" << std::endl;

		}
		write_mpiio(filename, ss.str().c_str(), ss.str().length(), comm);

	}
	BL_BENCH_COLLECTIVE_END(print_chain_freq, "print_edge_freqs", edge_freqs.size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(print_chain_freq, "print_chain_freqs", comm);

}



#endif // COMPACT_DBG_STATS_CPP
