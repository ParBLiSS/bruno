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
#ifndef COMPACT_DBG_IO_CPP
#define COMPACT_DBG_IO_CPP


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

	int count;
	MPI_Get_count(&stat, MPI_BYTE, &count);

	// if other I/O error, and MPI_ERROR is 0, check available space.
	ss << " MPI_Status error: [ err=" << stat.MPI_ERROR << " cnt=" << count << "]" << std::endl;

	return ss.str();
}


void write_mpiio(std::string const & filename, const char* data, size_t len, mxx::comm const & comm ) {
	// TODO: subcommunicator to work with only nodes that have data.

	/// MPI file handle
	MPI_File fh;

	bool has_data = len > 0;
	bool all_has_data = mxx::all_of(has_data, comm);
	mxx::comm subcomm = all_has_data ? comm.copy() : comm.split(has_data);

	if (has_data) {

	int res = MPI_File_open(subcomm, const_cast<char *>(filename.c_str()), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if (res != MPI_SUCCESS) {
		throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "open", res, subcomm));
	}

	res = MPI_File_set_size(fh, 0);
	if (res != MPI_SUCCESS) {
		throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "truncate", res, subcomm));
	}

	// ensure atomicity is turned off
	MPI_File_set_atomicity(fh, 0);

	// get the global offset.
	MPI_Offset global_offset = ::mxx::exscan(len, subcomm);


	int step = (0x1 << 30);
	int iterations = (len + step - 1) / step;

	//std::cout << "rank " << comm.rank() << " mpiio write offset is " << global_offset << " len " << len << " iterations " << iterations << ::std::endl;

	// get the maximum number of iterations
	iterations = ::mxx::allreduce(iterations, [](int const & x, int const & y){
		return (x >= y) ? x : y;
	}, subcomm);

#ifndef NDEBUG
	printf("rank %d subcomm rank %d write mpiio. len %ld offset %lld step %d iterations %d\n", comm.rank(), subcomm.rank(), len, global_offset, step, iterations);
#endif

	int remainder = len;
	int curr_step = step;
	MPI_Status stat;
	int count = 0;
	bool success = true;
	for (int i = 0; i < iterations; ++i) {
		curr_step = std::min(remainder, step);

		res = MPI_File_write_at_all( fh, global_offset, const_cast<char*>(data), curr_step, MPI_BYTE, &stat);  /// use bytes, no endian issues.

		success = ::mxx::all_of(res == MPI_SUCCESS, subcomm);
		if (!success) {
		  MPI_File_close(&fh);

		  if (res != MPI_SUCCESS)
		    throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "write", res, stat, subcomm));
		}

		res = MPI_Get_count(&stat, MPI_BYTE, &count);
    success = ::mxx::all_of(res == MPI_SUCCESS, subcomm);
		if (!success) {
      MPI_File_close(&fh);

      if (res != MPI_SUCCESS)
        throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string(filename, "write count", res, stat, subcomm));
		}

    success = ::mxx::all_of(count == curr_step, subcomm);
		if (!success) {
      MPI_File_close(&fh);

      if (count != curr_step) {
        std::stringstream ss;
        ss << "ERROR in mpiio: rank " << comm.rank() << " subcomm rank " << subcomm.rank() << " write error. request " << curr_step << " bytes got " << count << " bytes" << std::endl;

        throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
      }
		}

		global_offset += curr_step;
		data += curr_step;
		remainder -= curr_step;
	}

	// close the file when done.
	MPI_File_close(&fh);
	}

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


// DEPRECATED
//ListRankedChainNodeVecType to_compacted_chain(ChainGraphType const & chainmap,
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
//				std::cout << "rank " << comm.rank() << " printing " << std::endl << std::flush;

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

void print_chain_biedges(std::string const & filename,
		ChainGraphType const & chainmap,
		mxx::comm const & comm) {

	BL_BENCH_INIT(print_chain_biedges);

	//===  print chain nodes (1)
	if (comm.rank() == 0) printf("PRINT CHAINMAP Nodes\n");

	// print out.
	BL_BENCH_START(print_chain_biedges);
	{
		std::stringstream ss2;

		for (auto it = chainmap.cbegin(); it != chainmap.cend(); ++it) {
			ss2 << (*it) << std::endl;
		}
		write_mpiio(filename, ss2.str().c_str(), ss2.str().length(), comm);

	}
	//      std::cout << ss.str() << std::endl;
	BL_BENCH_COLLECTIVE_END(print_chain_biedges, "print chain biedges", chainmap.local_size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(print_chain_biedges, "print_chain_biedges", comm);

}

template <typename CHAIN_SUM>
void print_chain_summaries(std::string const & filename,
		std::vector<CHAIN_SUM> const & chains,
		mxx::comm const & comm) {

	BL_BENCH_INIT(print_chain_summaries);

	//===  print chain nodes (1)
	if (comm.rank() == 0) printf("PRINT CHAIN SUMMARIES\n");

	// print out.
	BL_BENCH_START(print_chain_summaries);
	{
		std::stringstream ss2;

		for (auto it = chains.cbegin(); it != chains.cend(); ++it) {
			ss2 << std::get<0>(*it) << " " << std::get<1>(*it) << " " << 
			std::get<2>(*it) << " " << std::get<3>(*it) << " " << 
			std::get<4>(*it) << " " << static_cast<size_t>(std::get<5>(*it)) << " " << static_cast<size_t>(std::get<6>(*it)) <<
			std::endl;
		}
		write_mpiio(filename, ss2.str().c_str(), ss2.str().length(), comm);

	}
	//      std::cout << ss.str() << std::endl;
	BL_BENCH_COLLECTIVE_END(print_chain_summaries, "print chain summaries", chains.size(), comm);

	BL_BENCH_REPORT_MPI_NAMED(print_chain_summaries, "", comm);

}




#endif // COMPACT_DBG_IO_CPP
