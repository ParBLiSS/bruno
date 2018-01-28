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
#ifndef COMPACT_DBG_BUILD_CPP
#define COMPACT_DBG_BUILD_CPP


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

	if (comm.rank() == 0) printf("PARSING and INSERT incrementally\n");

    using CharIterType = typename ::bliss::io::file_data::const_iterator;
	using SeqParserType = FileParser<CharIterType>;
	using SeqIterType = SplitSeqIterType<CharIterType, FileParser>;
	using KmerParser = ::bliss::debruijn::biedge::io::debruijn_kmer_simple_biedge_parser<KmerType>;
	using Iter = typename ::bliss::iterator::ContainerConcatenatingIterator<SeqIterType, KmerParser>;

	for (auto x : file_data) {

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
    	size_t block_size = (free_mem / (8 * sizeof(typename KmerParser::value_type)));  // number of elements that can be held in freemem
    	block_size = std::min(block_size, x.getRange().size());

    	if (comm.rank() == 0) std::cout << "estimate num elements=" << block_size << ", value_type size=" <<
    			sizeof(typename KmerParser::value_type) << " bytes" << std::endl;

		BL_BENCH_COLLECTIVE_END(build, "parse_setup", block_size, comm);

    	//=== copy into array incrementally
		BL_BENCH_START(build);
		idx.insert_incremental(start, endd, block_size);
		BL_BENCH_COLLECTIVE_END(build, "insert_incr", idx.local_size(), comm);

	}


	size_t total = idx.size();
	if (comm.rank() == 0) printf("PARSING and INSERT incremental DONE: total size after insert/rehash is %lu\n", total);

	BL_BENCH_REPORT_MPI_NAMED(build, "construct_incr", comm);
}
#endif




#endif // COMPACT_DBG_BUILD_CPP