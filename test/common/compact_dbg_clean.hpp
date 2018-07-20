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
 * compact_dbg_clean.cpp
 *
 * purpose:  routins for cleaning graph
 * 
 *  created: Jul 20, 2018
 *      Author: tony pan
 *
 */
#ifndef COMPACT_DBG_CLEAN_CPP
#define COMPACT_DBG_CLEAN_CPP

#include <mxx/comm.hpp>
#include "../common/compact_dbg_io.hpp"

//================= FILTER DEFINITION =====================

// spurious link is a chain node, with kmer frequency = left + right edges (result of deadends.)
// if overlap is 1 kmer, then above holds.  if overlap is 2 kmers or more, than we should see
//   kmer frequency / 2 > left or right.  detect this...
struct spurious_link_filter
{
	template <typename KM, typename Edges>
	bool operator()(::std::pair<KM, Edges> const &x) const
	{
		typename Edges::CountType l = 0, r = 0;
		for (int i = 0; i < Edges::maxEdgeCount; ++i)
		{
			l += x.second.get_in_edge_frequency(i);
			r += x.second.get_out_edge_frequency(i);
		}

		typename Edges::CountType me = x.second.get_self_frequency();
		return (me > (std::min(l, r) >> 1));
	}
};

// edge frequency filter is used to identify edges at some threshold.  place holder for relative frequency comparison
struct edge_freq_filter
{
	uint min_freq;
	edge_freq_filter(uint minf = 1) : min_freq(minf) {}

	template <typename count_type>
	bool operator()(count_type const &x) const
	{
		return x < min_freq;
	}
};

// deadend filter applies to chain summaries to select ones that are shorter than some length, k
struct deadend_filter
{
	uint max_len;
	edge_freq_filter pred;

	deadend_filter(uint length = std::numeric_limits<uint>::max(), uint minf = 1) : max_len(length), pred(minf) {}

	template <typename SUMMARY>
	bool operator()(SUMMARY const &x) const
	{
		return (std::get<4>(x) < max_len) &&
			   (pred(std::get<5>(x)) && // recall the length is already too small.
				pred(std::get<6>(x)));
	}
};

// bubble filter applies to chain summaries to select paths in bubble with similar lengths, at least k  (can't be smaller than k and still be a bubble)
// assume sorting by length. for now, just say difference is less than or equal to 1.
struct bubble_filter
{
	uint max_len;
	edge_freq_filter pred;

	bubble_filter(uint minf = 1) : pred(minf) {}

	// binary operator to operate on 2 paths
	template <typename SUMMARY>
	bool operator()(SUMMARY const &lhs, SUMMARY const &rhs) const
	{
		uint l = std::get<4>(lhs);
		uint r = std::get<4>(rhs);
		//printf("bubble lengths l %u r %u\n", l, r );
		return ((std::max(l, r) - std::min(l, r)) <= 1) &&
			   (pred(std::get<5>(lhs)) ||
				pred(std::get<6>(lhs)) ||
				pred(std::get<5>(rhs)) ||
				pred(std::get<6>(rhs)));
	}

	// unary operator to operate on 1 path - probably not the best for bubbles, so make it return true always.
	template <typename SUMMARY>
	bool operator()(SUMMARY const &x) const
	{
		return true;
		//				pred(std::get<5>(x)) ||
		//				pred(std::get<6>(x)) ;
	}
};


//================ LOGIC convenience function ====================

template <typename Index>
void remove_spurious_links(Index &idx, std::string const & out_prefix, 
	::mxx::comm const & comm, bool const & benchmark = true) {
	BL_BENCH_INIT(spurious);

	// predicate is not right yet...
	//==============================================================
	// spurious link.  if these are arising from deadend merging, then
	// we expect higher frequencies for the common kmers than the kmers closer to the branch
	// specifically, a sharp increase in frequency.
	// these can be filtered node by node before chain forming.
	spurious_link_filter spurious_link_filt;

#ifndef NDEBUG
	// find low freq
	BL_BENCH_START(spurious);
	auto spur_links = idx.get_map().find(spurious_link_filt);
	BL_BENCH_COLLECTIVE_END(spurious, "find_spur_links", spur_links.size(), comm);

	BL_BENCH_START(spurious);
	std::string spur_filename(out_prefix);
	spur_filename.append("_spurious.edges");
	print_graph_nodes(spur_filename, spur_links, comm);
	BL_BENCH_COLLECTIVE_END(spurious, "print_low_freq", spur_links.size(), comm);
#endif

	BL_BENCH_START(spurious);
	idx.get_map().erase_nodes(spurious_link_filt);
	idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
	BL_BENCH_COLLECTIVE_END(spurious, "erase_spurious", idx.local_size(), comm);

#ifndef NDEBUG
	BL_BENCH_START(spurious);
	std::string graph_filename(out_prefix);
	graph_filename.append(".graph.no_spurious.nodes");
	print_graph_edge_frequencies(graph_filename, idx, comm);
	BL_BENCH_COLLECTIVE_END(spurious, "print graph", idx.local_size(), comm);

	BL_BENCH_START(spurious);
	if (comm.rank() == 0) printf("rank 0 checking spurious links removed index\n");
	print_edge_histogram(idx, comm);
	check_index(idx, comm);
	BL_BENCH_COLLECTIVE_END(spurious, "histo", idx.local_size(), comm);
#endif

	BL_BENCH_REPORT_MPI_NAMED(spurious, "spurious_link", comm);
}




// clean recompact.
template <typename Index, typename ChainMap>
void clean_deadend_and_bubbles_recompact(Index & idx, ChainMap & chainmap, 
std::vector<size_t> const & threshes, ::std::string const & out_prefix,
::mxx::comm const & comm, bool const & benchmark = true) {
	BL_BENCH_INIT(clean_recompact);

	BL_BENCH_START(clean_recompact);
	edge_freq_filter edge_freq_filt(threshes[2] + 1);
	deadend_filter deadend_filt(2 * KmerType::size - 1, threshes[2] + 1);
	bubble_filter bubble_filt(threshes[2] + 1);

	//TODO: continue to verify that recompact is correct.

	ChainGraphType old_chains(comm);
	ChainGraphType new_chains(comm);
	chainmap.make_terminal_chain_graph(old_chains);
	BL_BENCH_COLLECTIVE_END(clean_recompact, "alloc", comm.size(), comm);
	
	// =============================================================
	// find deadends
	BL_BENCH_START(clean_recompact);
	auto deadends = ::bliss::debruijn::topology::find_deadends(idx, old_chains, deadend_filt);
	BL_BENCH_COLLECTIVE_END(clean_recompact, "find_deadends", deadends.size(), comm);

	// find bubbles
	BL_BENCH_START(clean_recompact);
	auto bubbles = ::bliss::debruijn::topology::find_bubbles(idx, old_chains, bubble_filt, comm);
	BL_BENCH_COLLECTIVE_END(clean_recompact, "find_bubbles", bubbles.size(), comm);

	if (!benchmark)
	{
		BL_BENCH_START(clean_recompact);
		std::string chain_deadend_filename(out_prefix);
		chain_deadend_filename.append(".chain.summary.deadend.");
		chain_deadend_filename.append(std::to_string(0));
		print_chain_summaries(chain_deadend_filename, deadends, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "print_deadend", deadends.size(), comm);

		BL_BENCH_START(clean_recompact);
		std::string chain_bubble_filename(out_prefix);
		chain_bubble_filename.append(".chain.summary.bubble.");
		chain_bubble_filename.append(std::to_string(0));
		print_chain_summaries(chain_bubble_filename, bubbles, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "print_bubbles", bubbles.size(), comm);
	}

	bool done_clean = (deadends.size() == 0) && (bubbles.size() == 0);
	done_clean = mxx::all_of(done_clean, comm);

	size_t iteration = 0;
	std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge>> branch_nodes;
	::bliss::debruijn::biedge::compact_simple_biedge r, l;
	std::vector<KmerType> modified;

	while (!done_clean)
	{
		if (comm.rank() == 0)
			printf("cleaning iteration %lu\n", iteration);

		modified.clear();

		BL_BENCH_START(clean_recompact);
		//---------- remove deadends.
		// extract the edges as nodes.
		branch_nodes.clear();
		for (size_t i = 0; i < deadends.size(); ++i)
		{
			if (std::get<5>(deadends[i]) > 0)
			{ //frequency > 0 -> has edge to branch
				r.setCharsAtPos(
					::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<1>(deadends[i]).getCharsAtPos(0, 1)]],
					0, 1);
				branch_nodes.emplace_back(std::get<0>(deadends[i]), r);

				modified.emplace_back(std::get<0>(deadends[i]));
				modified.emplace_back(std::get<1>(deadends[i]));
			}
			if (std::get<6>(deadends[i]) > 0)
			{ //frequency > 0 -> has edge to branch
				l.setCharsAtPos(
					::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<2>(deadends[i]).getCharsAtPos(KmerType::size - 1, 1)]],
					1, 1);
				branch_nodes.emplace_back(std::get<3>(deadends[i]), l);

				modified.emplace_back(std::get<2>(deadends[i]));
				modified.emplace_back(std::get<3>(deadends[i]));
			}
		}
		BL_BENCH_COLLECTIVE_END(clean_recompact, "make_deadend_list", branch_nodes.size(), comm);

#ifndef NDEBUG
		{
			BL_BENCH_START(clean_recompact);
			auto branch_n = idx.get_map().find_edges(branch_nodes);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "find_deadend_branches", branch_n.size(), comm);

			BL_BENCH_START(clean_recompact);
			std::string edge_filter_filename(out_prefix);
			edge_filter_filename.append(".branch.summary.deadend.");
			edge_filter_filename.append(std::to_string(iteration));
			print_graph_edges(edge_filter_filename, branch_n, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_deadend_branches", branch_n.size(), comm);
		}
#endif

		// remove the deadend edges with freq smaller than some threshold
		BL_BENCH_START(clean_recompact);
		idx.get_map().erase_edges(branch_nodes, edge_freq_filt); // also has to meet edge frequency requirements.
		BL_BENCH_COLLECTIVE_END(clean_recompact, "severe_deadends", branch_nodes.size(), comm);

#ifndef NDEBUG
		BL_BENCH_START(clean_recompact);
		{
			std::string graph_filename(out_prefix);
			graph_filename.append(".graph.no_deadend.");
			graph_filename.append(std::to_string(iteration));
			print_graph_edge_frequencies(graph_filename, idx, comm);
		}
		BL_BENCH_COLLECTIVE_END(clean_recompact, "print_graph", idx.local_size(), comm);

		BL_BENCH_START(clean_recompact);
		if (comm.rank() == 0)
			printf("rank 0 checking deadend-removed index\n");
		print_edge_histogram(idx, comm);
		// check_index(idx, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "histo", idx.local_size(), comm);
#endif

		//---------------------------
		// remove bubbles

		// extract the edges as nodes.
		BL_BENCH_START(clean_recompact);
		branch_nodes.clear();
		for (size_t i = 0; i < bubbles.size(); ++i)
		{
			r.setCharsAtPos(
				::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<1>(bubbles[i]).getCharsAtPos(0, 1)]],
				0, 1);
			branch_nodes.emplace_back(std::get<0>(bubbles[i]), r);
			modified.emplace_back(std::get<0>(bubbles[i]));
			modified.emplace_back(std::get<1>(bubbles[i]));

			l.setCharsAtPos(
				::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<2>(bubbles[i]).getCharsAtPos(KmerType::size - 1, 1)]],
				1, 1);
			branch_nodes.emplace_back(std::get<3>(bubbles[i]), l);
			modified.emplace_back(std::get<2>(bubbles[i]));
			modified.emplace_back(std::get<3>(bubbles[i]));
		}
		BL_BENCH_COLLECTIVE_END(clean_recompact, "make_bubble_list", branch_nodes.size(), comm);

#ifndef NDEBUG
		{
			BL_BENCH_START(clean_recompact);
			auto branch_n = idx.get_map().find_edges(branch_nodes);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "find_bubble_branches", branch_n.size(), comm);

			BL_BENCH_START(clean_recompact);
			std::string edge_filter_filename(out_prefix);
			edge_filter_filename.append(".branch.summary.bubble.");
			edge_filter_filename.append(std::to_string(iteration));
			print_graph_edges(edge_filter_filename, branch_n, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_bubble_branches", branch_n.size(), comm);
		}
#endif

		BL_BENCH_START(clean_recompact);
		idx.get_map().erase_edges(branch_nodes, edge_freq_filt);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "severe_bubbles", branch_nodes.size(), comm);

#ifndef NDEBUG
		BL_BENCH_START(clean_recompact);
		{
			std::string graph_filename(out_prefix);
			graph_filename.append(".graph.no_bubble.");
			graph_filename.append(std::to_string(iteration));
			print_graph_edge_frequencies(graph_filename, idx, comm);
		}
		BL_BENCH_COLLECTIVE_END(clean_recompact, "print graph", idx.local_size(), comm);

		BL_BENCH_START(clean_recompact);
		if (comm.rank() == 0)
			printf("rank 0 checking bubble removed index\n");
		print_edge_histogram(idx, comm);
		// check_index(idx, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "histo", idx.local_size(), comm);

		old_chains.print_stats("old");
#endif

		//---------- recompact  old and put inot new.
		BL_BENCH_START(clean_recompact);
		::bliss::debruijn::topology::recompact(idx, modified, old_chains, new_chains, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "recompact", new_chains.size(), comm);

#ifndef NDEBUG
		new_chains.print_stats("new_compact");
		{

			BL_BENCH_START(clean_recompact);
			std::string chain_biedge_filename(out_prefix);
			chain_biedge_filename.append(".debug.chainmap.termini.");
			chain_biedge_filename.append(std::to_string(iteration));
			print_chain_biedges(chain_biedge_filename, new_chains, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_chain_biedge", new_chains.local_size(), comm);
		}
#endif
		if (!benchmark)
		{
			// =============================================================
			// generate chain_summaries
			BL_BENCH_START(clean_recompact);
			auto summaries = new_chains.to_summarized_chains();
			BL_BENCH_COLLECTIVE_END(clean_recompact, "chain_summaries", summaries.size(), comm);

			BL_BENCH_START(clean_recompact);
			std::string chain_summary_filename(out_prefix);
			chain_summary_filename.append(".chain.summary.");
			chain_summary_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_summary_filename, summaries, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_chain_summaries", summaries.size(), comm);
		}

		//------------ recompact_finalize from old to new.
		BL_BENCH_START(clean_recompact);
		::bliss::debruijn::topology::recompact_finalize(old_chains, new_chains, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "recompact_finalize", old_chains.size(), comm);

#ifndef NDEBUG
		old_chains.print_stats("old_merged");
		{

			BL_BENCH_START(clean_recompact);
			std::string chain_biedge_filename(out_prefix);
			chain_biedge_filename.append(".debug.chainmap.finalized.");
			chain_biedge_filename.append(std::to_string(iteration));
			print_chain_biedges(chain_biedge_filename, old_chains, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_chain_biedge", old_chains.local_size(), comm);
		}
#endif
		if (!benchmark)
		{
			// =============================================================
			// generate chain_summaries
			BL_BENCH_START(clean_recompact);
			auto summaries = old_chains.to_summarized_chains();
			BL_BENCH_COLLECTIVE_END(clean_recompact, "chain_summaries", summaries.size(), comm);

			BL_BENCH_START(clean_recompact);
			std::string chain_summary_filename(out_prefix);
			chain_summary_filename.append(".chain.summary.finalized.");
			chain_summary_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_summary_filename, summaries, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_chain_summaries", summaries.size(), comm);
		}

		{
			// =========== remove isolated - don't remove cycles because we need to use the info to update chains.
			BL_BENCH_START(clean_recompact);
			{
				size_t before = idx.local_size();
				idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
				printf("rank %d ERASE ISOLATED %lu after cycle 1 iter %ld\n", comm.rank(), idx.local_size() - before, iteration);
			}
			BL_BENCH_COLLECTIVE_END(clean_recompact, "remove isolated", idx.local_size(), comm);

#ifndef NDEBUG

			BL_BENCH_START(clean_recompact);
			if (comm.rank() == 0)
				printf("rank 0 checking isolated removed index\n");
			print_edge_histogram(idx, comm);
			check_index(idx, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "histo", idx.local_size(), comm);
#endif
		}

		BL_BENCH_START(clean_recompact);
		// before making terminal_chain_graph, need to clean up new_chains of the cycles and isolated items.
		// else incorrect chains could be introduced.
		new_chains.clear();
		BL_BENCH_COLLECTIVE_END(clean_recompact, "setup_new_chains", old_chains.size(), comm);

#ifndef NDEBUG
		{

			BL_BENCH_START(clean_recompact);
			std::string chain_biedge_filename(out_prefix);
			chain_biedge_filename.append("debug.termini.nodes.");
			chain_biedge_filename.append(std::to_string(iteration));
			print_chain_biedges(chain_biedge_filename, old_chains, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_chain_biedge", old_chains.local_size(), comm);
		}
#endif
		if (!benchmark)
		{

			BL_BENCH_START(clean_recompact);
			auto summaries = old_chains.to_summarized_chains();
			BL_BENCH_COLLECTIVE_END(clean_recompact, "chain_summaries", summaries.size(), comm);

			BL_BENCH_START(clean_recompact);
			std::string chain_summary_filename(out_prefix);
			chain_summary_filename.append(".chain.summary.clean.");
			chain_summary_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_summary_filename, summaries, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_chain_summaries", summaries.size(), comm);
		}

		//---------- detect deadends and bubbles again.
		BL_BENCH_START(clean_recompact);
		deadends = ::bliss::debruijn::topology::find_deadends(idx, old_chains, deadend_filt);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "find_deadends", deadends.size(), comm);

		// find bubbles
		BL_BENCH_START(clean_recompact);
		bubbles = ::bliss::debruijn::topology::find_bubbles(idx, old_chains, bubble_filt, comm);
		BL_BENCH_COLLECTIVE_END(clean_recompact, "find_bubbles", bubbles.size(), comm);

		if (!benchmark)
		{
			BL_BENCH_START(clean_recompact);
			std::string chain_deadend_filename(out_prefix);
			chain_deadend_filename.append(".chain.summary.deadend.");
			chain_deadend_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_deadend_filename, deadends, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_deadend", deadends.size(), comm);

			BL_BENCH_START(clean_recompact);
			std::string chain_bubble_filename(out_prefix);
			chain_bubble_filename.append(".chain.summary.bubble.");
			chain_bubble_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_bubble_filename, bubbles, comm);
			BL_BENCH_COLLECTIVE_END(clean_recompact, "print_bubbles", bubbles.size(), comm);
		}

		done_clean = (deadends.size() == 0) && (bubbles.size() == 0);
		done_clean = mxx::all_of(done_clean, comm);

		++iteration;
	}  // while loop to repeatedly clean bubbles and deadends

	// don't erase the isolated - else the freq map computation would be incorrect.  (fix?  now removing entries from chainmap too.)

	// ----- use recompact_finalize to update complete chainmap
	BL_BENCH_START(clean_recompact);
	old_chains.unseparate_isolated_and_cycles();
	::bliss::debruijn::topology::recompact_finalize(chainmap, old_chains, comm);
	BL_BENCH_COLLECTIVE_END(clean_recompact, "recompact_finalize", chainmap.size(), comm);

#ifndef NDEBUG
	chainmap.print_stats("final_merge");
#endif

	// before making terminal_chain_graph, need to clean up new_chains of the cycles and isolated items.
	// else incorrect chains could be introduced.
	// done in recompact_finalize

	// remove the cycles and isolated.
	BL_BENCH_START(clean_recompact);
	auto cycle_kmers = chainmap.get_cycle_node_kmers();
	idx.get_map().erase_nodes(cycle_kmers);
	{
		size_t before = idx.local_size();
		idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
		printf("rank %d ERASE ISOLATED %lu after cycle 1 iter %ld\n", comm.rank(), idx.local_size() - before, iteration);
	}
	BL_BENCH_COLLECTIVE_END(clean_recompact, "remove cycles/isolated/etc", idx.local_size(), comm);


#ifndef NDEBUG
	BL_BENCH_START(clean_recompact);
	if (comm.rank() == 0)
		printf("rank 0 checking cleaned index\n");
	print_edge_histogram(idx, comm);
	check_index(idx, comm);
	BL_BENCH_COLLECTIVE_END(clean_recompact, "histo", idx.local_size(), comm);
#endif

	BL_BENCH_REPORT_MPI_NAMED(clean_recompact, "clean_recompact", comm);
}



// clean naive list rank.

template <typename Index, typename ChainMap>
void clean_deadend_and_bubbles_naive(Index &idx, ChainMap &chainmap, 
	::std::vector<size_t> const & threshes, ::std::string const & out_prefix,
	::mxx::comm const & comm, bool const & LRoptimized = true, bool const & benchmark = true) {
	BL_BENCH_INIT(clean_compact);

	BL_BENCH_START(clean_compact);

	edge_freq_filter edge_freq_filt(threshes[2] + 1);
	deadend_filter deadend_filt(2 * KmerType::size - 1, threshes[2] + 1);
	bubble_filter bubble_filt(threshes[2] + 1);

	ChainGraphType new_chains(comm);
	chainmap.make_terminal_chain_graph(new_chains);
	BL_BENCH_COLLECTIVE_END(clean_compact, "alloc", comm.size(), comm);

	// =============================================================
	// find deadends
	BL_BENCH_START(clean_compact);
	auto deadends = ::bliss::debruijn::topology::find_deadends(idx, new_chains, deadend_filt);
	BL_BENCH_COLLECTIVE_END(clean_compact, "find_deadends", deadends.size(), comm);

	// find bubbles
	BL_BENCH_START(clean_compact);
	auto bubbles = ::bliss::debruijn::topology::find_bubbles(idx, new_chains, bubble_filt, comm);
	BL_BENCH_COLLECTIVE_END(clean_compact, "find_bubbles", bubbles.size(), comm);

	if (!benchmark)
	{
		BL_BENCH_START(clean_compact);
		std::string chain_deadend_filename(out_prefix);
		chain_deadend_filename.append(".chain.summary.deadend.");
		chain_deadend_filename.append(std::to_string(0));
		print_chain_summaries(chain_deadend_filename, deadends, comm);
		BL_BENCH_COLLECTIVE_END(clean_compact, "print_deadend", deadends.size(), comm);

		BL_BENCH_START(clean_compact);
		std::string chain_bubble_filename(out_prefix);
		chain_bubble_filename.append(".chain.summary.bubble.");
		chain_bubble_filename.append(std::to_string(0));
		print_chain_summaries(chain_bubble_filename, bubbles, comm);
		BL_BENCH_COLLECTIVE_END(clean_compact, "print_bubbles", bubbles.size(), comm);
	}

	bool done_clean = (deadends.size() == 0) && (bubbles.size() == 0);
	done_clean = mxx::all_of(done_clean, comm);

	size_t iteration = 0;
	std::vector<std::pair<KmerType, ::bliss::debruijn::biedge::compact_simple_biedge>> branch_nodes;
	::bliss::debruijn::biedge::compact_simple_biedge r, l;

	while (!done_clean)
	{
		if (comm.rank() == 0)
			printf("cleaning iteration %lu\n", iteration);

		BL_BENCH_START(clean_compact);
		//---------- remove deadends.
		// extract the edges as nodes.
		branch_nodes.clear();
		for (size_t i = 0; i < deadends.size(); ++i)
		{
			if (std::get<5>(deadends[i]) > 0)
			{ // frequency > 0 -> has edge to branch
				r.setCharsAtPos(
					::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<1>(deadends[i]).getCharsAtPos(0, 1)]],
					0, 1);
				branch_nodes.emplace_back(std::get<0>(deadends[i]), r);
			}
			if (std::get<6>(deadends[i]) > 0)
			{ // frequency > 0 -> has edge to branch
				l.setCharsAtPos(
					::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<2>(deadends[i]).getCharsAtPos(KmerType::size - 1, 1)]],
					1, 1);
				branch_nodes.emplace_back(std::get<3>(deadends[i]), l);
			}
		}
		BL_BENCH_COLLECTIVE_END(clean_compact, "make_deadend_list", branch_nodes.size(), comm);

#ifndef NDEBUG
		{
			BL_BENCH_START(clean_compact);
			auto branch_n = idx.get_map().find_edges(branch_nodes);
			BL_BENCH_COLLECTIVE_END(clean_compact, "find_deadend_branches", branch_n.size(), comm);

			BL_BENCH_START(clean_compact);
			std::string edge_filter_filename(out_prefix);
			edge_filter_filename.append(".branch.summary.deadend.");
			edge_filter_filename.append(std::to_string(iteration));
			print_graph_edges(edge_filter_filename, branch_n, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_deadend_branches", branch_n.size(), comm);
		}
#endif

		// remove the deadend edges with freq smaller than some threshold
		BL_BENCH_START(clean_compact);
		idx.get_map().erase_edges(branch_nodes, edge_freq_filt); // also has to meet edge frequency requirements.
		{
			size_t before = idx.local_size();
			idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
			printf("rank %d ERASE ISOLATED %lu after deadends iter %ld\n", comm.rank(), idx.local_size() - before, iteration);
		}
		BL_BENCH_COLLECTIVE_END(clean_compact, "severe_deadends", branch_nodes.size(), comm);

#ifndef NDEBUG

		BL_BENCH_START(clean_compact);
		{
			std::string graph_filename(out_prefix);
			graph_filename.append(".graph.no_deadend.");
			graph_filename.append(std::to_string(iteration));
			print_graph_edge_frequencies(graph_filename, idx, comm);
		}
		BL_BENCH_COLLECTIVE_END(clean_compact, "print_graph", idx.local_size(), comm);

		BL_BENCH_START(clean_compact);
		if (comm.rank() == 0)
			printf("rank 0 checking deadend-removed index\n");
		print_edge_histogram(idx, comm);
		check_index(idx, comm);
		BL_BENCH_COLLECTIVE_END(clean_compact, "histo", idx.local_size(), comm);
#endif

		//---------------------------
		// remove bubbles

		// extract the edges as nodes.
		BL_BENCH_START(clean_compact);
		branch_nodes.clear();
		for (size_t i = 0; i < bubbles.size(); ++i)
		{
			r.setCharsAtPos(
				::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<1>(bubbles[i]).getCharsAtPos(0, 1)]],
				0, 1);
			branch_nodes.emplace_back(std::get<0>(bubbles[i]), r);

			l.setCharsAtPos(
				::bliss::common::DNA16::FROM_ASCII[KmerType::KmerAlphabet::TO_ASCII[std::get<2>(bubbles[i]).getCharsAtPos(KmerType::size - 1, 1)]],
				1, 1);
			branch_nodes.emplace_back(std::get<3>(bubbles[i]), l);
		}
		BL_BENCH_COLLECTIVE_END(clean_compact, "make_bubble_list", branch_nodes.size(), comm);

#ifndef NDEBUG
		{
			BL_BENCH_START(clean_compact);
			auto branch_n = idx.get_map().find_edges(branch_nodes);
			BL_BENCH_COLLECTIVE_END(clean_compact, "find_bubble_branches", branch_n.size(), comm);

			BL_BENCH_START(clean_compact);
			std::string edge_filter_filename(out_prefix);
			edge_filter_filename.append(".branch.summary.bubble.");
			edge_filter_filename.append(std::to_string(iteration));
			print_graph_edges(edge_filter_filename, branch_n, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_bubble_branches", branch_n.size(), comm);
		}
#endif

		BL_BENCH_START(clean_compact);
		idx.get_map().erase_edges(branch_nodes, edge_freq_filt);
		{
			size_t before = idx.local_size();
			idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
			printf("rank %d ERASE ISOLATED %lu after bubble iter %ld\n", comm.rank(), idx.local_size() - before, iteration);
		}
		BL_BENCH_COLLECTIVE_END(clean_compact, "severe_bubbles", branch_nodes.size(), comm);

#ifndef NDEBUG
		BL_BENCH_START(clean_compact);
		{
			std::string graph_filename(out_prefix);
			graph_filename.append(".graph.no_bubble.");
			graph_filename.append(std::to_string(iteration));
			print_graph_edge_frequencies(graph_filename, idx, comm);
		}
		BL_BENCH_COLLECTIVE_END(clean_compact, "print graph", idx.local_size(), comm);

		BL_BENCH_START(clean_compact);
		if (comm.rank() == 0)
			printf("rank 0 checking bubble removed index\n");
		print_edge_histogram(idx, comm);
		check_index(idx, comm);
		BL_BENCH_COLLECTIVE_END(clean_compact, "histo", idx.local_size(), comm);
#endif

		//---------- recompact
		// clear the chain map.

		// ==== make chain map
		BL_BENCH_START(clean_compact);
		chainmap.clear();
		chainmap.extract_chains(idx);
		//make_chain_map(idx, chainmap, comm);
		BL_BENCH_COLLECTIVE_END(clean_compact, "chainmap", chainmap.local_size(), comm);
		// == DONE == make chain map
#ifndef NDEBUG
		chainmap.print_stats("chainmap_reset");
#endif
		// ===== parallel list ranking for chain compaction
		{
			BL_BENCH_START(clean_compact);
			size_t iters = 0;
			if (LRoptimized)
				iters = chainmap.list_rank_min_update();
			else
				iters = chainmap.list_rank();
			//auto cycle_node_kmers = list_rank(chainmap, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "list_rank", iters, comm);
		} // == DONE == parallel list ranking for chain compaction

#ifndef NDEBUG
		chainmap.print_stats("chainmap_compact");
		{

			BL_BENCH_START(clean_compact);
			std::string chain_biedge_filename(out_prefix);
			chain_biedge_filename.append(".debug.chainmap.finalized.");
			chain_biedge_filename.append(std::to_string(iteration));
			print_chain_biedges(chain_biedge_filename, chainmap, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_chain_biedge", chainmap.local_size(), comm);
		}
#endif
		if (!benchmark)
		{
			// =============================================================
			// generate chain_summaries
			BL_BENCH_START(clean_compact);
			auto summaries = chainmap.to_summarized_chains();
			BL_BENCH_COLLECTIVE_END(clean_compact, "chain_summaries", summaries.size(), comm);

			BL_BENCH_START(clean_compact);
			std::string chain_summary_filename(out_prefix);
			chain_summary_filename.append(".chain.summary.finalized.");
			chain_summary_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_summary_filename, summaries, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_chain_summaries", summaries.size(), comm);
		}

		{
			// =========== remove cycles and isolated
			BL_BENCH_START(clean_compact);
			auto cycle_kmers = chainmap.get_cycle_node_kmers();
			idx.get_map().erase_nodes(cycle_kmers);
			{
				size_t before = idx.local_size();
				idx.erase_if(::bliss::debruijn::filter::graph::IsIsolated());
				printf("rank %d ERASE ISOLATED %lu after cycles iter %ld\n", comm.rank(), idx.local_size() - before, iteration);
			}
			BL_BENCH_COLLECTIVE_END(clean_compact, "remove cycles/isolated/etc", idx.local_size(), comm);

#ifndef NDEBUG
			BL_BENCH_START(clean_compact);
			std::string graph_filename(out_prefix);
			graph_filename.append(".graph.no_cycle.");
			graph_filename.append(std::to_string(iteration));
			print_graph_edge_frequencies(graph_filename, idx, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print graph", idx.local_size(), comm);

			BL_BENCH_START(clean_compact);
			if (comm.rank() == 0)
				printf("rank 0 checking cycle removed index\n");
			print_edge_histogram(idx, comm);
			check_index(idx, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "histo", idx.local_size(), comm);
#endif
		}

		BL_BENCH_START(clean_compact);
		new_chains.clear();
		chainmap.make_terminal_chain_graph(new_chains);
		BL_BENCH_COLLECTIVE_END(clean_compact, "get_new_chains", new_chains.size(), comm);

		// generate chain_summaries

		if (!benchmark)
		{
			BL_BENCH_START(clean_compact);
			auto summaries = new_chains.to_summarized_chains();
			BL_BENCH_COLLECTIVE_END(clean_compact, "chain_summaries", summaries.size(), comm);

			BL_BENCH_START(clean_compact);
			std::string chain_summary_filename(out_prefix);
			chain_summary_filename.append(".chain.summary.clean.");
			chain_summary_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_summary_filename, summaries, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_chain_summaries", summaries.size(), comm);
		}

		//---------- detect deadends and bubbles again.
		BL_BENCH_START(clean_compact);
		deadends = ::bliss::debruijn::topology::find_deadends(idx, new_chains, deadend_filt);
		BL_BENCH_COLLECTIVE_END(clean_compact, "find_deadends", deadends.size(), comm);

		// find bubbles
		BL_BENCH_START(clean_compact);
		bubbles = ::bliss::debruijn::topology::find_bubbles(idx, new_chains, bubble_filt, comm);
		BL_BENCH_COLLECTIVE_END(clean_compact, "find_bubbles", bubbles.size(), comm);

		if (!benchmark)
		{
			BL_BENCH_START(clean_compact);
			std::string chain_deadend_filename(out_prefix);
			chain_deadend_filename.append(".chain.summary.deadend.");
			chain_deadend_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_deadend_filename, deadends, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_deadend", deadends.size(), comm);

			BL_BENCH_START(clean_compact);
			std::string chain_bubble_filename(out_prefix);
			chain_bubble_filename.append(".chain.summary.bubble.");
			chain_bubble_filename.append(std::to_string(iteration));
			print_chain_summaries(chain_bubble_filename, bubbles, comm);
			BL_BENCH_COLLECTIVE_END(clean_compact, "print_bubbles", bubbles.size(), comm);
		}

		done_clean = (deadends.size() == 0) && (bubbles.size() == 0);
		done_clean = mxx::all_of(done_clean, comm);

		++iteration;
	}  // while loop to repeatedly clean bubbles and deadends

#ifndef NDEBUG
	BL_BENCH_START(clean_compact);
	if (comm.rank() == 0)
		printf("rank 0 checking cleaned index\n");
	print_edge_histogram(idx, comm);
	check_index(idx, comm);
	BL_BENCH_COLLECTIVE_END(clean_compact, "histo", idx.local_size(), comm);
#endif

	BL_BENCH_REPORT_MPI_NAMED(clean_compact, "clean_compact", comm);
}


#endif // COMPACT_DBG_CLEAN_CPP