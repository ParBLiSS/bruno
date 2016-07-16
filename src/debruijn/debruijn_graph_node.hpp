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
 * debruijn_graph_node.hpp
 *
 *  Created on: Jul 31, 2015
 *      Author: yongchao
 *
 *  Rewrote: June 22, 2016
 *      Author: tony pan
 *
 */

#ifndef DEBRUIJN_GRAPH_NODE_HPP_
#define DEBRUIJN_GRAPH_NODE_HPP_

#include "bliss-config.hpp"

#include <vector>
#include <utility>      // pair and utility functions.
#include <type_traits>

#include "utils/logging.h"

#include "common/alphabets.hpp"
#include "utils/bit_ops.hpp"

#include "debruijn/edge_iterator.hpp"


// forward declares
namespace bliss
{
  namespace debruijn
  {
    namespace graph
    {
      template <typename EdgeEncoding, typename COUNT = bool, typename DUMMY = void>
      class compact_multi_biedge;
    }
  }
}

#if defined(USE_MPI)

#include <mxx/datatypes.hpp>

namespace mxx {

	template<typename A, typename T, typename D>
	  struct datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<A, T, D> >;

	template<typename A, typename T, typename D>
	  struct datatype_builder<const ::bliss::debruijn::graph::compact_multi_biedge<A, T, D> >;

} // mxx

#endif

namespace bliss
{
  namespace debruijn
  {
    namespace graph
    {

      /**
       * @brief kmer metadata holding the incoming and outgoing edge counts in a de bruijn graph.
       * @details lower 4 integers are out edge counts, upper 4 integers are in edge counts.
       *        the out and in edges are relative to the kmer's orientation.
       * @tparam EdgeEncoding  alphabet used for edge encoding.
       */
      template<typename COUNT, typename DUMMY>
      class compact_multi_biedge<::bliss::common::DNA, COUNT, DUMMY> {

#if defined(USE_MPI)
    	  friend class ::mxx::datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<::bliss::common::DNA, COUNT, DUMMY> >;
#endif

          static_assert(!::std::is_signed<COUNT>::value &&
                        ::std::is_integral<COUNT>::value, "only supports unsigned integer types for count");


        public:
          using EdgeEncoding = bliss::common::DNA;
          static constexpr size_t maxEdgeCount = 4;

        protected:
          COUNT sat_add(COUNT const & a, COUNT const & b) {
            COUNT c = a + b;
            return (c < a) ? -1 : c;
          }

          void sat_incr(COUNT & target) {
            target = sat_add(target, 1);
          }


          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{     //DNA16    // bitvec
                 8,   // 0000 .     // 0x00000000
                 0,   // 0001 A     // 0x00000001
                 1,   // 0010 C     // 0x00000010
                 8,   // 0011       // 0x00000000
                 2,   // 0100 G     // 0x00000100
                 8,   // 0101       // 0x00000000
                 8,   // 0110       // 0x00000000
                 8,   // 0111       // 0x00000000
                 3,   // 1000 T     // 0x00001000
                 8,   // 1001       // 0x00000000
                 8,   // 1010       // 0x00000000
                 8,   // 1011       // 0x00000000
                 8,   // 1100       // 0x00000000
                 8,   // 1101       // 0x00000000
                 8,   // 1110       // 0x00000000
                 8    // 1111 N     // 0x00000000
          }};

          static constexpr std::array<uint8_t, 4> CHAR_TO_INDEX =
          {{     //DNA4    // bitvec
                 0,   // 00 A     // 0x00000001
                 1,   // 01 C     // 0x00000010
                 2,   // 10 G     // 0x00000100
                 3    // 11 T     // 0x00001000
          }};

          static constexpr std::array<uint8_t, 4> INDEX_TO_CHAR =
          {{     //DNA4    // bitvec
                 0,   // 00 A     // 0x00000001
                 1,   // 01 C     // 0x00000010
                 2,   // 10 G     // 0x00000100
                 3    // 11 T     // 0x00001000
          }};

        public:

          /// array of counts.  format:  [out A C G T; in A C G T], ordered for the canonical strand, not necessarily same as for the input kmer..
          std::array<COUNT, 2 * maxEdgeCount + 1> counts;  // public because of friending mxx datatype_builder does not work

          using Alphabet = EdgeEncoding;
          using CountType = COUNT;
          using EdgeInputType = bliss::debruijn::compact_simple_biedge;   // hardcoded to DNA16 because of N


          friend std::ostream& operator<<(std::ostream& ost, const compact_multi_biedge<EdgeEncoding, COUNT, DUMMY> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: counts in = [";
            for (size_t i = maxEdgeCount; i < 2*maxEdgeCount; ++i) ost << node.counts[i] << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << node.counts[i] << ",";
            ost << "], self = " << node.counts[2 * maxEdgeCount];
            return ost;
          }


          /*constructor*/
          compact_multi_biedge() { counts.fill(0); };

          /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
          ~compact_multi_biedge() {}


          /**
           * @brief increments the edge counts.  does not perform flipping.
           * @param exts   input edges, should be at most 1 for in and out.
           */
          inline void update(EdgeInputType edges)
          {
            // take care of out
            uint8_t out = edges.getData()[0] & 0xF;
            if (FROM_DNA16[out] < maxEdgeCount) sat_incr(counts[FROM_DNA16[out]]);

            // take care of in.
            uint8_t in = edges.getData()[0] >> 4;
            if (FROM_DNA16[in] < maxEdgeCount) sat_incr(counts[FROM_DNA16[in] + maxEdgeCount]);

            sat_incr(counts[2 * maxEdgeCount]);
          }


          inline void merge(compact_multi_biedge const & other) {
            for (int i = 0; i <= 2 * maxEdgeCount; ++i) {
              counts[i] = sat_add(counts[i], other.counts[i]);
            }
          }

          inline COUNT get_out_edge_frequency(uint8_t idx) const {
            if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return counts[CHAR_TO_INDEX[idx]];
          }

          inline COUNT get_in_edge_frequency(uint8_t idx) const {
            if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return counts[CHAR_TO_INDEX[idx] + maxEdgeCount];
          }



          inline COUNT get_self_frequency() const {
        	  return counts[2 * maxEdgeCount];
          }

          inline uint8_t get_out_edge_count() const {
            uint8_t count = 0;
            for (size_t i = 0; i < maxEdgeCount; ++i ) {
              count += (counts[i] > 0 ? 1 : 0);
            }
            return count;
          }

          inline uint8_t get_in_edge_count() const {
            uint8_t count = 0;
            for (size_t i = maxEdgeCount; i < 2 * maxEdgeCount; ++i ) {
              count += counts[i] > 0 ? 1 : 0;
            }
            return count;
          }


          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_out_neighbors(Kmer const & kmer, std::vector<std::pair<Kmer, CountType> > & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            CountType count;
            for (unsigned char i = 0; i < maxEdgeCount; ++i) {   // no gap character
              count = counts[i];
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_in_neighbors(Kmer const & kmer, std::vector<std::pair<Kmer, CountType> > & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            CountType count;
            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // no gap character
              count = counts[i + maxEdgeCount];
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextReverseFromChar(i);
              }
            }
          }
      };

      template<typename COUNT, typename DUMMY>
      class compact_multi_biedge<::bliss::common::DNA6, COUNT, DUMMY> {

          static_assert(!::std::is_signed<COUNT>::value &&
                        ::std::is_integral<COUNT>::value, "only supports unsigned integer types for count");

#if defined(USE_MPI)
          friend class ::mxx::datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<::bliss::common::DNA6, COUNT, DUMMY> >;
#endif

        public:
          using EdgeEncoding = bliss::common::DNA6;
          static constexpr size_t maxEdgeCount = 5;

        protected:
          COUNT sat_add(COUNT const & a, COUNT const & b) {
            COUNT c = a + b;
            return (c < a) ? -1 : c;
          }

          void sat_incr(COUNT & target) {
            target = sat_add(target, 1);
          }


          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{     //DNA16    // bitvec
                 8,   // 0000 .     // 0x00000000
                 0,   // 0001 A     // 0x00000001
                 1,   // 0010 C     // 0x00000010
                 8,   // 0011       // 0x00000000
                 2,   // 0100 G     // 0x00000100
                 8,   // 0101       // 0x00000000
                 8,   // 0110       // 0x00000000
                 8,   // 0111       // 0x00000000
                 3,   // 1000 T     // 0x00001000
                 8,   // 1001       // 0x00000000
                 8,   // 1010       // 0x00000000
                 8,   // 1011       // 0x00000000
                 8,   // 1100       // 0x00000000
                 8,   // 1101       // 0x00000000
                 8,   // 1110       // 0x00000000
                 4    // 1111 N     // 0x00000000
          }};

          static constexpr std::array<uint8_t, 8> CHAR_TO_INDEX =
          {{     //DNA16    // bitvec
                 8,   // 000 .     // 0x00000000
                 0,   // 001 A     // 0x00000001
                 8,   // 010       // 0x00000010
                 1,   // 011 C     // 0x00000000
                 3,   // 100 T     // 0x00000100
                 8,   // 101       // 0x00000000
                 2,   // 110 G     // 0x00000000
                 4    // 111 N     // 0x00000000
          }};

          static constexpr std::array<uint8_t, 5> INDEX_TO_CHAR =
          {{     //DNA6    // bitvec
                 1,   // 001 A     // 0x00000001
                 3,   // 011 C     // 0x00000010
                 6,   // 110 G     // 0x00000100
                 4,   // 100 T     // 0x00001000
				 7,   // 111 N
          }};


        public:

          /// array of counts.  format:  [out A C G T; in A C G T], ordered for the canonical strand, not necessarily same as for the input kmer..
          std::array<COUNT, 2 * maxEdgeCount + 1> counts;  // public because of friending mxx datatype_builder does not work

          using Alphabet = EdgeEncoding;
          using CountType = COUNT;
          using EdgeInputType = bliss::debruijn::compact_simple_biedge;   // hardcoded to DNA16 because of N


          friend std::ostream& operator<<(std::ostream& ost, const compact_multi_biedge<EdgeEncoding, COUNT, DUMMY> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: counts in = [";
            for (size_t i = maxEdgeCount; i < 2*maxEdgeCount; ++i) ost << node.counts[i] << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << node.counts[i] << ",";
            ost << "], self = " << node.counts[2 * maxEdgeCount];
            return ost;
          }


          /*constructor*/
          compact_multi_biedge() { counts.fill(0); };

          /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
          ~compact_multi_biedge() {}


          /**
           * @brief increments the edge counts.  does not perform flipping.
           * @param exts   input edges, should be at most 1 for in and out.
           */
          inline void update(EdgeInputType edges)
          {
            // take care of out
            uint8_t out = edges.getData()[0] & 0xF;
            if (FROM_DNA16[out] < maxEdgeCount) sat_incr(counts[FROM_DNA16[out]]);

            // take care of in.
            uint8_t in = edges.getData()[0] >> 4;
            if (FROM_DNA16[in] < maxEdgeCount) sat_incr(counts[FROM_DNA16[in] + maxEdgeCount]);

            sat_incr(counts[2 * maxEdgeCount]);
          }


          inline void merge(compact_multi_biedge const & other) {
            for (int i = 0; i <= 2 * maxEdgeCount; ++i) {
              counts[i] = sat_add(counts[i], other.counts[i]);
            }
          }

          inline COUNT get_out_edge_frequency(uint8_t idx) const {
            if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return counts[CHAR_TO_INDEX[idx]];
          }

          inline COUNT get_in_edge_frequency(uint8_t idx) const {
            if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return counts[CHAR_TO_INDEX[idx] + maxEdgeCount];
          }

          inline COUNT get_self_frequency() const {
        	  return counts[2 * maxEdgeCount];
          }

          inline uint8_t get_out_edge_count() const {
            uint8_t count = 0;
            for (size_t i = 0; i < maxEdgeCount; ++i ) {
              count += (counts[i] > 0 ? 1 : 0);
            }
            return count;
          }

          inline uint8_t get_in_edge_count() const {
            uint8_t count = 0;
            for (size_t i = maxEdgeCount; i < 2 * maxEdgeCount; ++i ) {
              count += counts[i] > 0 ? 1 : 0;
            }
            return count;
          }


          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_out_neighbors(Kmer const & kmer, std::vector<std::pair<Kmer, CountType> > & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            CountType count;
            for (unsigned char i = 0; i < maxEdgeCount; ++i) {   // no gap character
              count = counts[i];
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_in_neighbors(Kmer const & kmer, std::vector<std::pair<Kmer, CountType> > & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            CountType count;
            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // no gap character
              count = counts[i + maxEdgeCount];
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextReverseFromChar(i);
              }
            }
          }
      };


      template<typename COUNT, typename DUMMY>
      class compact_multi_biedge<::bliss::common::DNA16, COUNT, DUMMY> {

          static_assert(!::std::is_signed<COUNT>::value &&
                        ::std::is_integral<COUNT>::value, "only supports unsigned integer types for count");

#if defined(USE_MPI)
    	  friend class ::mxx::datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<::bliss::common::DNA16, COUNT, DUMMY> >;
#endif

        public:
          using EdgeEncoding = bliss::common::DNA16;
          static constexpr size_t maxEdgeCount = 5;

        protected:
          COUNT sat_add(COUNT const & a, COUNT const & b) {
            COUNT c = a + b;
            return (c < a) ? -1 : c;
          }

          void sat_incr(COUNT & target) {
            target = sat_add(target, 1);
          }


          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{     //DNA16    // bitvec
                 8,   // 0000 .     // 0x00000000
                 0,   // 0001 A     // 0x00000001
                 1,   // 0010 C     // 0x00000010
                 8,   // 0011       // 0x00000000
                 2,   // 0100 G     // 0x00000100
                 8,   // 0101       // 0x00000000
                 8,   // 0110       // 0x00000000
                 8,   // 0111       // 0x00000000
                 3,   // 1000 T     // 0x00001000
                 8,   // 1001       // 0x00000000
                 8,   // 1010       // 0x00000000
                 8,   // 1011       // 0x00000000
                 8,   // 1100       // 0x00000000
                 8,   // 1101       // 0x00000000
                 8,   // 1110       // 0x00000000
                 4    // 1111 N     // 0x00000000
          }};

          static constexpr std::array<uint8_t, 16> CHAR_TO_INDEX =
          {{     //DNA16    // bitvec
                 8,   // 0000 .     // 0x00000000
                 0,   // 0001 A     // 0x00000001
                 1,   // 0010 C     // 0x00000010
                 8,   // 0011       // 0x00000000
                 2,   // 0100 G     // 0x00000100
                 8,   // 0101       // 0x00000000
                 8,   // 0110       // 0x00000000
                 8,   // 0111       // 0x00000000
                 3,   // 1000 T     // 0x00001000
                 8,   // 1001       // 0x00000000
                 8,   // 1010       // 0x00000000
                 8,   // 1011       // 0x00000000
                 8,   // 1100       // 0x00000000
                 8,   // 1101       // 0x00000000
                 8,   // 1110       // 0x00000000
                 4    // 1111 N     // 0x00000000
          }};

          static constexpr std::array<uint8_t, 5> INDEX_TO_CHAR =
          {{     //DNA16    // bitvec
                 1,   // 0001 A     // 0x00000001
                 2,   // 0010 C     // 0x00000010
                 4,   // 0100 G     // 0x00000100
                 8,   // 1000 T     // 0x00001000
				 15   // 1111 N
          }};

        public:
          /// array of counts.  format:  [out A C G T; in A C G T], ordered for the canonical strand, not necessarily same as for the input kmer..
          std::array<COUNT, 2 * maxEdgeCount + 1> counts;  // public because of friending mxx datatype_builder does not work

          using Alphabet = EdgeEncoding;
          using CountType = COUNT;
          using EdgeInputType = bliss::debruijn::compact_simple_biedge;   // hardcoded to DNA16 because of N


          friend std::ostream& operator<<(std::ostream& ost, const compact_multi_biedge<EdgeEncoding, COUNT, DUMMY> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: counts in = [";
            for (size_t i = maxEdgeCount; i < 2*maxEdgeCount; ++i) ost << node.counts[i] << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << node.counts[i] << ",";
            ost << "], self = " << node.counts[2 * maxEdgeCount];
            return ost;
          }


          /*constructor*/
          compact_multi_biedge() { counts.fill(0); };

          /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
          ~compact_multi_biedge() {}


          /**
           * @brief increments the edge counts.  does not perform flipping.
           * @param exts   input edges, should be at most 1 for in and out.
           */
          inline void update(EdgeInputType edges)
          {
            // take care of out
            uint8_t out = edges.getData()[0] & 0xF;
            if (FROM_DNA16[out] < maxEdgeCount) sat_incr(counts[FROM_DNA16[out]]);

            // take care of in.
            uint8_t in = edges.getData()[0] >> 4;
            if (FROM_DNA16[in] < maxEdgeCount) sat_incr(counts[FROM_DNA16[in] + maxEdgeCount]);

            sat_incr(counts[2 * maxEdgeCount]);
          }


          inline void merge(compact_multi_biedge const & other) {
            for (int i = 0; i <= 2 * maxEdgeCount; ++i) {
              counts[i] = sat_add(counts[i], other.counts[i]);
            }
          }


          inline COUNT get_out_edge_frequency(uint8_t idx) const {
            if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return counts[CHAR_TO_INDEX[idx]];
          }

          inline COUNT get_in_edge_frequency(uint8_t idx) const {
            if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return counts[CHAR_TO_INDEX[idx] + maxEdgeCount];
          }


          inline COUNT get_self_frequency() const {
        	  return counts[2 * maxEdgeCount];
          }

          inline uint8_t get_out_edge_count() const {
            uint8_t count = 0;
            for (size_t i = 0; i < maxEdgeCount; ++i ) {
              count += (counts[i] > 0 ? 1 : 0);
            }
            return count;
          }

          inline uint8_t get_in_edge_count() const {
            uint8_t count = 0;
            for (size_t i = maxEdgeCount; i < 2 * maxEdgeCount; ++i ) {
              count += counts[i] > 0 ? 1 : 0;
            }
            return count;
          }


          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_out_neighbors(Kmer const & kmer, std::vector<std::pair<Kmer, CountType> > & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            CountType count;
            for (unsigned char i = 0; i < maxEdgeCount; ++i) {   // no gap character
              count = counts[i];
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_in_neighbors(Kmer const & kmer, std::vector<std::pair<Kmer, CountType> > & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            CountType count;
            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // no gap character
              count = counts[i + maxEdgeCount];
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextReverseFromChar(i);
              }
            }
          }
      };



      //=============== existence only
      //  EDGE frequency is accessed using alphabet characters by value.  all possible alphabet characters are present.

      template <typename DUMMY>
      class compact_multi_biedge<::bliss::common::DNA, bool, DUMMY> {

#if defined(USE_MPI)
    	  friend class ::mxx::datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<::bliss::common::DNA, bool, DUMMY> >;
#endif

        protected:
          ::bliss::utils::bit_ops::pop_cnt<uint8_t> pcnt;

          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{     //DNA16       // DNA
                 8,   // 0000 .     // --
                 0,   // 0001 A     // 00
                 1,   // 0010 C     // 01
                 8,   // 0011       // --
                 2,   // 0100 G     // 10
                 8,   // 0101       // --
                 8,   // 0110       // --
                 8,   // 0111       // --
                 3,   // 1000 T     // 11
                 8,   // 1001       // --
                 8,   // 1010       // --
                 8,   // 1011       // --
                 8,   // 1100       // --
                 8,   // 1101       // --
                 8,   // 1110       // --
                 8    // 1111 N     // --
          }}; // all unknowns are thrown away.

          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          uint8_t counts;

        public:
          using EdgeEncoding = bliss::common::DNA;
          using Alphabet = EdgeEncoding;
          using CountType = uint8_t;
          using EdgeInputType = bliss::debruijn::compact_simple_biedge;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;

          friend std::ostream& operator<<(std::ostream& ost, const compact_multi_biedge<EdgeEncoding, bool, DUMMY> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (size_t i = maxEdgeCount; i < 2*maxEdgeCount; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }

          /*constructor*/
          compact_multi_biedge() : counts(0) {};

          /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
          ~compact_multi_biedge() {}


          inline void update(EdgeInputType edges)
          {
            counts |= edges.getDataRef()[0];

            //std::cout << "DNA exist update" << std::endl;
//		  uint8_t in = edges.getDataRef()[0] >> 4;
//		  uint8_t out = edges.getDataRef()[0] & 0xF;
//		  printf("edge Alphabet %s val->index->bitvec: in: %u->%u->%u, out %u->%u->%u. \n",
//				 typeid(EdgeEncoding).name(),
//				 in, FROM_DNA16[in], 0x1 << (FROM_DNA16[in]),
//				 out, FROM_DNA16[out], 0x1 << (FROM_DNA16[out]));

          }

          inline void merge(compact_multi_biedge const & other) {
            counts |= other.counts;
          }

          /**
           * @brief get the number of frequency of a particular edge
           * @param idx 	ACGTN. then the rest.
           */
          inline uint8_t get_out_edge_frequency(uint8_t idx) const {
            if (idx >= maxEdgeCount) return 0;
            return (counts >> idx ) & 0x1;
          }

          inline uint8_t get_in_edge_frequency(uint8_t idx) const {
            if (idx >= maxEdgeCount) return 0;
            return (counts >> (idx+maxEdgeCount) ) & 0x1;
          }

          inline uint8_t get_self_frequency() const {
        	  return 1;
          }


          inline uint8_t get_out_edge_count() const {
            return pcnt(counts & static_cast<uint8_t>(0xF));
          }
          inline uint8_t get_in_edge_count() const {
            return pcnt(counts >> 4);
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_out_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // no gap character here.
              if (((counts >> i ) & 0x1) > 0) {
                neighbors.emplace_back(kmer);
                neighbors.back().nextFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_in_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // no gap character here.
              if (((counts >> (i+maxEdgeCount) ) & 0x1) > 0) {
                neighbors.emplace_back(kmer);
                neighbors.back().nextReverseFromChar(i);  // indices are same as alphabet chars.
              }
            }
          }


      };



      /**
       * compact multi biedge specializaiton for DNA 6.  DO NOT USE THIS FOR ASSEMBLY,
       * unless entire assembler uses the same lexicographic ordering as this edge type:
       * ACTG instead of ACGT.
       */
      template <typename DUMMY>
      class compact_multi_biedge<::bliss::common::DNA6, bool, DUMMY> {

#if defined(USE_MPI)
    	  friend class ::mxx::datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<::bliss::common::DNA6, bool, DUMMY> >;
#endif

      protected:
          ::bliss::utils::bit_ops::pop_cnt<uint8_t> pcnt;


          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{     //DNA16    // bitvec
                 8,   // 0000 .     // 0x00000000          8,
                 0,   // 0001 A     // 0x00000001          1,
                 1,   // 0010 C     // 0x00000010          3,
                 8,   // 0011       // 0x00000000          8,
                 2,   // 0100 G     // 0x00000100          6,
                 8,   // 0101       // 0x00000000          8,
                 8,   // 0110       // 0x00000000          8,
                 8,   // 0111       // 0x00000000          8,
                 3,   // 1000 T     // 0x00001000          4,
                 8,   // 1001       // 0x00000000          8,
                 8,   // 1010       // 0x00000000          8,
                 8,   // 1011       // 0x00000000          8,
                 8,   // 1100       // 0x00000000          8,
                 8,   // 1101       // 0x00000000          8,
                 8,   // 1110       // 0x00000000          8,
                 4    // 1111 N     // 0x00000000          7
          }};

          static constexpr std::array<uint8_t, 8> CHAR_TO_INDEX =
          {{     //DNA16    // bitvec
                 8,   // 000 .     // 0x00000000
                 0,   // 001 A     // 0x00000001
                 8,   // 010       // 0x00000010
                 1,   // 011 C     // 0x00000000
                 3,   // 100 T     // 0x00000100
                 8,   // 101       // 0x00000000
                 2,   // 110 G     // 0x00000000
                 4    // 111 N     // 0x00000000
          }};

          static constexpr std::array<uint8_t, 5> INDEX_TO_CHAR =
          {{     //DNA6    // bitvec
                 1,   // 001 A     // 0x00000001
                 3,   // 011 C     // 0x00000010
                 6,   // 110 G     // 0x00000100
                 4,   // 100 T     // 0x00001000
				 7,   // 111 N
          }};


          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          ::std::array<uint8_t, 2> counts;


        public:
          using EdgeEncoding = ::bliss::common::DNA6;
          using Alphabet = EdgeEncoding;
          using CountType = uint8_t;
          using EdgeInputType = bliss::debruijn::compact_simple_biedge;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;


          friend std::ostream& operator<<(std::ostream& ost, const compact_multi_biedge<EdgeEncoding, bool, DUMMY> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << (((node.counts[1] >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << (((node.counts[0] >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }


          /*constructor*/
          compact_multi_biedge() { counts.fill(0); };

          /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
          ~compact_multi_biedge() {}

          /**
           *
           * @param relative_strand
           * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
           */
          inline void update(EdgeInputType edges)
          {
              //std::cout << "DNA exist update" << std::endl;

              // take care of out
            uint8_t out = edges.getData()[0] & 0xF;
            counts[0] |= 0x1 << (FROM_DNA16[out]);

            // take care of in.
            uint8_t in = (edges.getData()[0] >> 4) & 0xF;
            counts[1] |= 0x1 << (FROM_DNA16[in]);

//		  printf("edge Alphabet %s val->index->bitvec: in: %u->%u->%u, out %u->%u->%u. \n",
//				 typeid(EdgeEncoding).name(),
//				 in, FROM_DNA16[in], 0x1 << (FROM_DNA16[in]),
//				 out, FROM_DNA16[out], 0x1 << (FROM_DNA16[out]));
          }

          inline void merge(compact_multi_biedge const & other) {
            counts[0] |= other.counts[0];
            counts[1] |= other.counts[1];
          }


          /**
           * get the number of frequency of a particular edge
           */
          inline uint8_t get_out_edge_frequency(uint8_t idx) const {
              if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return (counts[0] >> CHAR_TO_INDEX[idx]) & 0x1;
          }
          inline uint8_t get_in_edge_frequency(uint8_t idx) const {
              if (CHAR_TO_INDEX[idx] >= maxEdgeCount) return 0;
            return (counts[1] >> CHAR_TO_INDEX[idx]) & 0x1;
          }

          inline uint8_t get_self_frequency() const {
        	  return 1;
          }



          inline uint8_t get_out_edge_count() const {
            return pcnt(counts[0]);
          }
          inline uint8_t get_in_edge_count() const {
            return pcnt(counts[1]);
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_out_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            for (unsigned char i = 0; i < maxEdgeCount; ++i) { // go through all valid values.
              if (((counts[0] >> i) & 0x1) > 0) {  // gap should have frequency of 0 since its count is not incremented..
                neighbors.emplace_back(kmer);
                neighbors.back().nextFromChar(i);  // values are same as alphabet characters
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_in_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // go through all valid values.
              if (((counts[1] >> i) & 0x1) > 0) { // gap should have frequency of 0 since its count is not incremented..
                neighbors.emplace_back(kmer);
                neighbors.back().nextReverseFromChar(i);  // values are same as alphabet characters
              }
            }
          }


      };


      /*node trait class*/
      template <typename DUMMY>
      class compact_multi_biedge<::bliss::common::DNA16, bool, DUMMY> {

#if defined(USE_MPI)
    	  friend class ::mxx::datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<::bliss::common::DNA16, bool, DUMMY> >;
#endif

	  protected:
          ::bliss::utils::bit_ops::pop_cnt<uint8_t> pcnt;

          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{       //DNA16
                   16,    // 0000 .      // gap is sent to bitbucket.
                   1,     // 0001 A
                   2,     // 0010 C
                   3,     // 0011
                   4,     // 0100 G
                   5,     // 0101
                   6,     // 0110
                   7,     // 0111
                   8,     // 1000 T
                   9,     // 1001
                   10,    // 1010
                   11,    // 1011
                   12,    // 1100
                   13,    // 1101
                   14,    // 1110
                   15     // 1111 N
          }};

          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          ::std::array<uint16_t, 2> counts;

        public:
          using EdgeEncoding = ::bliss::common::DNA16;
          using Alphabet = EdgeEncoding;
          using CountType = uint16_t;
          using EdgeInputType = bliss::debruijn::compact_simple_biedge;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;

          friend std::ostream& operator<<(std::ostream& ost, const compact_multi_biedge<EdgeEncoding, bool, DUMMY> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << (((node.counts[1] >> i) & 0x1) == 1 ? 1 : 0);
            ost << "], out = [";
            for (size_t i = 0; i < maxEdgeCount; ++i) ost << (((node.counts[0] >> i) & 0x1) == 1 ? 1 : 0);
            ost << "]";
            return ost;
          }


          /*constructor*/
          compact_multi_biedge() { counts.fill(0); };

          /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
          ~compact_multi_biedge() {}

          /**
           *
           * @param relative_strand
           * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
           */
          inline void update(EdgeInputType edges)
          {
//
//             std::cout << "DNA exist update" << std::endl;
//
            // take care of out
            uint8_t out = edges.getData()[0] & 0xF;
            counts[0] |= (0x1 << (FROM_DNA16[out]));

            // take care of in.
            uint8_t in = (edges.getData()[0] >> 4) & 0xF;
            counts[1] |= (0x1 << (FROM_DNA16[in]));

//		  printf("edge Alphabet %s val->index->bitvec: in: %u->%u->%u, out %u->%u->%u. \n",
//				 typeid(EdgeEncoding).name(),
//				 in, FROM_DNA16[in], 0x1 << (FROM_DNA16[in]),
//				 out, FROM_DNA16[out], 0x1 << (FROM_DNA16[out]));

          }

          inline void merge(compact_multi_biedge const & other) {
            counts[0] |= other.counts[0];
            counts[1] |= other.counts[1];
          }

          /**
           * @brief get the number of frequency of a particular edge
           * @param idx 	ACGTN. then the rest.
           */
          inline uint8_t get_out_edge_frequency(uint8_t idx) const {
              if (idx >= maxEdgeCount) return 0;
            return (counts[0] >> idx) & 0x1;
          }
          inline uint8_t get_in_edge_frequency(uint8_t idx) const {
              if (idx >= maxEdgeCount) return 0;
            return (counts[1] >> idx) & 0x1;
          }

          inline uint8_t get_self_frequency() const {
        	  return 1;
          }


          inline uint8_t get_out_edge_count() const {
            return pcnt(counts[0]);
          }
          inline uint8_t get_in_edge_count() const {
            return pcnt(counts[1]);
          }


          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_out_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // each character in alphabet is tested.
              if (((counts[0] >> i) & 0x1) > 0) {  // gap will always have 0
                neighbors.emplace_back(kmer);
                neighbors.back().nextFromChar(i);   // then inserted.
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          template <typename Kmer>
          void get_in_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
            static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                          "kmer and edge should use the same alphabet.");
            neighbors.clear();

            for (unsigned char i = 0; i < maxEdgeCount; ++i) {   // each character in alphabet is tested.
              if (((counts[1] >> i) & 0x1) > 0) {  // gap will always have 0
                neighbors.emplace_back(kmer);
                neighbors.back().nextReverseFromChar(i);   // then inserted.
              }
            }
          }

      };

      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA, bool, DUMMY>::FROM_DNA16;
      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA6, bool, DUMMY>::FROM_DNA16;
      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA16, bool, DUMMY>::FROM_DNA16;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA, COUNT_TYPE, DUMMY>::FROM_DNA16;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA6, COUNT_TYPE, DUMMY>::FROM_DNA16;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA16, COUNT_TYPE, DUMMY>::FROM_DNA16;

      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 4> compact_multi_biedge<::bliss::common::DNA, COUNT_TYPE, DUMMY>::INDEX_TO_CHAR;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 5> compact_multi_biedge<::bliss::common::DNA6, COUNT_TYPE, DUMMY>::INDEX_TO_CHAR;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 5> compact_multi_biedge<::bliss::common::DNA16, COUNT_TYPE, DUMMY>::INDEX_TO_CHAR;

      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 4> compact_multi_biedge<::bliss::common::DNA, COUNT_TYPE, DUMMY>::CHAR_TO_INDEX;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 8> compact_multi_biedge<::bliss::common::DNA6, COUNT_TYPE, DUMMY>::CHAR_TO_INDEX;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_multi_biedge<::bliss::common::DNA16, COUNT_TYPE, DUMMY>::CHAR_TO_INDEX;





    }/*namespace graph*/
  }/*namespace debruijn*/
}/*namespace bliss*/

#if defined(USE_MPI)

namespace mxx {

template<typename A, typename T>
  struct datatype_builder<::bliss::debruijn::graph::compact_multi_biedge<A, T> > :
  public datatype_builder<decltype(::bliss::debruijn::graph::compact_multi_biedge<A, T>::counts) > {

    typedef datatype_builder<decltype(::bliss::debruijn::graph::compact_multi_biedge<A, T>::counts) > baseType;

    static MPI_Datatype get_type(){
      return baseType::get_type();
    }

    static size_t num_basic_elements() {
      return baseType::num_basic_elements();
    }
  };


template<typename A, typename T>
  struct datatype_builder<const ::bliss::debruijn::graph::compact_multi_biedge<A, T> > :
  public datatype_builder<decltype(::bliss::debruijn::graph::compact_multi_biedge<A, T>::counts) > {

    typedef datatype_builder<decltype(::bliss::debruijn::graph::compact_multi_biedge<A, T>::counts) > baseType;

    static MPI_Datatype get_type(){
      return baseType::get_type();
    }

    static size_t num_basic_elements() {
      return baseType::num_basic_elements();
    }
  };

} // mxx

#endif


#endif /* DEBRUIJN_GRAPH_NODE_HPP_ */
