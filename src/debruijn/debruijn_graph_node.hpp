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
 *      Author: tony pan
 *
 */

#ifndef DEBRUIJN_GRAPH_NODE_HPP_
#define DEBRUIJN_GRAPH_NODE_HPP_

#include "bliss-config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

//#if defined(USE_OPENMP)
//#include "omp.h"
//#endif

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.
#include <tuple>        // tuple and utility functions
#include <utility>      // pair and utility functions.
#include <type_traits>

#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "utils/bit_ops.hpp"

#include "debruijn/edge_iterator.hpp"

namespace bliss
{
	namespace debruijn
	{
		namespace graph
		{



      template <typename EdgeEncoding, typename COUNT = bool, typename DUMMY = void>
      class compact_edge;

			/**
			 * @brief kmer metadata holding the incoming and outgoing edge counts in a de bruijn graph.
			 * @details lower 4 integers are out edge counts, upper 4 integers are in edge counts.
			 *        the out and in edges are relative to the kmer's orientation.
			 * @tparam EdgeEncoding  alphabet used for edge encoding.
			 */
			template<typename COUNT, typename DUMMY>
			class compact_edge<::bliss::common::DNA, COUNT, DUMMY> {

			    static_assert(!::std::is_signed<COUNT>::value &&
			                  ::std::is_integral<COUNT>::value, "only supports unsigned integer types for count");


			  public:
          using EdgeEncoding = bliss::common::DNA;
          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;

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


          /// array of counts.  format:  [out A C G T; in A C G T], ordered for the canonical strand, not necessarily same as for the input kmer..
          std::array<COUNT, 2 * maxEdgeCount> counts;



			  public:

          using Alphabet = EdgeEncoding;
          using CountType = COUNT;
          using EdgeInputType = bliss::debruijn::SimpleBiEdgeType;   // hardcoded to DNA16 because of N


	        friend std::ostream& operator<<(std::ostream& ost, const compact_edge<EdgeEncoding, COUNT, DUMMY> & node)
	        {
	          // friend keyword signals that this overrides an externally declared function
	          ost << " dBGr node: counts in = [";
	          for (size_t i = maxEdgeCount; i < 2*maxEdgeCount; ++i) ost << node.counts[i] << ",";
	          ost << "], out = [";
	          for (size_t i = 0; i < maxEdgeCount; ++i) ost << node.counts[i] << ",";
	          ost << "]";
	          return ost;
	        }


				/*constructor*/
				compact_edge() { counts.fill(0); };

				/*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
				~compact_edge() {}


				/**
				 * @brief increments the edge counts.  does not perform flipping.
				 * @param exts   input edges, should be at most 1 for in and out.
				 */
        void update(EdgeInputType edges)
        {
          // take care of out
          uint8_t out = edges.getData()[0] & 0xF;
          sat_incr(counts[FROM_DNA16[out]]);

          // take care of in.
          uint8_t in = edges.getData()[0] >> 4;
          sat_incr(counts[FROM_DNA16[in] + maxEdgeCount]);
        }



        void merge(compact_edge const & other) {
          for (int i = 0; i < 2 * maxEdgeCount; ++i) {
            counts[i] = sat_add(counts[i], other.counts[i]);
          }
        }



        COUNT get_out_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;
          return counts[idx];
        }


				COUNT get_in_edge_frequency(uint8_t idx) const {
				  if (idx >= maxEdgeCount) return 0;
				  return counts[idx + maxEdgeCount];
				}

        uint8_t get_out_edge_count() const {
          uint8_t count = 0;
          for (size_t i = 0; i < maxEdgeCount; ++i ) {
        	  count += (counts[i] > 0 ? 1 : 0);
          }
          return count;
        }

        uint8_t get_in_edge_count() const {
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
          for (size_t i = 0; i < maxEdgeCount; ++i) {   // no gap character
            count = get_out_edge_frequency(i);
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
          for (size_t i = 0; i < maxEdgeCount; ++i) {  // no gap character
            count = get_in_edge_frequency(i);
            if (count > 0) {
              neighbors.emplace_back(kmer, count);
              neighbors.back().first.nextReverseFromChar(i);
            }
          }
        }

			};


      template <typename DUMMY>
      class compact_edge<::bliss::common::DNA, bool, DUMMY> {
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
          using EdgeInputType = bliss::debruijn::SimpleBiEdgeType;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;

          friend std::ostream& operator<<(std::ostream& ost, const compact_edge<EdgeEncoding, bool, DUMMY> & node)
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
        compact_edge() : counts(0) {};

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~compact_edge() {}


        void update(EdgeInputType edges)
        {
          counts |= edges.getDataRef()[0];

//          uint8_t in = edges.getDataRef()[0] >> 4;
//          uint8_t out = edges.getDataRef()[0] & 0xF;
//
//          printf("edge Alphabet %s val->index->bitvec: in: %u->%u->%u, out %u->%u->%u. \n",
//                 typeid(EdgeEncoding).name(),
//                 in, FROM_DNA16[in], 0x1 << (FROM_DNA16[in]),
//                 out, FROM_DNA16[out], 0x1 << (FROM_DNA16[out]));

        }

        void merge(compact_edge const & other) {
          counts |= other.counts;
        }


        /**
         * get the number of frequency of a particular edge
         */
        uint8_t get_out_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;
          return (counts >> idx ) & 0x1;
        }

        uint8_t get_in_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;
          return (counts >> (idx+maxEdgeCount) ) & 0x1;
        }

        uint8_t get_out_edge_count() const {
          return pcnt(counts & static_cast<uint8_t>(0xF));
        }
        uint8_t get_in_edge_count() const {
          return pcnt(counts >> 4);
        }

        // construct a new kmer from a known edge, if that edge's count is non-zero
        template <typename Kmer>
        void get_out_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
          static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                        "kmer and edge should use the same alphabet.");
          neighbors.clear();

          for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // no gap character here.
            if (get_out_edge_frequency(i) > 0) {
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

          for (size_t i = 0; i < maxEdgeCount; ++i) {  // no gap character here.
            if (get_in_edge_frequency(i) > 0) {
              neighbors.emplace_back(kmer);
              neighbors.back().nextReverseFromChar(i);
            }
          }
        }


      };



      /*node trait class*/
      template <typename DUMMY>
      class compact_edge<::bliss::common::DNA6, bool, DUMMY> {
        protected:
          ::bliss::utils::bit_ops::pop_cnt<uint8_t> pcnt;


          // convert DNA16 value to bit array index.  not strictly needed here but have it for consistency
          // any unmapped gets shifted out to bit 8.
          static constexpr std::array<uint8_t, 16> FROM_DNA16 =
          {{       //DNA16      // DNA5
            8,     // 0000 .    // 000   // has a mapping, but throw away gap and invalid characters.
            1,     // 0001 A    // 001
            3,     // 0010 C    // 011
            8,     // 0011      // ---
            6,     // 0100 G    // 110
            8,     // 0101      // ---
            8,     // 0110      // ---
            8,     // 0111      // ---
            4,     // 1000 T    // 100
            8,     // 1001      // ---
            8,     // 1010      // ---
            8,     // 1011      // ---
            8,     // 1100      // ---
            8,     // 1101      // ---
            8,     // 1110      // ---
            7      // 1111 N    // 111
          }};

          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          ::std::array<uint8_t, 2> counts;



        public:
          using EdgeEncoding = ::bliss::common::DNA6;
          using Alphabet = EdgeEncoding;
          using CountType = uint8_t;
          using EdgeInputType = bliss::debruijn::SimpleBiEdgeType;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;


          friend std::ostream& operator<<(std::ostream& ost, const compact_edge<EdgeEncoding, bool, DUMMY> & node)
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
        compact_edge() { counts.fill(0); };

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~compact_edge() {}

        /**
         *
         * @param relative_strand
         * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
         */
        void update(EdgeInputType edges)
        {
          // take care of out
          uint8_t out = edges.getData()[0] & 0xF;
          counts[0] |= 0x1 << (FROM_DNA16[out]);

          // take care of in.
          uint8_t in = (edges.getData()[0] >> 4) & 0xF;
          counts[1] |= 0x1 << (FROM_DNA16[in]);

//          printf("edge Alphabet %s val->index->bitvec: in: %u->%u->%u, out %u->%u->%u. \n",
//                 typeid(EdgeEncoding).name(),
//                 in, FROM_DNA16[in], 0x1 << (FROM_DNA16[in]),
//                 out, FROM_DNA16[out], 0x1 << (FROM_DNA16[out]));
        }

        void merge(compact_edge const & other) {
          counts[0] |= other.counts[0];
          counts[1] |= other.counts[1];
        }


        /**
         * get the number of frequency of a particular edge
         */
        uint8_t get_out_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;
          return (counts[0] >> idx) & 0x1;
        }
        uint8_t get_in_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;
          return (counts[1] >> idx) & 0x1;
        }

        uint8_t get_out_edge_count() const {
          return pcnt(counts[0]);
        }
        uint8_t get_in_edge_count() const {
          return pcnt(counts[1]);
        }

        // construct a new kmer from a known edge, if that edge's count is non-zero
        template <typename Kmer>
        void get_out_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
          static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                        "kmer and edge should use the same alphabet.");
          neighbors.clear();

          for (unsigned char i = 0; i < maxEdgeCount; ++i) { // go through all valid values.
            if (get_out_edge_frequency(i) > 0) {  // gap should have frequency of 0 since its count is not incremented..
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

          for (size_t i = 0; i < maxEdgeCount; ++i) {  // go through all valid values.
            if (get_in_edge_frequency(i) > 0) { // gap should have frequency of 0 since its count is not incremented..
              neighbors.emplace_back(kmer);
              neighbors.back().nextReverseFromChar(i);
            }
          }
        }


      };


      /*node trait class*/
      template <typename DUMMY>
      class compact_edge<::bliss::common::DNA16, bool, DUMMY> {
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
          using EdgeInputType = bliss::debruijn::SimpleBiEdgeType;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE;

          friend std::ostream& operator<<(std::ostream& ost, const compact_edge<EdgeEncoding, bool, DUMMY> & node)
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
        compact_edge() { counts.fill(0); };

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~compact_edge() {}

        /**
         *
         * @param relative_strand
         * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
         */
        void update(EdgeInputType edges)
        {
          // take care of out
          uint8_t out = edges.getData()[0] & 0xF;
          counts[0] |= (0x1 << (FROM_DNA16[out]));

          // take care of in.
          uint8_t in = (edges.getData()[0] >> 4) & 0xF;
          counts[1] |= (0x1 << (FROM_DNA16[in]));

//          printf("edge Alphabet %s val->index->bitvec: in: %u->%u->%u, out %u->%u->%u. \n",
//                 typeid(EdgeEncoding).name(),
//                 in, FROM_DNA16[in], 0x1 << (FROM_DNA16[in]),
//                 out, FROM_DNA16[out], 0x1 << (FROM_DNA16[out]));

        }

        void merge(compact_edge const & other) {
          counts[0] |= other.counts[0];
          counts[1] |= other.counts[1];
        }


        /**
         * get the number of frequency of a particular edge
         */
        uint8_t get_out_edge_frequency(uint8_t idx) const {
          return (counts[0] >> idx) & 0x1;
        }
        uint8_t get_in_edge_frequency(uint8_t idx) const {
          return (counts[1] >> idx) & 0x1;
        }

        uint8_t get_out_edge_count() const {
          return pcnt(counts[0]);
        }
        uint8_t get_in_edge_count() const {
          return pcnt(counts[1]);
        }


        // construct a new kmer from a known edge, if that edge's count is non-zero
        template <typename Kmer>
        void get_out_neighbors(Kmer const & kmer, std::vector<Kmer> & neighbors) {
          static_assert(std::is_same<typename Kmer::KmerAlphabet, Alphabet>::value,
                        "kmer and edge should use the same alphabet.");
          neighbors.clear();

          for (unsigned char i = 0; i < maxEdgeCount; ++i) {  // gap is not counted
            if (get_out_edge_frequency(i) > 0) {  // gap will always have 0
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

          for (size_t i = 0; i < maxEdgeCount; ++i) {  // gap is not counted.
            if (get_in_edge_frequency(i) > 0) {  // gap will always have 0
              neighbors.emplace_back(kmer);
              neighbors.back().nextReverseFromChar(i);
            }
          }
        }

      };

      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_edge<::bliss::common::DNA, bool, DUMMY>::FROM_DNA16;
      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_edge<::bliss::common::DNA6, bool, DUMMY>::FROM_DNA16;
      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_edge<::bliss::common::DNA16, bool, DUMMY>::FROM_DNA16;
      template <typename COUNT_TYPE, typename DUMMY>
      constexpr std::array<uint8_t, 16> compact_edge<::bliss::common::DNA, COUNT_TYPE, DUMMY>::FROM_DNA16;



		}/*namespace graph*/
	}/*namespace debruijn*/
}/*namespace bliss*/




#endif /* DEBRUIJN_GRAPH_NODE_HPP_ */
