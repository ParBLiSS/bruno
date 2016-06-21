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
 * de_bruijn_node_trait.hpp
 *
 *  Created on: Jul 31, 2015
 *      Author: yongchao
 *      Author: tony pan
 *
 */

#ifndef DE_BRUIJN_NODE_TRAIT_HPP_
#define DE_BRUIJN_NODE_TRAIT_HPP_

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

namespace bliss
{
	namespace de_bruijn
	{
		namespace node
		{
      static constexpr unsigned char SENSE = 0;
      static constexpr unsigned char ANTI_SENSE = 1;

      // utilities operating on nodes
      template <typename Kmer, typename EdgeType>
      class node_utils {
        public:

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_out_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<Kmer> & neighbors) {
            neighbors.clear();

            for (size_t i = 0; i < EdgeType::maxHalfEdgeCount; ++i) {
              if (edge.get_edge_frequency(i) > 0) {
                neighbors.emplace_back(kmer);
                neighbors.back().nextFromChar(Kmer::KmerAlphabet::FROM_INDEX[i]);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_in_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<Kmer> & neighbors) {
            neighbors.clear();

            for (size_t i = 0; i < EdgeType::maxHalfEdgeCount; ++i) {
              if (edge.get_edge_frequency(i + EdgeType::maxHalfEdgeCount) > 0) {
                neighbors.emplace_back(kmer);
                neighbors.back().nextReverseFromChar(Kmer::KmerAlphabet::FROM_INDEX[i]);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_out_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<std::pair<Kmer, typename EdgeType::CountType> > & neighbors) {
            neighbors.clear();

            typename EdgeType::CountType count;
            for (size_t i = 0; i < EdgeType::maxHalfEdgeCount; ++i) {
              count = edge.get_edge_frequency(i);
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextFromChar(Kmer::KmerAlphabet::FROM_INDEX[i]);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_in_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<std::pair<Kmer, typename EdgeType::CountType> > & neighbors) {
            neighbors.clear();

            typename EdgeType::CountType count;

            for (size_t i = 0; i < EdgeType::maxHalfEdgeCount; ++i) {
              count = edge.get_edge_frequency(i + EdgeType::maxHalfEdgeCount);
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextReverseFromChar(Kmer::KmerAlphabet::FROM_INDEX[i]);
              }
            }
          }


      };

			/**
			 * @brief kmer metadata holding the incoming and outgoing edge counts in a de bruijn graph.
			 * @details lower 4 integers are out edge counts, upper 4 integers are in edge counts.
			 *        the out and in edges are relative to the kmer's orientation.
			 * @tparam EdgeEncoding  alphabet used for edge encoding.
			 */
			template<typename EdgeEncoding, typename COUNT = uint32_t>
			class edge_counts {

			    static_assert(!::std::is_signed<COUNT>::value &&
			                  ::std::is_integral<COUNT>::value, "only supports unsigned integer types for count");

			  protected:
			    COUNT sat_add(COUNT const & a, COUNT const & b) {
			      COUNT c = a + b;
			      return (c < a) ? -1 : c;
			    }

			    void sat_incr(COUNT & target) {
			      target = sat_add(target, 1);
			    }


			  public:

          using Alphabet = EdgeEncoding;
          using CountType = COUNT;
          using EdgeInputType = bliss::common::Kmer<2, ::bliss::common::DNA16, uint8_t>;   // hardcoded to DNA16 because of N

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE * 2;
          static constexpr size_t maxHalfEdgeCount = EdgeEncoding::SIZE;

	        friend std::ostream& operator<<(std::ostream& ost, const edge_counts<EdgeEncoding, COUNT> & node)
	        {
	          // friend keyword signals that this overrides an externally declared function
	          ost << " dBGr node: counts in = [";
	          for (size_t i = maxHalfEdgeCount; i < maxEdgeCount; ++i) ost << node.counts[i] << ",";
	          ost << "], out = [";
	          for (size_t i = 0; i < maxHalfEdgeCount; ++i) ost << node.counts[i] << ",";
	          ost << "]";
	          return ost;
	        }


			    /// array of counts.  format:  [out A C G T; in A C G T], ordered for the canonical strand, not necessarily same as for the input kmer..
			    std::array<COUNT, maxEdgeCount> counts;

				/*constructor*/
				edge_counts() { counts.fill(0); };

				/*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
				~edge_counts() {}


				/**
				 * @brief increments the edge counts.  does not perform flipping.
				 * @param exts   input edges, should be at most 1 for in and out.
				 */
        void update(EdgeInputType edges)
        {
          // take care of out
          uint8_t out = edges.getData()[0] & 0xF;
          sat_incr(counts[EdgeEncoding::TO_INDEX[out]]);

          // take care of in.
          uint8_t in = (edges.getData()[0] >> 4) & 0xF;
          sat_incr(counts[EdgeEncoding::TO_INDEX[in] + maxHalfEdgeCount]);
        }



        void merge(edge_counts const & other) {
          for (int i = 0; i < maxEdgeCount; ++i) {
            counts[i] = sat_add(counts[i], other.counts[i]);
          }
        }





				COUNT get_edge_frequency(uint8_t idx) const {
				  if (idx >= maxEdgeCount) return 0;

				  return counts[idx];
				}

        uint8_t get_out_edge_count() const {
          uint8_t count = 0;
          for (size_t i = 0; i < maxHalfEdgeCount; ++i ) {
        	  count += (counts[i] > 0 ? 1 : 0);
          }
          return count;
        }

        uint8_t get_in_edge_count() const {
            uint8_t count = 0;
            for (size_t i = maxHalfEdgeCount; i < maxEdgeCount; ++i ) {
          	  count += counts[i] > 0 ? 1 : 0;
            }
          return count;
        }

			};


      /*node trait class*/
      template<typename EdgeEncoding>
      class edge_exists;

      template <>
      class edge_exists<::bliss::common::DNA> {
        public:
          using EdgeEncoding = bliss::common::DNA;
          using Alphabet = EdgeEncoding;
          using CountType = uint8_t;
          using EdgeInputType = bliss::common::Kmer<2, ::bliss::common::DNA16, uint8_t>;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE * 2;
          static constexpr size_t maxHalfEdgeCount = EdgeEncoding::SIZE;

          friend std::ostream& operator<<(std::ostream& ost, const edge_exists<EdgeEncoding> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (size_t i = maxHalfEdgeCount; i < maxEdgeCount; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxHalfEdgeCount; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }


          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          uint8_t counts;

        /*constructor*/
        edge_exists() : counts(0) {};

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~edge_exists() {}


        void update(EdgeInputType edges)
        {
          counts |= edges.getDataRef()[0];
        }

        void merge(edge_exists const & other) {
          counts |= other.counts;
        }


        /**
         * get the number of frequency of a particular edge
         */
        uint8_t get_edge_frequency(uint8_t idx) const {
          if (idx >= 8) return 0;

          return (counts >> idx) & 0x1;
        }

        uint8_t get_out_edge_count() const {
          uint8_t count = (counts & 0x1);
          count += ((counts >> 1) & 0x1);
          count += ((counts >> 2) & 0x1);
          count += ((counts >> 3) & 0x1);
          return count;
        }
        uint8_t get_in_edge_count() const {
          uint8_t count = ((counts >> 4) & 0x1);
          count += ((counts >> 5) & 0x1);
          count += ((counts >> 6) & 0x1);
          count += ((counts >> 7) & 0x1);
          return count;
        }


      };



      /*node trait class*/
      template<>
      class edge_exists<::bliss::common::DNA6> {
        public:
          using EdgeEncoding = ::bliss::common::DNA6;
          using Alphabet = EdgeEncoding;
          using CountType = uint8_t;
          using EdgeInputType = bliss::common::Kmer<2, ::bliss::common::DNA16, uint8_t>;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE * 2;
          static constexpr size_t maxHalfEdgeCount = EdgeEncoding::SIZE;


          friend std::ostream& operator<<(std::ostream& ost, const edge_exists<EdgeEncoding> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (size_t i = 0; i < maxHalfEdgeCount; ++i) ost << (((node.counts[1] >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxHalfEdgeCount; ++i) ost << (((node.counts[0] >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }


          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          ::std::array<uint8_t, 2> counts;

        /*constructor*/
        edge_exists() { counts.fill(0); };

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~edge_exists() {}

        /**
         *
         * @param relative_strand
         * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
         */
        void update(EdgeInputType edges)
        {
          // take care of out
          uint8_t out = edges.getData()[0] & 0xF;
          counts[0] |= 0x1 << (EdgeEncoding::TO_INDEX[out]);

          // take care of in.
          uint8_t in = (edges.getData()[0] >> 4) & 0xF;
          counts[1] |= 0x1 << (EdgeEncoding::TO_INDEX[in]);
        }

        void merge(edge_exists const & other) {
          counts[0] |= other.counts[0];
          counts[1] |= other.counts[1];
        }


        /**
         * get the number of frequency of a particular edge
         */
        uint8_t get_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;

          return (counts[(idx / maxHalfEdgeCount)] >> (idx % maxHalfEdgeCount)) & 0x1;
        }

        uint8_t get_out_edge_count() const {
          uint8_t count = 0;
          for (size_t i = 0; i < maxHalfEdgeCount; ++i) {
            count += ((counts[0] >> i) & 0x1);
          }
          return count;
        }
        uint8_t get_in_edge_count() const {
          uint8_t count = 0;
          for (size_t i = 0; i < maxHalfEdgeCount; ++i) {
            count += ((counts[1] >> i) & 0x1);
          }
          return count;
        }

      };



      /*node trait class*/
      template<>
      class edge_exists<::bliss::common::DNA16> {
        public:
          using EdgeEncoding = ::bliss::common::DNA16;
          using Alphabet = EdgeEncoding;
          using CountType = uint8_t;
          using EdgeInputType = bliss::common::Kmer<2, ::bliss::common::DNA16, uint8_t>;

          static constexpr size_t maxEdgeCount = EdgeEncoding::SIZE * 2;
          static constexpr size_t maxHalfEdgeCount = EdgeEncoding::SIZE;


          friend std::ostream& operator<<(std::ostream& ost, const edge_exists<EdgeEncoding> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (size_t i = 0; i < maxHalfEdgeCount; ++i) ost << (((node.counts[1] >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (size_t i = 0; i < maxHalfEdgeCount; ++i) ost << (((node.counts[0] >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }


          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          ::std::array<uint16_t, 2> counts;

        /*constructor*/
        edge_exists() { counts.fill(0); };

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~edge_exists() {}

        /**
         *
         * @param relative_strand
         * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
         */
        void update(EdgeInputType edges)
        {
          // take care of out
          uint8_t out = edges.getData()[0] & 0xF;
          counts[0] |= 0x1 << (EdgeEncoding::TO_INDEX[out]);

          // take care of in.
          uint8_t in = (edges.getData()[0] >> 4) & 0xF;
          counts[1] |= 0x1 << (EdgeEncoding::TO_INDEX[in]);
        }

        void merge(edge_exists const & other) {
          counts[0] |= other.counts[0];
          counts[1] |= other.counts[1];
        }


        /**
         * get the number of frequency of a particular edge
         */
        uint8_t get_edge_frequency(uint8_t idx) const {
          if (idx >= maxEdgeCount) return 0;

          return (counts[(idx / maxHalfEdgeCount)] >> (idx % maxHalfEdgeCount)) & 0x1;
        }

        uint8_t get_out_edge_count() const {
          uint8_t count = 0;
          for (size_t i = 0; i < maxHalfEdgeCount; ++i) {
            count += ((counts[0] >> i) & 0x1);
          }
          return count;
        }
        uint8_t get_in_edge_count() const {
          uint8_t count = 0;
          for (size_t i = 0; i < maxHalfEdgeCount; ++i) {
            count += ((counts[1] >> i) & 0x1);
          }
          return count;
        }

      };



		}/*namespace node*/
	}/*namespace de_bruijn*/
}/*namespace bliss*/




#endif /* DE_BRUIJN_NODE_TRAIT_HPP_ */
