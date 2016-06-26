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
 * debruijn_chain_operations.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_CHAIN_OPERATIONS_HPP_
#define DEBRUIJN_CHAIN_OPERATIONS_HPP_

#include <tuple>
#include <utility>
#include "debruijn/debruijn_chain_node.hpp"

namespace bliss {
  namespace debruijn {

  template <typename KMER>
      struct lex_less;

    namespace operation {

      static constexpr char IN = -1;
      static constexpr char OUT = 1;

      namespace chain {


        // k-mer is the target k-mer ON SAME STRAND as the canonical key kmer that this is used with.
        //    i.e. if key kmer is canonical, Kmer in field is on same strand (whether Kmer itself is canonical or not)
        // char is an indicator of edge orientation.  + means OUT edge, - means IN edge.
        template <typename Kmer>
        using terminus_update_md = ::std::pair<Kmer, char>;


        /// update the terminal de bruijn chain nodes to indicate it as such.  used on neighbors of branch points.
        /// return 1 if updated, 0 if not updated
        template <typename Kmer>
        struct terminus_update {
            size_t operator()(::bliss::debruijn::simple_biedge<Kmer> & x,
                              terminus_update_md<Kmer> const & y)  {
              assert(y.second != 0);  // second entry should not be 0.

              if (y.second < 0) {  // IN edge

                if ((std::get<2>(x) == 0) || (y.first != std::get<0>(x))) {
                  std::cout << "y: " << bliss::utils::KmerUtils::toASCIIString(y.first) <<
                      ", y ?= x[in] " << std::get<2>(x) << ": " << bliss::utils::KmerUtils::toASCIIString(std::get<0>(x)) <<
                      ", x[out] " << std::get<3>(x) << ": " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(x)) << std::endl;
                }

                assert(std::get<2>(x) == 1);   // just to make sure that the edge we're trying to update is not already "cut".
                assert(y.first == std::get<0>(x));   // just to make sure that we are talking about the same kmer.

                std::get<2>(x) = 0;   // mark as terminus
                return 1;
              } else if (y.second > 0) {  // OUT edge
                if ((std::get<3>(x) == 0) || (y.first != std::get<1>(x))) {
                  std::cout << "y: " << bliss::utils::KmerUtils::toASCIIString(y.first) <<
                      ", x[in] " << std::get<2>(x) << ": " <<  bliss::utils::KmerUtils::toASCIIString(std::get<0>(x)) <<
                      ", y ?= x[out] " << std::get<3>(x) << ": " << bliss::utils::KmerUtils::toASCIIString(std::get<1>(x)) <<  std::endl;
                }

                assert(std::get<3>(x) == 1);   // just to make sure that the edge we're trying to update is not already "cut".
                assert(y.first == std::get<1>(x));   // just to make sure that we are talking about the same kmer.

                std::get<3>(x) = 0;   // mark as terminus
                return 1;
              }
              return 0;
            }
        };


        // k-mer is the target k-mer ON SAME STRAND as the canonical key kmer that this is used with.
        //    i.e. if key kmer is canonical, Kmer in field is on same strand (whether Kmer itself is canonical or not)
        // int is the distance between Kmer and key kmer on the key kmer's strand.  - means Kmer is a terminal node, + means it's not..
        // char is indicator for edge orientation.   + means OUT edge, - means IN edge.
        template <typename Kmer>
        using chain_update_md = ::std::tuple<Kmer, int, char>;


        /// update the internal chains.  used on chains only.
        // the update metadata is on the same strand as the kmer key (canonical).  the orientation is adjusted accordingly by lex_less
        template <typename Kmer>
        struct chain_update {

            /// update the chain node.  return 1 if updated, 0 if not.
            size_t operator()(::bliss::debruijn::simple_biedge<Kmer> & x,
                              chain_update_md<Kmer> const & y)  {

              assert(std::get<2>(y) != 0);   // orientation needs to be 1 or -1

              // ***
              //  received an update for the current kmer.  check to see if we are jumping sufficiently

              auto dist = abs(std::get<1>(y));

              assert(dist > 0);   // update distance should be larger than 0.


              if (std::get<2>(y) < 0) {  // IN edge
                if ( ( std::get<2>(x) == 0) ||
                    ! ( ( ( std::get<2>(x) < 0 ) && ((dist < -(std::get<2>(x)) ) || ((dist == -(std::get<2>(x))) && (std::get<0>(y) == std::get<0>(x)))) ) ||
                        ( ( std::get<2>(x) > 0 ) && ((dist != std::get<2>(x) ) || ((dist == std::get<2>(x)) && (std::get<0>(y) == std::get<0>(x)))) )) )
                {

                  std::cout << "NEW:\tin dist " << std::get<1>(y) << " kmer: " << std::get<0>(y) << std::endl;
                  std::cout << "\tin dist " << std::get<2>(x) << " kmer: " << std::get<0>(x) << std::endl;
                  std::cout << "\tout dist " << std::get<3>(x) << " kmer: " << std::get<1>(x) << std::endl;
                }


                // some checks
                // current node not be a terminus. and if it's negative (finished), the sent update should point to same.
                // k-mer sending update is between current and end, and should have finished sooner.
                // terminus (x[3] == 0) should not get an update itself, since nothing is to the right of it to send it update.
                assert( std::get<2>(x) != 0);
                assert( ( ( std::get<2>(x) < 0 ) && ((dist < -(std::get<2>(x)) ) || ((dist == -(std::get<2>(x))) && (std::get<0>(y) == std::get<0>(x)))) ) ||
                        ( ( std::get<2>(x) > 0 ) && ((dist != std::get<2>(x) ) || ((dist == std::get<2>(x)) && (std::get<0>(y) == std::get<0>(x)))) ));
                // if current node not pointing to terminus, update distance should be larger than current distance,
                // or if equal, should have same k-mer.
                // cannot be smaller.

                // update out edge IF new distance is larger than current.
                //                if ((std::get<2>(x) > 0) && (dist > std::get<2>(x))) {
                //                  std::get<2>(x) = std::get<1>(y);   // note that if update indicates finished, it's propagated.
                //                  std::get<0>(x) = std::get<0>(y);
                //                }
                if (std::get<2>(x) > 0) {
                  if (dist > std::get<2>(x)) std::get<0>(x) = std::get<0>(y);
                  if (dist >= std::get<2>(x)) {
                    std::get<2>(x) = std::get<1>(y);   // note that if update indicates finished, it's propagated.
                    return 1;
                  }
                }

              } else if (std::get<2>(y) > 0) {  // OUT edge

                if ( (std::get<3>(x) == 0) ||
                    !(( ( std::get<3>(x) < 0 ) && ((dist < -(std::get<3>(x)) ) || ((dist == -(std::get<3>(x))) && (std::get<0>(y) == std::get<1>(x)))) ) ||
                        ( ( std::get<3>(x) > 0 ) && ((dist != std::get<3>(x) ) || ((dist == std::get<3>(x)) && (std::get<0>(y) == std::get<1>(x)))) ))
                ) {
                  std::cout << "NEW:\tout dist " << std::get<1>(y) << " kmer: " << std::get<0>(y) << std::endl;
                  std::cout << "\tin dist " << std::get<2>(x) << " kmer: " << std::get<0>(x) << std::endl;
                  std::cout << "\tout dist " << std::get<3>(x) << " kmer: " << std::get<1>(x) << std::endl;
                }

                // some checks
                // current node not be a terminus. and if it's negative (finished), the sent update should point to same.
                // k-mer sending update is between current and end, and should have finished sooner.
                // terminus (x[3] == 0) should not get an update itself, since nothing is to the right of it to send it update.
                assert( std::get<3>(x) != 0);
                assert( ( ( std::get<3>(x) < 0 ) && ((dist < -(std::get<3>(x)) ) || ((dist == -(std::get<3>(x))) && (std::get<0>(y) == std::get<1>(x)))) ) ||
                        ( ( std::get<3>(x) > 0 ) && ((dist != std::get<3>(x) ) || ((dist == std::get<3>(x)) && (std::get<0>(y) == std::get<1>(x)))) ));
                // if current node not pointing to terminus, update distance should be larger than current distance,
                // or if equal, should have same k-mer.
                // cannot be smaller.
                //            	if (std::get<1>(y) < 0) {
                //            		printf("out terminal old out = %d, new out = %d\n", std::get<3>(x), std::get<1>(y));
                //            		std::cout << "old out kmer " << std::get<1>(x) << " new out kmer " << std::get<0>(y) << std::endl;
                //            	}

                // update out edge IF new distance is larger than current.
                if (std::get<3>(x) > 0) {
                  if (dist > std::get<3>(x)) std::get<1>(x) = std::get<0>(y);
                  if (dist >= std::get<3>(x)) {
                    std::get<3>(x) = std::get<1>(y);   // note that if update indicates finished, it's propagated.
                    return 1;
                  }
                }

              }
              return 0;
            }
        };



        template <typename KMER>
        struct chain_node_to_char_transform {
            // convert a chain node to a compacted node (for print out).
            inline ::bliss::debruijn::chain::compacted_chain_node<KMER> operator()(::std::pair<KMER, ::bliss::debruijn::simple_biedge<KMER> > const & x) const {

            	// APPROACH:
            	// 1. choose an endpoint.  compare the 2 end points and choose the smaller.
            	//    note that we want to eventually print, so we want to maintain a 5' to 3' order, regardless of whether the
            	//    chosen endpoint kmer is canonical or not.
            	//    also because of this, look at only Left to current to right, or
            	//        revcomp(right) to revcomp(current) to revcomp(left), as current is canonical.
            	//    we therefore compare Left (incoming) kmer to reverse complement of Right (out going) kmer.
            	//    this in effect chooses the strand to print, each corresponding to 1 of 2 end points of chain.


            	// 2. if we choose Left, then we get the last character in the current kmer
            	//    else we get the last character in the reverse complement of the current kmer  (same as complement of first char)
            	// 3. we also choose the corresponding distance.

              // if this is end node, L is set to node k-mer.  else left edge k-mer (and right k-mer revcomp)
              KMER L = (std::get<2>(x.second) == 0) ? x.first : std::get<0>(x.second);
              KMER R = ((std::get<3>(x.second) == 0) ? x.first : std::get<1>(x.second)).reverse_complement();

              // now compare L and R, the smaller is the "representative" of a chain
              // choosing one also chooses a strand
              bool sense = (L < R);

              // once the representative is chosen, the distance is set too.
              int dist = sense ? abs(std::get<2>(x.second)) : abs(std::get<3>(x.second));

              // recall that left and right are on same strand as kmer.
              // so if we have sense (L/IN), we want the last char from the kmer as the next
              //    if we have antisense (R/OUT), we want the complement of first char from the kmer as the next
              char ch = sense ? (x.first.getData()[0] & ((0x1 << KMER::bitsPerChar) - 1)) :
            		  KMER::KmerAlphabet::TO_COMPLEMENT[(x.first >> (KMER::size - 1)).getData()[0] & ((0x1 << KMER::bitsPerChar) - 1)];
//
//              auto md = x.second;
//              if ((std::get<2>(md) > -2) || (std::get<3>(md) > -2)) {
//				  std::cout << "\tin dist " << std::get<2>(md) << " kmer: " << std::get<0>(md) << std::endl;
//				  std::cout << "\tkmer: " << x.first << std::endl;
//				  std::cout << "\tout dist " << std::get<3>(md) << " kmer: " << std::get<1>(md) << std::endl;
//				  std::cout << "\tchar = " << KMER::KmerAlphabet::TO_ASCII[ch] << std::endl;
//              }


              // store the canonical end, and dist.
              return ::bliss::debruijn::chain::compacted_chain_node<KMER>(sense ? L : R, dist, ch);
            }
        };

        template <typename KMER>
        struct chain_node_less {
            inline bool operator()(::bliss::debruijn::chain::compacted_chain_node<KMER> const & x, ::bliss::debruijn::chain::compacted_chain_node<KMER> const & y) {
              int8_t cmp = std::get<0>(x).compare(std::get<0>(y));
              return (cmp < 0) ? true :
                  ((cmp == 0) ? (std::get<1>(x) < std::get<1>(y)) : false);
            }
        };

        template <typename KMER>
        struct print_chain_node {
            inline void operator()(::bliss::debruijn::chain::compacted_chain_node<KMER> const & x) {
              if (std::get<1>(x) == 0) printf("\n[CHAIN] %s", bliss::utils::KmerUtils::toASCIIString(std::get<0>(x)).c_str());
              else printf("%c", KMER::KmerAlphabet::TO_ASCII[std::get<2>(x)]);
              //else printf("\n%d %c", std::get<1>(x), KMER::KmerAlphabet::TO_ASCII[std::get<2>(x)]);
            }
        };


      } // namespace chain

    } //namespace operation


    // for lex_less with kmer edges.
    namespace transform {
      /// standard companion function (used by lex_less) to compact_simple_biedge to get reverse complement of edge.
      template <typename KMER>
      inline ::bliss::debruijn::operation::chain::terminus_update_md<KMER>
      reverse_complement(::bliss::debruijn::operation::chain::terminus_update_md<KMER> const & x) {
        return ::bliss::debruijn::operation::chain::terminus_update_md<KMER>(x.first.reverse_complement(), -(x.second));
      }


      template <typename KMER>
      inline ::bliss::debruijn::operation::chain::chain_update_md<KMER>
      reverse_complement(::bliss::debruijn::operation::chain::chain_update_md<KMER> const & x) {
        return ::bliss::debruijn::operation::chain::chain_update_md<KMER>(std::get<0>(x).reverse_complement(), std::get<1>(x), -(std::get<2>(x)));
      }
    }  // ns transform



  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_CHAIN_OPERATIONS_HPP_
