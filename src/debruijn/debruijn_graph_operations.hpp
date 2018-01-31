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
 * debruijn_graph_operations.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_GRAPH_OPERATIONS_HPP_
#define DEBRUIJN_GRAPH_OPERATIONS_HPP_

#include <tuple>
#include <utility>  // pair
#include "debruijn/debruijn_graph_node.hpp"

namespace bliss {
  namespace debruijn {

  // forward declare
  template <typename KMER>
      struct lex_less;

    namespace operation {

      namespace graph {

      /// update DBG node and edges.  effect is only entries already existing in the DBG are updated.
      /// return 1 if updated, 0 if not updated
      template <typename Kmer>
      struct compact_multi_biedge_update {

    	  template <typename Count>
          size_t operator()(::bliss::debruijn::graph::compact_multi_biedge<typename Kmer::KmerAlphabet, Count> & x,
        		  ::bliss::debruijn::biedge::compact_simple_biedge const & y) const {
        	  x.update(y);
        	  return 1;
          }
      };



        template <typename KMER>
        struct print_graph_node {
        	std::ostream & os;

        	print_graph_node(std::ostream & _os) : os(_os) {};

        	template <typename Alphabet, typename CountType, typename DUMMY>
            inline void operator()(std::pair<KMER, ::bliss::debruijn::graph::compact_multi_biedge<Alphabet, CountType, DUMMY> > const & x) {
            	os << bliss::utils::KmerUtils::toASCIIString(x.first) << "\t";

            	// N is not counted....
        		os << static_cast<size_t>(x.second.get_in_edge_frequency(Alphabet::FROM_ASCII['A'])) << "\t";
        		os << static_cast<size_t>(x.second.get_in_edge_frequency(Alphabet::FROM_ASCII['C'])) << "\t";
        		os << static_cast<size_t>(x.second.get_in_edge_frequency(Alphabet::FROM_ASCII['G'])) << "\t";
        		os << static_cast<size_t>(x.second.get_in_edge_frequency(Alphabet::FROM_ASCII['T'])) << "\t";
        		os << static_cast<size_t>(x.second.get_out_edge_frequency(Alphabet::FROM_ASCII['A'])) << "\t";
        		os << static_cast<size_t>(x.second.get_out_edge_frequency(Alphabet::FROM_ASCII['C'])) << "\t";
        		os << static_cast<size_t>(x.second.get_out_edge_frequency(Alphabet::FROM_ASCII['G'])) << "\t";
        		os << static_cast<size_t>(x.second.get_out_edge_frequency(Alphabet::FROM_ASCII['T'])) << "\t";
            	os << static_cast<size_t>(x.second.get_self_frequency()) << std::endl;
            }
        };

        template <typename KMER>
        struct print_graph_node_fasta {
          std::ostream & os;

          print_graph_node_fasta(std::ostream & _os) : os(_os) {};

          template <typename Alphabet, typename CountType, typename DUMMY>
            inline void operator()(std::pair<KMER, ::bliss::debruijn::graph::compact_multi_biedge<Alphabet, CountType, DUMMY> > const & x) {
              os << ">" << bliss::utils::KmerUtils::toASCIIString(x.first) << std::endl;
              os << bliss::utils::KmerUtils::toASCIIString(x.first) << std::endl;
            }
        };

      } // namespace graph

    } //namespace operation




  } //namespace debruijn
} //namespace bliss




#endif // DEBRUIJN_GRAPH_OPERATIONS_HPP_
