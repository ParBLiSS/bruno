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
 * debruijn_graph_loader.hpp
 *
 *  Created on: Aug 17, 2015
 *      Author: tony pan
 *      Author: yongchao
 */

#ifndef DEBRUIJN_GRAPH_LOADER_HPP_
#define DEBRUIJN_GRAPH_LOADER_HPP_

#include "bliss-config.hpp"

#include <utility>      // pair and utility functions.
#include <type_traits>

#include "utils/logging.h"

#include "common/alphabets.hpp"
#include "common/sequence.hpp"

#include "utils/kmer_utils.hpp"

#include "iterators/transform_iterator.hpp"
#include "iterators/zip_iterator.hpp"

#include "io/sequence_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "index/quality_score_iterator.hpp"
#include "debruijn/edge_iterator.hpp"

#include "utils/benchmark_utils.hpp"

namespace bliss
{
  namespace debruijn
  {
    //    	// generic class to create a directed de Bruijn graph
    // reuses the kmer index class
    // only need to supply correct de bruijn parser.

    /**
     * @brief  transform from k+2-mer to <kmer, compact_simple_biedge>
     *
     */
    template <typename KmerType>
    struct k2mer_to_edge {
        using Alphabet = typename KmerType::KmerAlphabet;

        static_assert(::std::is_same<Alphabet, ::bliss::common::DNA16>::value ||
                      ::std::is_same<Alphabet, ::bliss::common::DNA6>::value ||
                       ::std::is_same<Alphabet, ::bliss::common::RNA6>::value,
                        "Only alphabets with explicit gap character is supported in order to represent empty edges.");

        using EdgeAlphabet = ::bliss::common::DNA16;
        using edge_type = ::bliss::debruijn::compact_simple_biedge;

        using K2merType = ::bliss::common::Kmer<KmerType::size + 2, Alphabet, typename KmerType::KmerWordType>;

        ::std::pair<KmerType, edge_type> operator()(K2merType const & k2mer) {
          // make the k-mer from k+2-mer  - center k characters
          K2merType tmp = k2mer >> 1;
          KmerType kmer(tmp.getData());  // copy the content.  sanitize is called by constructor.

          //std::cout << "k2mer " << bliss::utils::KmerUtils::toASCIIString(k2mer) << std::endl;

          // extract the biedge.
          edge_type edge;
          unsigned char right = k2mer.getData()[0] & ((0x1 << K2merType::bitsPerChar) - 1); // right char
          unsigned char left  = (k2mer >> (K2merType::size - 1)).getData()[0] & ((0x1 << K2merType::bitsPerChar) - 1); // left char

          edge.getDataRef()[0] = (EdgeAlphabet::FROM_ASCII[Alphabet::TO_ASCII[left]] << 4) |
              EdgeAlphabet::FROM_ASCII[Alphabet::TO_ASCII[right]];

          return std::make_pair(kmer, edge);
        }
    };

    template <typename KmerType>
    struct first_k2mer_to_edge {
        using Alphabet = typename KmerType::KmerAlphabet;

        static_assert(::std::is_same<Alphabet, ::bliss::common::DNA16>::value ||
                      ::std::is_same<Alphabet, ::bliss::common::DNA6>::value ||
                       ::std::is_same<Alphabet, ::bliss::common::RNA6>::value,
                        "Only alphabets with explicit gap character is supported in order to represent empty edges.");

        using EdgeAlphabet = ::bliss::common::DNA16;
        using edge_type = ::bliss::debruijn::compact_simple_biedge;

        using K2merType = ::bliss::common::Kmer<KmerType::size + 2, Alphabet, typename KmerType::KmerWordType>;

        ::std::pair<KmerType, edge_type> operator()(K2merType const & k2mer) {
          // make the k-mer from first k of k+2-mer
          K2merType tmp = k2mer >> 2;
          KmerType kmer(tmp.getData());  // copy the content.  sanitize is called by constructor.

          //std::cout << "k2mer " << bliss::utils::KmerUtils::toASCIIString(k2mer) << std::endl;

          // extract the biedge.
          edge_type edge;
          unsigned char right = ((k2mer.getData()[0]) >> K2merType::bitsPerChar) & ((0x1 << K2merType::bitsPerChar) - 1); // right char
          // no left

          edge.getDataRef()[0] = EdgeAlphabet::FROM_ASCII[Alphabet::TO_ASCII[right]];

          return std::make_pair(kmer, edge);
        }
    };


    template <typename KmerType>
    struct last_k2mer_to_edge {
        using Alphabet = typename KmerType::KmerAlphabet;

        static_assert(::std::is_same<Alphabet, ::bliss::common::DNA16>::value ||
                      ::std::is_same<Alphabet, ::bliss::common::DNA6>::value ||
                       ::std::is_same<Alphabet, ::bliss::common::RNA6>::value,
                        "Only alphabets with explicit gap character is supported in order to represent empty edges.");

        using EdgeAlphabet = ::bliss::common::DNA16;
        using edge_type = ::bliss::debruijn::compact_simple_biedge;

        using K2merType = ::bliss::common::Kmer<KmerType::size + 2, Alphabet, typename KmerType::KmerWordType>;

        ::std::pair<KmerType, edge_type> operator()(K2merType const & k2mer) {
          // make the k-mer from last k of k+2-mer
          KmerType kmer(k2mer.getData());  // copy the content.  sanitize is called by constructor.

          //std::cout << "k2mer " << bliss::utils::KmerUtils::toASCIIString(k2mer) << std::endl;

          // extract the biedge.
          edge_type edge;
          // no right
          unsigned char left  = (k2mer >> (K2merType::size - 2)).getData()[0] & ((0x1 << K2merType::bitsPerChar) - 1); // left char

          edge.getDataRef()[0] = (EdgeAlphabet::FROM_ASCII[Alphabet::TO_ASCII[left]] << 4);

          return std::make_pair(kmer, edge);
        }
    };



    /**
     * @brief  generate de bruijn graph nodes (k-mers) and edges (compacted left and right edges) from reads
     * @details  internally, read a k+2-mer.  if k+2-mer is the beginning or end, generate additional shifted version.
     *          convert the k+2-mer into <kmer, compacted_edge> pairs, then insert into output vector.
     * @tparam TupleType  value type of outputIt, or the generated type.  supplied so that we can use template template param with Index.
     */
    template <typename KmerType>
    struct debruijn_graph_parser {

        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUT.

        k2mer_to_edge<KmerType> transformer;
        using K2merType = typename k2mer_to_edge<KmerType>::K2merType;
        using Alphabet = typename KmerType::KmerAlphabet;
        using value_type = std::pair<KmerType, typename k2mer_to_edge<KmerType>::edge_type >;

        template <typename SeqType, typename OutputIt>
        OutputIt operator()(SeqType & read, OutputIt output_iter) {
          static_assert(std::is_same<value_type, typename ::std::iterator_traits<OutputIt>::value_type>::value,
                        "output type and output container value type are not the same");

          // filter out EOL characters
          using CharIter = bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>;
          // converter from ascii to alphabet values
          using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
          // kmer generation iterator
          using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, K2merType>;

          bliss::index::kmer::NotEOL neol;
          CharIter char_start(neol, read.seq_begin, read.seq_end);
          CharIter char_end(neol, read.seq_end);

          // count number of valid characters.
          size_t i = 0;
          for (CharIter char_iter = char_start; (char_iter != char_end) && (i < K2merType::size); ++char_iter, ++i) {};


          // if too few characters, do nothing
          if (i < KmerType::size) return output_iter;

          // have some stuff to process.
          BaseCharIterator base_start(char_start, bliss::common::ASCII2<Alphabet>());
          BaseCharIterator base_end(char_end, bliss::common::ASCII2<Alphabet>());

          // if exactly k characters.
          // note that if record is split (i.e. FASTA), there is overlap of at least k+1, so we should not encounter this.
          if (i == KmerType::size) {
            // make the kmer
            KmerType kmer;
            for (; base_start != base_end; ++base_start) {
              kmer.nextFromChar(*base_start);
            }

            // no edge.
            *output_iter = std::make_pair(kmer, typename k2mer_to_edge<KmerType>::edge_type());
            ++output_iter;

            return output_iter;
          }

          // if exactly k+1 characters
          if (i == (KmerType::size + 1)) {
            // create the k+2mer
            K2merType k2mer;
            for (; base_start != base_end; ++base_start) {
              k2mer.nextFromChar(*base_start);
            }

            // transform (left edge empty)
            *output_iter = transformer(k2mer);
            ++output_iter;

            // shift and transform (right edge empty)
            *output_iter = transformer(k2mer << 1);
            ++output_iter;

            return output_iter;
          }


          // else >= k+2 characters.

          //== set up the kmer generating iterators.
          KmerIter start(base_start, true);
          KmerIter end(base_end, false);

          // process the first if the seq_begin is at the beginning of the record.
          if (! read.is_record_truncated_at_begin()) {
            *output_iter = transformer((*start) >> 1);
            ++output_iter;
          }

          K2merType k2mer;  // temp storage  (for the last k2mer to be processed)
          // now transform and copy all.  can't use std::transform because saving k2mer, which is probably cheaper than copying KmerIter.
          for (; start != end; ++start) {
            k2mer = *start;
            *output_iter = transformer(k2mer);
            ++output_iter;
          }

          // if seq_end is at then end of the record,
          if (! read.is_record_truncated_at_end()) {
            *output_iter = transformer(k2mer << 1);
            ++output_iter;
          }

          return output_iter;

        }
    };

    /**
     * @brief 				parse the first k-mer of a read.
     * @tparam KmerType       output value type of this parser.  not necessarily the same as the map's final storage type.
     */
    template <typename KmerType>
    struct FirstKmerParser {

        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
    	using value_type = KmerType;
    	using kmer_type = KmerType;

      /**
       * @brief generate kmers from 1 sequence.  result inserted into output_iter, which may be preallocated.
       * @param read          sequence object, which has pointers to the raw byte array.
       * @param output_iter   output iterator pointing to insertion point for underlying container.
       * @return new position for output_iter
       * @tparam SeqType      type of sequence.  inferred.
       * @tparam OutputIt     output iterator type, inferred.
       */
    	template <typename SeqType, typename OutputIt>
    	OutputIt operator()(SeqType & read, OutputIt output_iter) {

    		static_assert(std::is_same<KmerType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
    						"output type and output container value type are not the same");

    		using Alphabet = typename KmerType::KmerAlphabet;

    		/// converter from ascii to alphabet values
    		using BaseCharIterator = bliss::iterator::transform_iterator<bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>,
    				bliss::common::ASCII2<Alphabet> >;

    		/// kmer generation iterator
    		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

    		static_assert(std::is_same<typename std::iterator_traits<KmerIterType>::value_type,
    				KmerType>::value,
    				"input iterator and output iterator's value types differ");

    		//== filtering iterator
    		bliss::index::kmer::NotEOL neol;
    		bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType> eolstart(neol, read.seq_begin, read.seq_end);
    		bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType> eolend(neol, read.seq_end);

    		//== set up the kmer generating iterators.
    		KmerIterType start(BaseCharIterator(eolstart, bliss::common::ASCII2<Alphabet>()), true);
    		KmerIterType end(BaseCharIterator(eolend, bliss::common::ASCII2<Alphabet>()), false);

    //    printf("First: pos %lu kmer %s\n", read.id.id, bliss::utils::KmerUtils::toASCIIString(*start).c_str());
    		if (start != end) {
    			*output_iter = *start;
    			++output_iter;
    		}

    		return output_iter;

    	}
    };

    template <typename KmerType>
    struct LastKmerParser {

        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
    	using value_type = KmerType;
    	using kmer_type = KmerType;

      /**
       * @brief generate kmers from 1 sequence.  result inserted into output_iter, which may be preallocated.
       * @param read          sequence object, which has pointers to the raw byte array.
       * @param output_iter   output iterator pointing to insertion point for underlying container.
       * @return new position for output_iter
       * @tparam SeqType      type of sequence.  inferred.
       * @tparam OutputIt     output iterator type, inferred.
       */
    	template <typename SeqType, typename OutputIt>
    	OutputIt operator()(SeqType & read, OutputIt output_iter) {

    		static_assert(std::is_same<KmerType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
    						"output type and output container value type are not the same");

    		using Alphabet = typename KmerType::KmerAlphabet;

    		/// converter from ascii to alphabet values
    		using BaseCharIterator = bliss::iterator::transform_iterator<bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>,
    				bliss::common::ASCII2<Alphabet> >;

    		/// kmer generation iterator
    		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

    		static_assert(std::is_same<typename std::iterator_traits<KmerIterType>::value_type,
    				KmerType>::value,
    				"input iterator and output iterator's value types differ");

    		// then compute and store into index (this will generate kmers and insert into index)

    		//== filtering iterator
    		bliss::index::kmer::NotEOL neol;
    		bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType> eolstart(neol, read.seq_begin, read.seq_end);
    		bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType> eolend(neol, read.seq_end);

    		size_t count = std::distance(eolstart, eolend);
    		if (count < KmerType::size) return output_iter;

    		//== move to last.
    		std::advance(eolstart, count - KmerType::size );
    		assert(eolstart != eolend);

    		//== set up the kmer generating iterators.
    		KmerIterType start(BaseCharIterator(eolstart, bliss::common::ASCII2<Alphabet>()), true);
    		KmerIterType end(BaseCharIterator(eolend, bliss::common::ASCII2<Alphabet>()), false);

    //    printf("First: pos %lu kmer %s\n", read.id.id, bliss::utils::KmerUtils::toASCIIString(*start).c_str());
    		assert(start != end);

			*output_iter = *start;
			++output_iter;

    		return output_iter;

    	}
    };


//    template <typename KmerType>
//    struct debruijn_graph_parser {
//
//        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUT.
//        using edge_type = ::bliss::debruijn::compact_simple_biedge;
//        using value_type = std::pair<KmerType, edge_type>;
//
//        using Alphabet = typename KmerType::KmerAlphabet;
//
//        template <typename SeqType, typename OutputIt>
//        OutputIt operator()(SeqType & read, OutputIt output_iter) {
//          static_assert(std::is_same<value_type, typename ::std::iterator_traits<OutputIt>::value_type>::value,
//                        "output type and output container value type are not the same");
//
//          // filter out EOL characters
//          using CharIter = bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>;
//          // converter from ascii to alphabet values
//          using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
//          // kmer generation iterator
//          using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
//
//          // handles the ascii to DNA16 mapping internally.
//          using EdgeIterType = bliss::debruijn::iterator::biedge_generating_iterator<CharIter, ::bliss::common::DNA16>;
//
//
//          // combine kmer iterator and position iterator to create an index iterator type.
//          using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, EdgeIterType>;
//
//          static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
//                        value_type>::value,
//                        "Generating iterator and output container value types differ");
//
//          // then compute and store into index (this will generate kmers and insert into index)
//          if (read.seq_begin == read.seq_end) return output_iter;
//
//          //== set up the kmer generating iterators.
//          bliss::index::kmer::NotEOL neol;
//          KmerIter start(BaseCharIterator(CharIter(neol, read.seq_begin, read.seq_end), bliss::common::ASCII2<Alphabet>()), true);
//          KmerIter end(BaseCharIterator(CharIter(neol, read.seq_end), bliss::common::ASCII2<Alphabet>()), false);
//
//
//          // set up edge iterator.  returns the left and right chars around the kmer IN THE READ.
//          EdgeIterType edge_start(CharIter(neol, read.seq_begin, read.seq_end), CharIter(neol, read.seq_end), KmerType::size);
//          EdgeIterType edge_end (CharIter(neol, read.seq_end));
//
//
//          // ==== set up the zip iterators
//          KmerIndexIterType node_start (start, edge_start);
//          KmerIndexIterType node_end(end, edge_end);
//
//          return ::std::copy(node_start, node_end, output_iter);
//
//        }
//    };


    // deprecated for now until this can be made compatible with the k+2-mer approach.

//    /*generate de Brujin graph nodes and edges, which each node associated with base quality scores*/
//    template <typename KmerType, typename QualType=double,
//        template<typename> class QualityEncoder = bliss::index::Illumina18QualityScoreCodec>
//    struct debruijn_graph_quality_parser {
//
//        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUTIT.
//        using edge_type = ::bliss::debruijn::compact_simple_biedge;
//        using value_type = std::pair<KmerType, std::pair<edge_type, QualType> >;
//
//        using Alphabet = typename KmerType::KmerAlphabet;
//
//        template <typename SeqType, typename OutputIt>
//        OutputIt operator()(SeqType & read, OutputIt output_iter) {
//
//          static_assert(std::is_same<value_type, typename ::std::iterator_traits<OutputIt>::value_type>::value,
//                        "output type and output container value type are not the same");
//
//
//          // filter out EOL characters
//          using CharIter = bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>;
//          // converter from ascii to alphabet values
//          using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
//          // kmer generation iterator
//          using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
//
//          // handles the ascii to DNA16 mapping internally.
//          using EdgeIterType = bliss::debruijn::iterator::biedge_generating_iterator<CharIter, bliss::common::DNA16>;
//
//          // also remove eol from quality score
//          using QualIterType =
//              bliss::index::QualityScoreGenerationIterator<bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>, KmerType::size, QualityEncoder<QualType> >;
//
//          // combine kmer iterator and position iterator to create an index iterator type.
//          using KmerInfoIterType = bliss::iterator::ZipIterator<EdgeIterType, QualIterType>;
//
//          using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, KmerInfoIterType>;
//
//          static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
//                        value_type>::value,
//                        "generating iterator and output iterator's value types differ");
//
//
//          // then compute and store into index (this will generate kmers and insert into index)
//          if (read.seq_begin == read.seq_end || read.qual_begin == read.qual_end) return output_iter;
//
//          //== set up the kmer generating iterators.
//          bliss::index::kmer::NotEOL neol;
//          KmerIter start(BaseCharIterator(CharIter(neol, read.seq_begin, read.seq_end), bliss::common::ASCII2<Alphabet>()), true);
//          KmerIter end(BaseCharIterator(CharIter(neol, read.seq_end), bliss::common::ASCII2<Alphabet>()), false);
//
//
//          // set up edge iterator.  returns the left and right chars around the kmer IN THE READ.
//          EdgeIterType edge_start(CharIter(neol, read.seq_begin, read.seq_end), CharIter(neol, read.seq_end), KmerType::size);
//          EdgeIterType edge_end (CharIter(neol, read.seq_end));
//
//
//          QualIterType qual_start(CharIter(neol, read.qual_begin, read.qual_end));
//          QualIterType qual_end(CharIter(neol, read.qual_end));
//
//          KmerInfoIterType info_start(edge_start, qual_start);
//          KmerInfoIterType info_end(edge_end, qual_end);
//
//
//
//          // ==== set up the zip iterators
//          KmerIndexIterType node_start(start, info_start);
//          KmerIndexIterType node_end(end, info_end);
//
//          return ::std::copy(node_start, node_end, output_iter);
//
//        }
//    };


  } /* namespace debruijn */
} /* namespace bliss */

#endif /* DEBRUIJN_GRAPH_LOADER_HPP_ */
