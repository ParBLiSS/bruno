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
 * debruijn_biedge_loader.hpp
 *
 *  Created on: Aug 4, 2015
 *      Author: yongchao
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_BIEDGE_LOADER_HPP_
#define DEBRUIJN_BIEDGE_LOADER_HPP_

#include <iterator>

// #include "common/alphabets.hpp"
//#include "common/kmer.hpp"
//#include "io/kmer_parser.hpp"
//#include "utils/filter_utils.hpp"

namespace bliss
{

  namespace debruijn
  {


    namespace biedge
    {
      /// simple biedge type. packed in and out edge character (adjacent to current k-mer.)  e.g. x-Kmer-y would be have xy as the biedge representation.  simple because it's not compound.
      /// NOTE: simple biedge encodes 2 characters from the alphabet, not the same as compact biedge.
      ///   uses DNA alphabet and is the most informative of the alphabet.
      using compact_simple_biedge = bliss::common::Kmer<2, ::bliss::common::DNA16, uint8_t>;

      // sanity check.
      static_assert(sizeof(compact_simple_biedge) == 1, "size of compact simple biedge is larger than 1 byte.");

      /// Kmer representing one edge between 2 nodes.
      template <typename KmerType>
      using kmer_simple_edge = bliss::common::Kmer<KmerType::size + 1, typename KmerType::KmerAlphabet, typename KmerType::KmerWordType>;

      /// node from Kmer and biedge
      template <typename KmerType>
      using compact_simple_biedge_kmer_node = ::std::pair<KmerType, compact_simple_biedge>;


      // TODO: turn these into class member functions.

      /**
       * @brief function to get IN edge source kmer given central kmer and biedge.
       * @note  converts edge character from DNA16 to that of the KmerType.  there is potentially information loss if KmerType has an
       *        alphabet with fewer characters.
       */
      template <typename KmerType>
      KmerType get_in_kmer(KmerType const & kmer, compact_simple_biedge const & edge) {
        KmerType kk(kmer);

        using Alpha = typename KmerType::KmerAlphabet;

        kk.nextReverseFromChar(Alpha::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[(edge.getData()[0] >> (compact_simple_biedge::bitsPerChar)) ] ]);

        return kk;
      }
      template <typename KmerType>
      KmerType get_in_kmer(compact_simple_biedge_kmer_node<KmerType> const & edge) {
        return get_in_kmer(edge.first, edge.second);
      }

      /**
       * @brief function to get OUT edge source kmer given central kmer and biedge.
       * @note  converts edge character from DNA16 to that of the KmerType.  there is potentially information loss if KmerType has an
       *        alphabet with fewer characters.
       */
      template <typename KmerType>
      KmerType get_out_kmer(KmerType const & kmer, compact_simple_biedge const & edge) {
        KmerType kk(kmer);

        using Alpha = typename KmerType::KmerAlphabet;

        kk.nextFromChar(Alpha::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[edge.getData()[0] & 0xF]]);

        return kk;
      }
      template <typename KmerType>
      KmerType get_out_kmer(compact_simple_biedge_kmer_node<KmerType> const & edge) {
        return get_out_kmer(edge.first, edge.second);
      }


      /**
       * @brief function to get IN edge k+1 mer given central kmer and biedge.
       * @note  converts edge character from DNA16 to that of the KmerType.  there is potentially information loss if KmerType has an
       *        alphabet with fewer characters.
       */
      template <typename KmerType>
      kmer_simple_edge<KmerType> get_in_edge_k1mer(KmerType const & kmer, compact_simple_biedge const & edge) {

        using Alpha = typename KmerType::KmerAlphabet;

        kmer_simple_edge<KmerType> kk;
        std::memcpy(kk.getDataRef(), kmer.getData(), KmerType::nWords * sizeof(typename KmerType::KmerWordType));

        kmer_simple_edge<KmerType> kk2;
        kk2.nextFromChar(Alpha::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[(edge.getData()[0] >> (compact_simple_biedge::bitsPerChar))]]);
        kk2 <<= KmerType::size;   // shift by k

        return kk | kk2;
      }
      template <typename KmerType>
      kmer_simple_edge<KmerType> get_in_edge_k1mer(compact_simple_biedge_kmer_node<KmerType> const & edge) {
        return get_in_edge_k1mer(edge.first, edge.second);
      }


      /**
       * @brief function to get OUT edge k+1 mer given central kmer and biedge.
       * @note  converts edge character from DNA16 to that of the KmerType.  there is potentially information loss if KmerType has an
       *        alphabet with fewer characters.
       */
      template <typename KmerType>
      kmer_simple_edge<KmerType> get_out_edge_k1mer(KmerType const & kmer, compact_simple_biedge const & edge) {

        using Alpha = typename KmerType::KmerAlphabet;

        kmer_simple_edge<KmerType> kk;
        std::memcpy(kk.getDataRef(), kmer.getData(), KmerType::nWords * sizeof(typename KmerType::KmerWordType));

        kk.nextFromChar(Alpha::FROM_ASCII[::bliss::common::DNA16::TO_ASCII[edge.getData()[0] & 0xF]]);

        return kk;
      }
      template <typename KmerType>
      kmer_simple_edge<KmerType> get_out_edge_k1mer(compact_simple_biedge_kmer_node<KmerType> const & edge) {
        return get_out_edge_k1mer(edge.first, edge.second);
      }


    // for doing reverse complement of the biedge.
    namespace transform {

      /**
       * @brief extract the biedges (without corresponding kmers), if kmer ordering is deterministic and implicit (e.g. ordered by generation from file)
       * @details   the biedges can now be used to represent presence and absence of edges.
       *            in fact, once kmers are read, the edges can be constructed directly.
       *
       *            USED BY std::transform, unary operator form
       */
      struct extract_biedge {

          template <typename KmerType>
          compact_simple_biedge operator()(compact_simple_biedge_kmer_node<KmerType> const & edge) {
            return edge.second;
          }

      };


      /**
       * @brief given kmer vector and biedge vector, merge them into single node of kmer (center) + biedge.
       * @details  construct kmer + 2 edges from previously saved biedge and newly reconstructed kmers.
       *
       *            assumes ordering of kmers and biedge are the same (e.g. order of construction from file.).
       *
       *          USED by std::transform, binary operator form
       */
      struct merge_kmer_biedge {

          template <typename KmerType>
          compact_simple_biedge_kmer_node<KmerType> operator()(KmerType const & kmer, compact_simple_biedge const & edge) {
            return std::make_pair(kmer, edge);
          }

      };

    }




    namespace iterator
    {

      // CURRENTLY NOT USED, but may want to resurrect.

      // careful with the use of enable_if.  evaluation should occur at function call time,
      //   i.e. class template params will be evaluated with no substitution.
      // instead, a function should declare a template parameter proxy for the class template parameter.
      //   then enable_if evaluates using the proxy.
      //  e.g. template< class c = C; typename std::enable_if<std::is_same<c, std::string>::value, int>::type x = 0>

        template<typename IT, typename Alphabet = bliss::common::DNA16>
        class biedge_generating_iterator;

      /**
       * @class   biedge_generating_iterator
       * @brief   given a k-mer position, retrieve its left and right bases, including the dummy bases at both ends of read
       * @details specializations of this class uses a byte to manage the edge info.
       *          upper 4 bits holds the left base (encoded.  correspond to in edge)
       *          lower 4 bits holds the right base (encoded, correspond to out edge)
       *
       *          no reverse complement or reordering is applied.
       *
       *          edge iterator should be valid for std::distance(_data_end - _data_start - k + 1) iterations.
       *
       *          using DNA16 because it is the most encompassing, also it supports 'N'
       *
       *          directly generates left and right edges
       *
       *          NOTE: simple biedge encodes 2 characters from the alphabet, not the same as compact biedge.
       *
       *          TODO: rework this - if try to get _left without calling ++ first, fails.
       */
      template<typename IT>
      class biedge_generating_iterator<IT, bliss::common::DNA16> :
        public ::std::iterator<::std::forward_iterator_tag, ::bliss::debruijn::biedge::compact_simple_biedge >
      {
        protected:

            // curr position = last char in k-mer
            IT _curr;
            //previous position = char before k-mer
            IT _left;
            //next position, char after k-mer
            IT _right;

            /*data*/
            mutable IT _data_start;
            mutable IT _data_end;

          public:

            using Alphabet = bliss::common::DNA16;
            using self_type = biedge_generating_iterator<IT, Alphabet>; /*define edge iterator type*/
            using edge_type = ::bliss::debruijn::biedge::compact_simple_biedge; //type to represent an edge

            // accessors
            IT const & getBaseIterator() const { return _curr; }
            IT const & getLeftIterator() const { return (_left == _data_end ? _data_start : _left); }
            IT const & getRightIterator() const { return _right; }

            //constructor
            biedge_generating_iterator(IT data_start, IT data_end, const uint32_t k)
              : _curr (data_start), _left(data_end), _right(data_start), _data_start(data_start), _data_end(data_end)
            {
              /*compute the offset*/
              ::std::advance(_curr, k - 1);
              _right = _curr;
              ::std::advance(_right, 1);
            }
            biedge_generating_iterator(IT data_end)
              : _curr(data_end), _left(data_end), _right(data_end), _data_start(data_end), _data_end(data_end)
            {}

            /// copy constructor
            biedge_generating_iterator(const self_type& Other)
              : _curr (Other._curr), _left(Other._left), _right(Other._right),
                _data_start(Other._data_start), _data_end(Other._data_end)
            { /*do nothing*/ }


            /// copy assignment iterator
            self_type& operator=(const self_type& Other)
            {
              _curr = Other._curr;
              _left = Other._left;
              _right = Other._right;
              _data_start = Other._data_start;
              _data_end = Other._data_end;

              return *this;
            }

            /// move constructor
            biedge_generating_iterator(self_type&& Other)
              : _curr(std::move(Other._curr)),
				_left(std::move(Other._left)),
				_right(std::move(Other._right)),
                _data_start(std::move(Other._data_start)),
				_data_end(std::move(Other._data_end))
            { /*do nothing*/ }


            /// move assignment iterator
            self_type& operator=(self_type&& Other)
            {
              _curr = std::move(Other._curr);
              _left = std::move(Other._left);
              _right = std::move(Other._right);
              _data_start = std::move(Other._data_start);
              _data_end = std::move(Other._data_end);

              return *this;
            }

            biedge_generating_iterator() : _data_start(), _data_end(_data_start) {};

            /// increment to next matching element in base iterator
            self_type& operator++()
            {  // if _curr at end, subsequent calls should not move _curr.
               // on call, if not at end, need to move first then evaluate.
              if (_curr == _data_end){  // if at end, don't move it.
                return *this;
              }

              /*save the previous position*/
              if (_left == _data_end) _left = _data_start;
              else ++_left;

              /*move forward by 1*/
              ++_curr;

              /*ensure that _right does not exceed _end*/
              if(_right != _data_end){
                ++_right;
              }
              return *this;
            }

            /**
             * post increment.  make a copy then increment that.
             */
            self_type operator++(int)
            {
              self_type output(*this);
              this->operator++();
              return output;
            }

            /// compare 2 filter iterators
            inline bool operator==(const self_type& rhs)
            {
              return _curr == rhs._curr;
            }

            /// compare 2 filter iterators
            inline bool operator!=(const self_type& rhs)
            {
              return _curr != rhs._curr;
            }

            /// dereference operator. _curr is guaranteed to be valid.  out edge is at MSB,
            inline edge_type operator*() const
            {
              edge_type edges = edge_type();

              /*using four bits to represent an edge*/
              if(_left != _data_end && _right != _data_end){
                /*internal k-mer node*/
                edges.getDataRef()[0] = (Alphabet::FROM_ASCII[*_left] << 4) | Alphabet::FROM_ASCII[*_right];
              } else if(_left == _data_end && _right != _data_end){  /*the left-most k-mer node*/
                edges.getDataRef()[0] = Alphabet::FROM_ASCII[*_right];
              } else if(_left != _data_end && _right == _data_end){  /*the right-most k-mer node*/
                edges.getDataRef()[0] = Alphabet::FROM_ASCII[*_left] << 4;
              }

              /*if(_left == _end && _right == _end)*/
              return edges;
            }
        };

      template<typename IT>
      using DNA16_biedge_generating_iterator = biedge_generating_iterator<IT, bliss::common::DNA16>;


    } // iterator


    namespace io {

        /**
         * @brief  generate kmer and biedge pair (simple biedge.)
         * @details  performs the same function as debruijn_graph_parser, but does not start with k+2 mer then transform to kmer + simple biedge.
         *           doing it via k2mer requires shifts at beginning and end of the sequence fragment.  this approach does not require that as a separate step.
         *
         *           some questions:  if input is less than k, what happens?  start and end should then be the same and no output should be generated.
         *                            if there is overlap, does to stop the parsing at the last position of the overlap?
         *                              need info about the valid range of the block being parsed.
         */
        template <typename KmerType>
        class debruijn_kmer_simple_biedge_parser {

        public:

          /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUT.
          using edge_type = ::bliss::debruijn::biedge::compact_simple_biedge;
          using value_type = std::pair<KmerType, edge_type>;
          using kmer_type = KmerType;
          static constexpr size_t window_size = kmer_type::size + 2;

        protected:
          using Alphabet = typename KmerType::KmerAlphabet;


          // filter out EOL characters
          template <typename SeqType>
          using CharIter = bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>;

          // converter from ascii to alphabet values
          template <typename SeqType>
          using BaseCharIterator = bliss::iterator::transform_iterator<CharIter<SeqType>, bliss::common::ASCII2<Alphabet> >;

          // kmer generation iterator
          template <typename SeqType>
          using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator<SeqType>, KmerType>;

          // handles the ascii to DNA16 mapping internally.
          template <typename SeqType>
          using EdgeIterType = bliss::debruijn::biedge::iterator::biedge_generating_iterator<CharIter<SeqType>, ::bliss::common::DNA16>;

          ::bliss::partition::range<size_t> valid_range;



        public:

          // combine kmer iterator and position iterator to create an index iterator type.
          template <typename SeqType>
          using iterator_type = bliss::iterator::ZipIterator<KmerIter<SeqType>, EdgeIterType<SeqType> >;

          debruijn_kmer_simple_biedge_parser(::bliss::partition::range<size_t> const & _valid_range) : valid_range(_valid_range) {};

          // default constructor for purpose of copying
          debruijn_kmer_simple_biedge_parser() {};

          debruijn_kmer_simple_biedge_parser(debruijn_kmer_simple_biedge_parser const & other) : valid_range(other.valid_range) {};

          debruijn_kmer_simple_biedge_parser(debruijn_kmer_simple_biedge_parser && other) : valid_range(std::move(other.valid_range)) {};

          debruijn_kmer_simple_biedge_parser& operator=(debruijn_kmer_simple_biedge_parser const & other) {
        	  valid_range = other.valid_range;
            return *this;
          }

          debruijn_kmer_simple_biedge_parser& operator=(debruijn_kmer_simple_biedge_parser && other) {
        	  valid_range = std::move(other.valid_range);
            return *this;
          }


          template <typename SeqType>
          iterator_type<SeqType> begin(SeqType const & read, size_t const & window = window_size) const {
              static_assert(std::is_same<typename std::iterator_traits<iterator_type<SeqType> >::value_type,
                            value_type>::value,
                            "Generating iterator value type differs from expected");

              typename SeqType::IteratorType seq_begin;
              typename SeqType::IteratorType seq_end;
              bool longer_than_window = false;

              std::tie(seq_begin, seq_end, longer_than_window) =
            		  ::bliss::index::kmer::KmerParser<KmerType>::get_valid_iterator_range(read, valid_range, window);

              //== set up the kmer generating iterators.
              bliss::utils::file::NotEOL neol;


              if (longer_than_window) {
                  KmerIter<SeqType> start(BaseCharIterator<SeqType>(
                		  CharIter<SeqType>(neol, seq_begin, seq_end),
                		  bliss::common::ASCII2<Alphabet>()),
                		  true);

                  // set up edge iterator.  returns the left and right chars around the kmer IN THE READ.
                  EdgeIterType<SeqType> edge_start(
                		  CharIter<SeqType>(neol, seq_begin, seq_end),
    					  CharIter<SeqType>(neol, seq_end),
    					  KmerType::size);

                  // ==== set up the zip iterators
				  iterator_type<SeqType> node_start(start, edge_start);

				  // advance one if needed.
	//			  std::cout << "partition " << valid_range << " seq " << read << " truncated ? " << (read.is_record_truncated_at_begin() ? "y" : "n") << std::endl;
				  if (read.is_record_truncated_at_begin()) {
	//				  std::cout << " truncated.  move forward 1." << std::endl;
					  ++node_start;
				  }

				  return node_start;
              } else {
                  KmerIter<SeqType> end(BaseCharIterator<SeqType>(
                		  CharIter<SeqType>(neol, seq_end),
                		  bliss::common::ASCII2<Alphabet>()),
                		  false);

                  // set up edge iterator.  returns the left and right chars around the kmer IN THE READ.
                  EdgeIterType<SeqType> edge_end(CharIter<SeqType>(neol, seq_end));

				  // ==== set up the zip iterators
				  return iterator_type<SeqType>(end, edge_end);
              }
          }

          template <typename SeqType>
          iterator_type<SeqType> end(SeqType const & read, size_t const & window = window_size) const {
              static_assert(std::is_same<typename std::iterator_traits<iterator_type<SeqType> >::value_type,
                            value_type>::value,
                            "Generating iterator value type differs from expected");

              typename SeqType::IteratorType seq_begin;
              typename SeqType::IteratorType seq_end;
              bool longer_than_window = false;

              // if sequence is split across partitions, and this record is truncated at its end,
              // then we want to iterate up to end position in the overlap at window-1
              // The iterator can stop with kmer trailing end just inside the overlap region (valid end)
              // and left edge is at last valid region position while right edge is at the window-1 position.
              size_t seq_overlap_window = (read.is_record_truncated_at_end() ? window - 1 : window);

              std::tie(seq_begin, seq_end, longer_than_window) =
            		  ::bliss::index::kmer::KmerParser<KmerType>::get_valid_iterator_range(read, valid_range, seq_overlap_window);

              //== set up the kmer generating iterators.
              bliss::utils::file::NotEOL neol;

              KmerIter<SeqType> end(BaseCharIterator<SeqType>(
            		  CharIter<SeqType>(neol, seq_end),
            		  bliss::common::ASCII2<Alphabet>()),
            		  false);

              // set up edge iterator.  returns the left and right chars around the kmer IN THE READ.
              EdgeIterType<SeqType> edge_end(CharIter<SeqType>(neol, seq_end));


              // ==== set up the zip iterators
              return iterator_type<SeqType>(end, edge_end);
          }




          template <typename SeqType, typename OutputIt, typename Predicate = ::bliss::filter::TruePredicate>
          OutputIt operator()(SeqType const & read, OutputIt output_iter, Predicate const & pred = Predicate()) {

//            std::cout << " first raw dist = " << read.seq_size() << std::endl;

            static_assert(std::is_same<value_type, typename ::std::iterator_traits<OutputIt>::value_type>::value,
                          "output type and output container value type are not the same");

//            static_assert(std::is_same<typename std::iterator_traits<iterator_type<SeqType> >::value_type,
//                          value_type>::value,
//                          "Generating iterator and output container value types differ");
//
//
//            // then compute and store into index (this will generate kmers and insert into index)
//            if (read.seq_begin == read.seq_end) return output_iter;
//
//            //== set up the kmer generating iterators.
//            bliss::utils::file::NotEOL neol;
//            KmerIter<SeqType> start(BaseCharIterator<SeqType>(CharIter<SeqType>(neol, read.seq_begin, read.seq_end), bliss::common::ASCII2<Alphabet>()), true);
//            KmerIter<SeqType> end(BaseCharIterator<SeqType>(CharIter<SeqType>(neol, read.seq_end), bliss::common::ASCII2<Alphabet>()), false);
//
//
//            // set up edge iterator.  returns the left and right chars around the kmer IN THE READ.
//            EdgeIterType<SeqType> edge_start(CharIter<SeqType>(neol, read.seq_begin, read.seq_end),
//            		CharIter<SeqType>(neol, read.seq_end),
//					KmerType::size);
//            EdgeIterType<SeqType> edge_end(CharIter<SeqType>(neol, read.seq_end));
//
//
//            // ==== set up the zip iterators
//            iterator_type<SeqType> node_start(start, edge_start);
//            iterator_type<SeqType> node_end(end, edge_end);
//
////            std::cout << " raw dist = " << read.seq_size() << std::endl;
////            std::cout << " kmer iter dist = " << std::distance(start, end) << std::endl;
////            std::cout << " kmerbiedge iter dist = " << std::distance(node_start, node_end) << std::endl;
//
//            if (read.is_record_truncated_at_begin()) ++node_start;
//            // for truncated at end, we are already doing overlap of size k+1.
//
//            ::bliss::partition::range<size_t> seq_range(read.seq_global_offset(), read.seq_global_offset() + read.seq_size());
//            if (seq_range.contains(valid_range.end)) {
//              // seq_range contains overlap.
//
//                // not checking by end iterator at valid_range.end, since the NonEOLIter is a filter iterator that may skip over that pos.
//                int64_t valid_dist = valid_range.end - seq_range.start;
//
//              for (auto it = node_start; it != node_end; ++it, ++output_iter) {
//                // check tail of window -> transform iterator, get base -> non EOL iterator, get base -> seq raw char iter.
//                if (std::distance(read.seq_begin,
//                                  it.get_second_iterator().getLeftIterator().getBaseIterator()) >= valid_dist) {
//                  // NOTE: last position (valid_range.end) should be the trailing position of the kmer in the kmer-biedge,
//                  //       so that the incoming edge character is the last valid character.  This is to avoid double counting of edge at process boundary.
//                  //       in other words, if the node is considered a k2mer, then the last valid char should sit
//                  //         at trailing end of the k2mer at left of process boundary, at right side of boundary, 1 position is skipped.
//                  //       this does not apply to real sequence ends.
//                  //
//                  //       also note that since we use canonical, frequency accounting is affected by canonicalization.
//                  // when reconstructing node from kmer and edge vectors, the last kmer should sit just outside the valid range
//
//                  break;
//                }
//
//
//				*output_iter = *it;
//              }
//
//              return output_iter;
//
//            } else {
//            	return ::std::copy(node_start, node_end, output_iter);
//            }
            iterator_type<SeqType> istart = begin(read, window_size);
            iterator_type<SeqType> iend = end(read, window_size);

            return std::copy(istart, iend, output_iter);
          }
      };
        template <typename KmerType>
        constexpr size_t debruijn_kmer_simple_biedge_parser<KmerType>::window_size;


    // deprecated for now until this can be made compatible with the k+2-mer approach.

  //    /*generate de Brujin graph nodes and edges, which each node associated with base quality scores*/
  //    template <typename KmerType, typename QualType=double,
  //        template<typename> class QualityEncoder = bliss::index::Illumina18QualityScoreCodec>
  //    struct debruijn_graph_quality_parser {
  //
  //        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUTIT.
  //        using edge_type = ::bliss::debruijn::biedge::compact_simple_biedge;
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
  //          using EdgeIterType = bliss::debruijn::biedge::iterator::biedge_generating_iterator<CharIter, bliss::common::DNA16>;
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
      } // io
    } // biedge

    namespace transform
    {
      /// standard companion function (used by lex_less) to compact_simple_biedge to get reverse complement of edge.
      template <typename Kmer = ::bliss::debruijn::biedge::compact_simple_biedge>
      inline ::bliss::debruijn::biedge::compact_simple_biedge
      reverse_complement(::bliss::debruijn::biedge::compact_simple_biedge const & x) {
        return x.reverse_complement();
      }
    }


  } // debruijn
} // bliss



#endif /* DEBRUIJN_BIEDGE_LOADER_HPP_ */
