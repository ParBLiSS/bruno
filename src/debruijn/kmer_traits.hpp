/*
 * Copyright 2018 Georgia Institute of Technology
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
 * kmer_traits.hpp
 *
 *  Created on: Feb 5, 2018
 *      Author: tony pan
 *
 */

#ifndef KMER_TRAITS_HPP_
#define KMER_TRAITS_HPP_

#include "bliss-config.hpp"

// #include "utils/logging.h"
// #include "utils/transform_utils.hpp"

#include "common/kmer.hpp"

// this file is used for testing for palindromes.  k-palindrome can cause frequency to be split between two sides of k+2-mer,
// while k+1-mer may aggregate to one side if only k+1-mer is palindromic in k+2-mer, while k+2-mers with 2 k1 palindromes
// will count correctly...

namespace bliss
{
  namespace common
  {

    namespace kmer
    {
        template <typename KMER>
        struct kmer_traits;
        
        template <unsigned int KMER_SIZE, typename WORD_TYPE>
        struct kmer_traits<::bliss::common::Kmer<KMER_SIZE, ::bliss::common::DNA, WORD_TYPE> > {
            using KmerType = ::bliss::common::Kmer<KMER_SIZE, ::bliss::common::DNA, WORD_TYPE>;

            // check to see if kmer and its rc are identical
            static inline bool is_rc_palindrome(KmerType const & k) {
                if ((KMER_SIZE & 0x1) > 0) return false;  // odd DNA.  cannot be palindromic 
                return (k == k.reverse_complement());
            }
            // check to see if k+1-mer and its rc are identical.  append c to k on the right
            static inline bool is_k1_rc_palindrome(KmerType const & k, unsigned char const & r) {
                if ((KMER_SIZE & 0x1) == 0) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (k.getCharsAtPos(KMER_SIZE - 1, 1) != ::bliss::common::DNA::to_complement(r)) return false;   // two sides are not the same so not a rc_palindrome.
                ::bliss::common::Kmer<KMER_SIZE, ::bliss::common::DNA, WORD_TYPE> kk = k.reverse_complement();
                kk.nextReverseFromChar(::bliss::common::DNA::to_complement(r));
                return (k == kk);
            }
            // check to see if k+1-mer and its rc are identical.  append c to k on the left.
            static inline bool is_k1_rc_palindrome(unsigned char const & l, KmerType const & k) {
                if ((KMER_SIZE & 0x1) == 0) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (k.getCharsAtPos(0, 1) != ::bliss::common::DNA::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                ::bliss::common::Kmer<KMER_SIZE, ::bliss::common::DNA, WORD_TYPE> kk = k.reverse_complement();
                kk.nextFromChar(::bliss::common::DNA::to_complement(l));
                return (k == kk);
            }
            // check to see if k+2-mer and its rc are identical.  append c and k to both side.
            static inline bool is_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
                if ((KMER_SIZE & 0x1) > 0) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (r != ::bliss::common::DNA::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                return (k == k.reverse_complement());
            }
            // check to see if k+2-mer contains a k+1-mer palindrome.
            static inline unsigned char has_k1_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
            //     if ((KMER_SIZE & 0x1) == 0) return 0;  // even DNA.  so k-1-mer cannot be palindromic
            //     unsigned char res = 0;
            //     if (k.getCharsAtPos(KMER_SIZE - 1, 1) != ::bliss::common::DNA::to_complement(r)) res |= 1;   // two sides are not the same so not a rc_palindrome.
            //     if (k.getCharsAtPos(0, 1) != ::bliss::common::DNA::to_complement(l)) res |= 2;   // two sides are not the same so not a rc_palindrome.
            //     if ()
            // 
            //  TODO: this could be done with 1 revcomp check.
                return (is_k1_rc_palindrome(l, k) ? 2 : 0) | (is_k1_rc_palindrome(k, r) ? 1 : 0);
            }
            // check to see if k+2-mer and its rc are identical.  append c and k to both side.
            static inline unsigned char is_k_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
                if ((KMER_SIZE & 0x1) > 0) return false;  // even DNA.  so k-1-mer cannot be palindromic
                unsigned char res = ((k == k.reverse_complement()) ? 3 : 0);
                if (r != ::bliss::common::DNA::to_complement(l)) res &= 1;   // two sides are not the same so not a rc_palindrome.
                return res;
            }
            static inline unsigned char has_k1_rc_palindrome(KmerType const & k, bliss::debruijn::biedge::compact_simple_biedge & e) {
                return has_k1_rc_palindrome(bliss::common::DNA::FROM_ASCII[bliss::common::DNA16::TO_ASCII[e.getData()[0] >> 4]], 
                    k, 
                    bliss::common::DNA::FROM_ASCII[bliss::common::DNA16::TO_ASCII[e.getData()[0] & 0xf]]);
            }
        };

        template <unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
        struct kmer_traits<::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> > {
            using KmerType = ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;

            // check to see if kmer and its rc are identical
            static inline bool is_rc_palindrome(KmerType const & k) {
                return (k == k.reverse_complement());
            }
            // check to see if k+1-mer and its rc are identical.  append c to k on the right
            static inline bool is_k1_rc_palindrome(KmerType const & k, unsigned char const & r) {
                if (k.getCharsAtPos(KMER_SIZE - 1, 1) != ALPHABET::to_complement(r)) return false;   // two sides are not the same so not a rc_palindrome.
                ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> kk = k.reverse_complement();
                kk.nextReverseFromChar(ALPHABET::to_complement(r));
                return (k == kk);
            }
            // check to see if k+1-mer and its rc are identical.  append c to k on the left
            static inline bool is_k1_rc_palindrome(unsigned char const & l, KmerType const & k) {
                if (k.getCharsAtPos(0, 1) != ALPHABET::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> kk = k.reverse_complement();
                kk.nextFromChar(ALPHABET::to_complement(l));
                return (k == kk);
            }
            // check to see if k+2-mer and its rc are identical.  append c and k to both side.
            static inline bool is_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
                if (r != ALPHABET::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                return (k == k.reverse_complement());
            }
            // check to see if k+2-mer contains a k+1-mer palindrome.
            static inline unsigned char has_k1_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
            //     if ((KMER_SIZE & 0x1) == 0) return 0;  // even DNA.  so k-1-mer cannot be palindromic
            //     unsigned char res = 0;
            //     if (k.getCharsAtPos(KMER_SIZE - 1, 1) != ::bliss::common::DNA::to_complement(r)) res |= 1;   // two sides are not the same so not a rc_palindrome.
            //     if (k.getCharsAtPos(0, 1) != ::bliss::common::DNA::to_complement(l)) res |= 2;   // two sides are not the same so not a rc_palindrome.
            //     if ()
            //  TODO: this could be done with 1 revcomp check.
                return (is_k1_rc_palindrome(l, k) ? 2 : 0) | (is_k1_rc_palindrome(k, r) ? 1 : 0);
            }
            static inline unsigned char is_k_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
                unsigned char res = ((k == k.reverse_complement()) ? 3 : 0);
                if (r != ::bliss::common::DNA::to_complement(l)) res &= 1;   // two sides are not the same so not a rc_palindrome.
                return res;
            }
            static inline unsigned char has_k1_rc_palindrome(KmerType const & k, bliss::debruijn::biedge::compact_simple_biedge & e) {
                return has_k1_rc_palindrome(ALPHABET::FROM_ASCII[bliss::common::DNA16::TO_ASCII[e.getData()[0] >> 4]], 
                    k, 
                    ALPHABET::FROM_ASCII[bliss::common::DNA16::TO_ASCII[e.getData()[0] & 0xf]]);
            }
        };

    } /*namespace kmer*/
  }/*namespace common*/
}/*namespace bliss*/




#endif /* KMER_TRAITS_HPP_ */