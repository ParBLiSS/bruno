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
        
        template <unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
        struct kmer_traits<::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> > {
            using KmerType = ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;

	    // ============  BASE CHECKES

	    // k-1-mer palindrome checker.
	    static inline bool is_kminus1_rc_palindrome_high(KmerType const & k) {
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
	        ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> kk = k.reverse_complement_shift(1);  // reverse complement and left shift 1.
                return (k.masked_equal_MSK_1(kk));  // masked compare excluding 1 most significant character.
     	    }
	    static inline bool is_kminus1_rc_palindrome_low(KmerType const & k) {
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
	        ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> kk = k.reverse_complement_shift(-1);  // reverse complement and right shift 1.
                return (k.masked_equal_LSK_1(kk));  // masked compare excluding 1 most significant character.
     	    }

            // check to see if kmer and its rc are identical
            static inline bool is_rc_palindrome(KmerType const & k) {
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) > 0)) return false;  // odd DNA.  cannot be palindromic 
                return (k == k.reverse_complement());
            }

	    // ============ K1 and k1pair checks based on k-1 and k palindrome checks.

            // check to see if k+1-mer and its rc are identical.  
            static inline bool is_k1_rc_palindrome(KmerType const & k, unsigned char const & r, bool palindromic_kminus1_low) {
                assert((r < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (k.getCharsAtPos(KMER_SIZE - 1, 1) != ALPHABET::to_complement(r)) return false;   // two sides are not the same so not a rc_palindrome.
		return palindromic_kminus1_low;
            }
            // check to see if k+1-mer and its rc are identical. 
            static inline bool is_k1_rc_palindrome(unsigned char const & l, KmerType const & k, bool palindromic_kminus1_high) {
                assert((l < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (k.getCharsAtPos(0, 1) != ALPHABET::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
		return palindromic_kminus1_high;
            }
            // two k-mers form k1pair rc_palindrome when each, when extended with an edge in the same direction, forms the other kmer.  rc(x[2..k]c)=y, and rc(y[2..k]d) = x.  c=rc(y[1]), d=rc(x[1])
                //   so rc(x[2..k]) = y[2..k], and vice versa, so x[2..k] and y[2..k] are palindromes.
            static inline bool is_k1pair_rc_palindrome(KmerType const & k, unsigned char const & r, bool palindromic_kminus1_low) {
                assert((r < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (k.getCharsAtPos(KMER_SIZE - 1, 1) != r) return false;   // two sides are not the same so not a rc_palindrome.
		return palindromic_kminus1_low;
            }
            // two k-mers form k1pair rc_palindrome when each, when extended with an edge in the same direction, forms the other kmer.  rc(x[2..k]c)=y, and rc(y[2..k]d) = x.  c=rc(y[1]), d=rc(x[1])
                //   so rc(x[2..k]) = y[2..k], and vice versa, so x[2..k] and y[2..k] are palindromes.
            static inline bool is_k1pair_rc_palindrome(unsigned char const & l, KmerType const & k, bool palindromic_kminus1_high) {
                assert((l < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                if (k.getCharsAtPos(0, 1) != l) return false;   // two sides are not the same so not a rc_palindrome.
		return palindromic_kminus1_high;
            }

	    // ==============  combination check for k1 and k1pair palindrome

            static inline unsigned char is_k1_k1pair_rc_palindrome(KmerType const & k, unsigned char const & r, bool palindromic_kminus1_low) {
                assert((r < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return 0;  // even DNA.  so k-1-mer cannot be palindromic
                
                if (!palindromic_kminus1_low) return 0;  // masked compare excluding 1 most significant character.
                return 
                 ((k.getCharsAtPos(KMER_SIZE - 1, 1) == ALPHABET::to_complement(r)) ? 1 : 0) |   // rc palindrome.
                 ((k.getCharsAtPos(KMER_SIZE - 1, 1) == r) ? 2 : 0);   // paired palindrome.
            }
            static inline unsigned char is_k1_k1pair_rc_palindrome(unsigned char const & l, KmerType const & k, bool palindromic_kminus1_high) {
                assert((l < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return 0;  // even DNA.  so k-1-mer cannot be palindromic

                if (!palindromic_kminus1_high) return 0;  // masked compare excluding 1 most significant character.
                return 
                 ((k.getCharsAtPos(0, 1) == ALPHABET::to_complement(l)) ? 1 : 0) |   // rc palindrome.
                 ((k.getCharsAtPos(0, 1) == l) ? 2 : 0);   // paired palindrome.
            }

	    // =============== k2 checks

            // check to see if k+2-mer and its rc are identical.  
            static inline bool is_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r, bool palindromic_k) {
                assert((r < ALPHABET::SIZE) && (l < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) > 0)) return false;  // odd DNA.  so k2-mer cannot be palindromic
                if (r != ALPHABET::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                return palindromic_k;
            }

	    
	    // =============== conveience versions
 
            // check to see if k+1-mer and its rc are identical.  
            static inline bool is_k1_rc_palindrome(KmerType const & k, unsigned char const & r) {
                // if (std::is_saddme<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                // if (k.getCharsAtPos(KMER_SIZE - 1, 1) != ALPHABET::to_complement(r)) return false;   // two sides are not the same so not a rc_palindrome.
                // kk.nextReverseFromChar(ALPHABET::to_complement(r));
                // return (k == kk)
                return is_k1_rc_palindrome(k, r, is_kminus1_rc_palindrome_low(k));  // masked compare excluding 1 most significant character.
            }
            // check to see if k+1-mer and its rc are identical. 
            static inline bool is_k1_rc_palindrome(unsigned char const & l, KmerType const & k) {
                // if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                // if (k.getCharsAtPos(0, 1) != ALPHABET::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                // kk.nextFromChar(ALPHABET::to_complement(l));
                //return (k == kk);
                return is_k1_rc_palindrome(l, k, is_kminus1_rc_palindrome_high(k));  // masked compare excluding 1 most significant character.
            }
            // two k-mers form k1pair rc_palindrome when each, when extended with an edge in the same direction, forms the other kmer.  rc(x[2..k]c)=y, and rc(y[2..k]d) = x.  c=rc(y[1]), d=rc(x[1])
                //   so rc(x[2..k]) = y[2..k], and vice versa, so x[2..k] and y[2..k] are palindromes.
            static inline bool is_k1pair_rc_palindrome(KmerType const & k, unsigned char const & r) {
                // if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                // if (k.getCharsAtPos(KMER_SIZE - 1, 1) != r) return false;   // two sides are not the same so not a rc_palindrome.
                // kk.nextReverseFromChar(ALPHABET::to_complement(r));
                // return (k == kk);
                return is_k1pair_rc_palindrome(k, r, is_kminus1_rc_palindrome_low(k));  // masked compare excluding 1 most significant character.
            }
            // two k-mers form k1pair rc_palindrome when each, when extended with an edge in the same direction, forms the other kmer.  rc(x[2..k]c)=y, and rc(y[2..k]d) = x.  c=rc(y[1]), d=rc(x[1])
                //   so rc(x[2..k]) = y[2..k], and vice versa, so x[2..k] and y[2..k] are palindromes.
            static inline bool is_k1pair_rc_palindrome(unsigned char const & l, KmerType const & k) {
                // if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                // if (k.getCharsAtPos(0, 1) != l) return false;   // two sides are not the same so not a rc_palindrome.
                // kk.nextFromChar(ALPHABET::to_complement(l));
                // return (k == kk);
                return is_k1pair_rc_palindrome(l, k, is_kminus1_rc_palindrome_high(k));  // masked compare excluding 1 most significant character.
            }


            static inline unsigned char is_k1_k1pair_rc_palindrome(KmerType const & k, unsigned char const & r) {
                //if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return 0;  // even DNA.  so k-1-mer cannot be palindromic
                //
                //if (!is_kminus1_rc_palindrome_low(k)) return 0;  // masked compare excluding 1 most significant character.
                //return 
                // ((k.getCharsAtPos(KMER_SIZE - 1, 1) == ALPHABET::to_complement(r)) ? 1 : 0)) |   // rc palindrome.
                // ((k.getCharsAtPos(KMER_SIZE - 1, 1) == r) ? 2 : 0));   // paired palindrome.
		return is_k1_k1pair_rc_palindrome(k, r, is_kminus1_rc_palindrome_low(k));
            }
            static inline unsigned char is_k1_k1pair_rc_palindrome(unsigned char const & l, KmerType const & k) {
                //if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) == 0)) return 0;  // even DNA.  so k-1-mer cannot be palindromic
		//
                //if (!is_kminus1_rc_palindrome_high(k)); return 0;  // masked compare excluding 1 most significant character.
                //return 
                // ((k.getCharsAtPos(0, 1) == ALPHABET::to_complement(l)) ? 1 : 0)) |   // rc palindrome.
                // ((k.getCharsAtPos(0, 1) == l) ? 2 : 0));   // paired palindrome.
		return is_k1_k1pair_rc_palindrome(l, k, is_kminus1_rc_palindrome_high(k));
            }

            // check to see if k+2-mer and its rc are identical.  
            static inline bool is_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
                //if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) > 0)) return false;  // odd DNA.  so k2-mer cannot be palindromic
                //if (r != ALPHABET::to_complement(l)) return false;   // two sides are not the same so not a rc_palindrome.
                //return (k == k.reverse_complement());
		return is_k2_rc_palindrome(l, k, r, is_rc_palindrome(k));
            }

            // // check to see if k+2-mer contains a k+1-mer palindrome.  3: both are k1 palindrome. 2: in edge.  1. out edge.
            // static inline unsigned char has_k1_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
            // //  TODO: this could be done with 1 revcomp check.
            //     return (is_k1_rc_palindrome(l, k) ? 2 : 0)) | (is_k1_rc_palindrome(k, r) ? 1 : 0));
            // }
            // check to see if k-mer and k+2-mer are rc_palindromic.  3: kmer and k2mer.  1. kmer only.
            static inline unsigned char is_k_k2_rc_palindrome(unsigned char const & l, KmerType const & k, unsigned char const & r) {
                assert((r < ALPHABET::SIZE) && (l < ALPHABET::SIZE) && "character in wrong alphabet" );
                if (std::is_same<ALPHABET, ::bliss::common::DNA>::value && ((KMER_SIZE & 0x1) > 0)) return false;  // even DNA.  so k-1-mer cannot be palindromic
                unsigned char res = ((k == k.reverse_complement()) ? 3 : 0);
                if (r != ALPHABET::to_complement(l)) res &= 1;   // two sides are not the same so not a rc_palindrome.
                return res;
            }
            
            // static inline unsigned char has_k1_rc_palindrome(KmerType const & k, bliss::debruijn::biedge::compact_simple_biedge & e) {
            //     return has_k1_rc_palindrome(bliss::common::DNA::FROM_ASCII[bliss::common::DNA16::TO_ASCII[e.getData()[0] >> 4]], 
            //         k, 
            //         bliss::common::DNA::FROM_ASCII[bliss::common::DNA16::TO_ASCII[e.getData()[0] & 0xf]]);
            // }

        };
    } /*namespace kmer*/
  }/*namespace common*/
}/*namespace bliss*/




#endif /* KMER_TRAITS_HPP_ */
