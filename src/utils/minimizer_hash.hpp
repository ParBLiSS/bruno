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

/**
 * @file    minimizer_hash.hpp
 * @ingroup bliss::hash
 * @author  tpan
 * @brief   collections of hash functions defined for kmers, specifically to compute minimizers for use as hash values
 * @details support the following:  raw bits directly extracted; std::hash version; murmurhash; and farm hash
 *
 *          assuming the use is in a distributed hash table with total buckets N,  N = p * t * l,
 *          where p is number of processes, t is number of threads, and l is number of local buckets.
 *
 *          then the key is assigned to a bucket via hash(key) % N.
 *            the process assignment is: hash(key) / (t * l)
 *            the thread assignment is:  (hash(key) / l) % t
 *            the local bucket assignment is: hash(key) % l;
 *
 *          there are unlikely to be more than 2^64 local buckets, so we can limit the hash(key) % l to be the lower 64bit of hash(key).
 *          this also means that if the hash key is 64 bits, then no shifting or bit masking is needed, which improves performance for local hashtable lookup.
 *
 *          l is a variable that is decided by the local hash table based on number of entries.
 *          we should instead look at the first ceil(log (p*t)) bits of the hash(key).  let's call this "prefix".
 *
 *          process assignment is then:  prefix / t
 *          thread assignment is then:   prefix % t.
 *
 *          prefix for process assignment can be changed to use ceil(log(p)) bits, let's call this "pre-prefix".
 *
 *          process assignment is then:  pre-prefix % p.
 *
 *          2 functions are sufficient then:  prefix_hash(), and suffix_hash().
 *            we restrict our hash() functions to return 64 bits, and make suffix hash() to be the same as hash().
 *
 *
 *          namespace bliss::hash::kmer has a generic hash function that can work with kmer, kmer xor rev comp, computed or provided.
 *          the generic hash function also allows customization via bliss::hash::kmer::detail::{std,identity,murmur,farm} specializations
 *
 *          as stated above, 2 versions for each specialization: hash() and hash_prefix().  the specialization is especially for identity and murmur hashes,
 *          as murmur hash produces 128 bit value, and identity hash uses the original kmer.
 *
 */
#ifndef MMIZER_HASH_HPP_
#define MMIZER_HASH_HPP_

#ifdef __SSE4_1__
#include <x86intrin.h>
#endif

#include <tuple>  // for hash - std::pair
#include <exception>  // for hash - std::system_error
#include <algorithm>
#include <type_traits>  // enable_if

#include "common/alphabets.hpp"
#include "common/kmer.hpp"

#include "utils/transform_utils.hpp"

// includ the murmurhash code.
#ifndef _MURMURHASH3_H_
#include <smhasher/MurmurHash3.cpp>
#endif

// and farm hash
#ifndef FARM_HASH_H_
#include <farmhash/src/farmhash.cc>
#endif





namespace bliss {

  namespace kmer {

    namespace hash {
    

    
      namespace minimizer {
        // currently supports only DNA and DNA16-mers (not DNA5), and limits to 16 bit output (64k processors) due to intrinsics for horizontal min.
        // possible to contruct 32 bit output - the min 16-bit tells the starting position in 32 bit - more work is needed, though
        // can alternatively use logarithmic algo for hmin.
        // TODO:  1. Need to work with reverse complement.  - don't - that'd make it specific for k-mers.
        //        2. set up as a transform.
        //        3. optimize the code.
        class kmer_minimizer {
          protected:
            mutable uint8_t temp[16];
            const __m128i mask, ones;
          public:

            kmer_minimizer() : 
              mask(_mm_setr_epi16(0x0, 0x0, 0x0, 0xFFFF, 0x0, 0x0, 0x0, 0xFFFF)),
			        ones(_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)) {}


            // leverages the 2 lanes for offsets, each lane loads a 1 byte offset pair - this cuts the number of shift round to half - 1 .
            // the last 2 bytes in each 64bit processing needs to be reprocessed because of the shift op introduces 0.
            // each iteration processes 27 (DNA) or 15 (DNA16) charaters.  left lane and right lane are processed separately, until the end.
            //  due to the shifts and offset copies, good for up to
            //  DNA 31-mers (low lane: first 27-mers. high lane, 27-mers shifted by 4, which covers the 1/2 short shift.)
            //  DNA16 15-mers (low lane: first 13-mers. high lane, 13-mers shifted by 2, which covers the 1/2 short shift.)
            template <typename Kmer,
              typename std::enable_if<((Kmer::nBits + Kmer::bitsPerChar) <= 64), int>::type = 1>
            uint16_t operator()(Kmer const & kmer) const {
              using WordType = typename Kmer::KmerWordType;

              // first make a copy. then set the padding bits to 1
              Kmer k0 = kmer;
              k0.getDataRef()[Kmer::nWords - 1] |= getMostSignificantBitsMask<WordType>(::bliss::common::UnpaddedStreamTraits<WordType, Kmer::nBits>::padBits);
              unsigned char const * data = reinterpret_cast<unsigned char const *>(k0.getData());
              
              // next start iterating.  3 uint16_t each iteration, 2 offsets at a time, so shifting by 6 bytes at a time.
              __m128i curr, total;

//		std::cout << "total init: " << _mm_cvtsi128_si64(total) << std::endl;
              // 31-mer or smaller. the highest 2 bytes in each 64bit lane will not correspond to any usable bits.
              // we can iterate just once.
              // at exactly 32-mer, we'd need another iteration, and rely on the 0xFF of the rest of the register, so use the other version

              // load the data.
              memset(temp, 0xFF, 16);  // for the case where there is partial fill
              memcpy(&(temp[0]), data, std::min(8U, Kmer::nBytes));
              memcpy(&(temp[8]), data + 1, std::min(8U, Kmer::nBytes - 1U));
              // load the values
              curr = _mm_lddqu_si128(reinterpret_cast<__m128i *>(temp));

//		std::cout << "curr low: " << *(reinterpret_cast<size_t*>(temp)) << " hi: " << *(reinterpret_cast<size_t*>(temp + 8)) << std::endl;
              // we load into the 2 lanes 64 bits that are shifted by 1 byte because 
              // we have to shift in zeros, which makes the results of the last short inaccurate.
              //  but we can reprocess this byte later.  we avoid more complex shift.

              // accumulate the minimum for the shifts.
              total = curr;
              // DNA type
              switch (Kmer::bitsPerChar) {
                case 2: 
                  curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 2 chars
                  total = _mm_min_epu16(total, curr);
                  curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 3 chars
                  total = _mm_min_epu16(total, curr);
                case 4:
                  curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char
                  total = _mm_min_epu16(total, curr);
                default:
                  break;
              }
              // ignore the highest short in each lane. by setting it to max
              total = _mm_or_si128(total, mask);
              
              // hmin
              __m128i res = _mm_minpos_epu16(total);

//	    std::cout << "kmer: " << kmer << " copy " << k0 << " minimizer " << (_mm_cvtsi128_si32(res) & 0xFFFF) << " at pos " << (_mm_cvtsi128_si32(res) >> 16) << std::endl;
// 		printf("minimizer %u at pos %u\n", (_mm_cvtsi128_si32(res) & 0xFFFF), (_mm_cvtsi128_si32(res) >> 16));
              // lower 32 bit has the position and value (lowest 16 bit). should always be unsigned
              return static_cast<uint32_t>(_mm_cvtsi128_si32(res)) & static_cast<uint16_t>(0xFFFF);  
            }


            // good for 32-mer or larger:  loads 6 byte shifted into upper lane and do 7 shifts (fewer mem access)
            template <typename Kmer,
              typename std::enable_if<(Kmer::nBits >= 64), int>::type = 1>
            uint16_t operator()(Kmer const & kmer) const {
              using WordType = typename Kmer::KmerWordType;

              // first make a copy. then set the padding bits to 1
              Kmer k0 = kmer;
              k0.getDataRef()[Kmer::nWords - 1] |= getMostSignificantBitsMask<WordType>(::bliss::common::UnpaddedStreamTraits<WordType, Kmer::nBits>::padBits);
              unsigned char const * data = reinterpret_cast<unsigned char const *>(k0.getData());
              
              // next start iterating.  3 uint16_t each iteration, 2 offsets at a time, so shifting by 6 bytes at a time.
              __m128i curr, total;
              total = ones;

              // 32-mer or larger. the highest 2 bytes in each 64bit lane will not correspond to any usable bits.
              for (unsigned int i = 0; i < Kmer::nBytes; ) {
                // load the data.
                memset(temp, 0xFF, 16);  // for the case where there is partial fill
                memcpy(&(temp[0]), data + i, std::min(8U, Kmer::nBytes - i));
                i+=6;
                if (i < Kmer::nBytes) {
                  memcpy(&(temp[8]), data + i, std::min(8U, Kmer::nBytes - i));
                  i+=6;
                }
                // load the values
                curr = _mm_lddqu_si128(reinterpret_cast<__m128i *>(temp));

                // we load into the 2 lanes 64 bits that are shifted by 1 byte because 
                // we have to shift in zeros, which makes the results of the last short inaccurate.
                //  but we can reprocess this byte later.  we avoid more complex shift.

                // accumulate the minimum for the shifts.
                total = _mm_min_epu16(total, curr);  // unshifted.
                switch (Kmer::bitsPerChar) {
                  case 2:
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);

                  case 4:
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 3 times (12bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 3 times (12bits)
                    total = _mm_min_epu16(total, curr);
                  case 8:
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 1 times (8bits)
                    total = _mm_min_epu16(total, curr);
                    break;
                  default:
                    break;
                }
              }
              // ignore the last short in each lane. by setting it to max
              total = _mm_or_si128(total, mask);
              
              // hmin
              __m128i res = _mm_minpos_epu16(total);
              
              // lower 32 bit has the position and value (lowest 16 bit). should always be unsigned
              return static_cast<uint32_t>(_mm_cvtsi128_si32(res)) & static_cast<uint16_t>(0xFFFF);  
            }
        // DNA5 version is not implemented YET.

        };

        // currently supports only DNA and DNA16-mers (not DNA5), and limits to 16 bit output (64k processors) due to intrinsics for horizontal min.
        // possible to contruct 32 bit output - the min 16-bit tells the starting position in 32 bit - more work is needed, though
        // can alternatively use logarithmic algo for hmin.
        // TODO:  1. Need to work with reverse complement.  - don't - that'd make it specific for k-mers.
        //        2. set up as a transform.
        //        3. optimize the code.
        class kmolecule_minimizer {
          protected:
            mutable uint8_t temp[16];
            const __m128i mask, ones;
          public:

            kmolecule_minimizer() : 
              mask(_mm_setr_epi16(0x0, 0x0, 0x0, 0xFFFF, 0x0, 0x0, 0x0, 0xFFFF)),
			        ones(_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)) {}


            // leverages the 2 lanes for offsets, each lane loads a 1 byte offset pair - this cuts the number of shift round to half - 1 .
            // the last 2 bytes in each 64bit processing needs to be reprocessed because of the shift op introduces 0.
            // each iteration processes 27 (DNA) or 15 (DNA16) charaters.  left lane and right lane are processed separately, until the end.
            //  due to the shifts and offset copies, good for up to
            //  DNA 31-mers (low lane: first 27-mers. high lane, 27-mers shifted by 4, which covers the 1/2 short shift.)
            //  DNA16 15-mers (low lane: first 13-mers. high lane, 13-mers shifted by 2, which covers the 1/2 short shift.)


            // good for 32-mer or larger:  loads 6 byte shifted into upper lane and do 7 shifts (fewer mem access)
            template <typename Kmer>
            uint16_t operator()(Kmer const & kmer) const {
              using WordType = typename Kmer::KmerWordType;

              // first make a copy. then set the padding bits to 1
              Kmer k0 = kmer;
              Kmer krc;
              kmer.reverse_complement(krc);

              k0.getDataRef()[Kmer::nWords - 1] |= getMostSignificantBitsMask<WordType>(::bliss::common::UnpaddedStreamTraits<WordType, Kmer::nBits>::padBits);
              unsigned char const * data = reinterpret_cast<unsigned char const *>(k0.getData());
              krc.getDataRef()[Kmer::nWords - 1] |= getMostSignificantBitsMask<WordType>(::bliss::common::UnpaddedStreamTraits<WordType, Kmer::nBits>::padBits);
              unsigned char const * data_rc = reinterpret_cast<unsigned char const *>(krc.getData());
              
              // next start iterating.  3 uint16_t each iteration, 2 offsets at a time, so shifting by 6 bytes at a time.
              __m128i curr, total;
              total = ones;

              // 32-mer or larger. the highest 2 bytes in each 64bit lane will not correspond to any usable bits.
              for (unsigned int i = 0; i < Kmer::nBytes; i += 6) {
                // load the data.
                memset(temp, 0xFF, 16);  // for the case where there is partial fill
                memcpy(&(temp[0]), data + i, std::min(8U, Kmer::nBytes - i));
                memcpy(&(temp[8]), data_rc + i, std::min(8U, Kmer::nBytes - i));

                // load the values
                curr = _mm_lddqu_si128(reinterpret_cast<__m128i *>(temp));

                // we load into the 2 lanes 64 bits that are shifted by 1 byte because 
                // we have to shift in zeros, which makes the results of the last short inaccurate.
                //  but we can reprocess this byte later.  we avoid more complex shift.

                // accumulate the minimum for the shifts.
                total = _mm_min_epu16(total, curr);  // unshifted.
                switch (Kmer::bitsPerChar) {
                  case 2:
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 7 times (14bits)
                    total = _mm_min_epu16(total, curr);

                  case 4:
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 3 times (12bits)
                    total = _mm_min_epu16(total, curr);
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 3 times (12bits)
                    total = _mm_min_epu16(total, curr);
                  case 8:
                    curr = _mm_srli_epi64(curr, Kmer::bitsPerChar);  // shift right 1 char, 1 times (8bits)
                    total = _mm_min_epu16(total, curr);
                    break;
                  default:
                    break;
                }
              }
              // ignore the last short in each lane. by setting it to max
              total = _mm_or_si128(total, mask);
              
              // hmin
              __m128i res = _mm_minpos_epu16(total);
              
              // lower 32 bit has the position and value (lowest 16 bit). should always be unsigned
              return static_cast<uint32_t>(_mm_cvtsi128_si32(res)) & static_cast<uint16_t>(0xFFFF);  
            }
        // DNA5 version is not implemented YET.

        };


      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type.  max is 64bits.
       */
      template<typename KMER, bool Prefix = false, typename MMIZER=::bliss::kmer::hash::minimizer::kmolecule_minimizer>
      class cpp_std {
        protected:

          ::std::hash<uint16_t> op;
          MMIZER mini;

        public:
          static constexpr uint8_t batch_size = 1;

          // if nBits is more than 32, then default prefix should be 32.
          static constexpr unsigned int default_init_value = 16U;

          cpp_std(const unsigned int prefix_bits = default_init_value)  {};

          /// operator to compute hash
          inline size_t operator()(const KMER & kmer) const {
            if (Prefix) {
              throw ::std::logic_error("ERROR: minimizer with cpp_std hash currently does not support prefix.");
            } else {
              return op(mini(kmer));
            }
          }

      };
      template<typename KMER, bool Prefix, typename MMIZER>
      constexpr uint8_t cpp_std<KMER, Prefix, MMIZER>::batch_size;

      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER, bool Prefix = false, typename MMIZER=::bliss::kmer::hash::minimizer::kmolecule_minimizer>
      class identity {

        protected:
          unsigned int bits;
          MMIZER mini;

        public:
          static constexpr uint8_t batch_size = 1;

          static constexpr unsigned int default_init_value = 16U;
          static constexpr unsigned int suffix_bits = 48U;

          /// constructor
          identity(const unsigned int prefix_bits = default_init_value) : bits(std::min(KMER::nBits, prefix_bits)) {
          };

          /// operator to compute hash value
          inline uint64_t operator()(const KMER & kmer) const {
            if (Prefix) {
              throw ::std::logic_error("ERROR: minimizer with identity hash currently does not support prefix.");
            } else {
              return mini(kmer);
            }
            // to few bits to try to split it up.
          }
      };
      template<typename KMER, bool Prefix, typename MMIZER>
      constexpr uint8_t identity<KMER, Prefix, MMIZER>::batch_size;

      /**
       * @brief Kmer specialization for MurmurHash.  generated hash is 128 bit.
       *
       * TODO: move KMER type template param to operator.
       * TODO: change h to member variable.
       */
      template <typename KMER, bool Prefix = false, typename MMIZER=::bliss::kmer::hash::minimizer::kmolecule_minimizer>
      class murmur {


        protected:
          static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;
          uint32_t seed;
          MMIZER mini;

        public:
          static constexpr uint8_t batch_size = 1;

          static const unsigned int default_init_value = 24U;  // allow 16M processors.  but it's ignored here.

          murmur(const unsigned int prefix_bits = default_init_value, uint32_t const & _seed = 42 ) : seed(_seed) {};

          inline uint64_t operator()(const KMER & kmer) const
          {
            uint16_t minimizer = mini(kmer);


            // produces 128 bit hash.
            uint64_t h[2];
            // let compiler optimize out all except one of these.
            if (sizeof(void*) == 8)
              MurmurHash3_x64_128(reinterpret_cast<char *>(&minimizer), 2, seed, h);
            else if (sizeof(void*) == 4)
              MurmurHash3_x86_128(reinterpret_cast<char *>(&minimizer), 2, seed, h);
            else
              throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

            // use the upper 64 bits.
            if (Prefix)
              return h[1];
            else
              return h[0];
          }

      };
      template<typename KMER, bool Prefix, typename MMIZER>
      constexpr uint8_t murmur<KMER, Prefix, MMIZER>::batch_size;

      /**
       * @brief  Kmer hash, returns the least significant NumBits from murmur hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       *		IMPORTANT: upper and lower halves of the hash values are not independent, when farm hash is used with a given seed.
       *		therefore, different seed must be applied.  to adhere to existing api, we generate a different seed for "prefix" version.
       * @tparam Prefix:
       */
      template <typename KMER, bool Prefix = false, typename MMIZER=::bliss::kmer::hash::minimizer::kmolecule_minimizer>
      class farm {

        protected:
          static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;
          size_t shift;
          uint32_t seed;
          MMIZER mini;

        public:
          static constexpr uint8_t batch_size = 1;

          static const unsigned int default_init_value = 24U;   // this allows 16M processors.

          farm(const unsigned int prefix_bits = default_init_value, uint32_t const & _seed = 42 ) : shift(64U - std::min(prefix_bits, 64U)), seed(_seed)  {
          };

          /// operator to compute hash.  64 bit again.
          inline uint64_t operator()(const KMER & kmer) const {
            uint16_t minimizer = mini(kmer);

            if (Prefix)
              return ::util::Hash64WithSeed(reinterpret_cast<char*>(&minimizer), 2, (seed << 1) - 1);
            else
              return ::util::Hash64WithSeed(reinterpret_cast<char*>(&minimizer), 2, seed);
          }

      };
      template<typename KMER, bool Prefix, typename MMIZER>
      constexpr uint8_t farm<KMER, Prefix, MMIZER>::batch_size;

      } // namespace minimizer
    } // namespace hash
  } // namespace kmer
} // namespace bliss



#endif /* MMIZER_HASH_HPP_ */
