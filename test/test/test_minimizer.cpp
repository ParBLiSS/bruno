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

// include google test
#include <gtest/gtest.h>
//#include <boost/concept_check.hpp>

#include <tuple>  // pair

// include classes to test
#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "utils/kmer_utils.hpp"
#include "utils/logging.h"
#include "utils/minimizer_hash.hpp"



// include files to test

//TESTS: Hash functions.  test boundary cases - the number of unique outputs should be relatively close to number of unique inputs.

template <typename T>
class KmerMinimizerTest : public ::testing::Test {
  protected:

    static std::vector<T> kmers;
    static constexpr size_t iterations = 100;

    using MINI_MER = bliss::common::Kmer<16/(T::bitsPerChar), typename T::KmerAlphabet, uint16_t>;

  public:
    static void SetUpTestCase()
    {
      T kmer;

      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {
          kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }

      kmers.resize(iterations);
      for (size_t i = 0; i < iterations; ++i) {
    	  kmers[i] = kmer;
          kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }
    }

    static void TearDownTestCase() {
    	std::vector<T>().swap(kmers);
    }

  protected:

    std::tuple<bool, int64_t, int64_t> test_minimizer(T const & kmer, MINI_MER const & mini_mer) {
        std::string kmstr(bliss::utils::KmerUtils::toASCIIString(kmer));
        std::string mmstr(bliss::utils::KmerUtils::toASCIIString(mini_mer));

        size_t kmlen = kmstr.length();
        size_t mmlen = mmstr.length();
        //std::cout << kmstr << " " << mmstr << std::endl;

        // compare.
        bool result = true;
        int comp = 0;
        int64_t pos = -1;
        int64_t lower_pos = -1;
        for (size_t i = 0; i <= (kmlen - mmlen); ++i) {  // kmlen - mmlen + 1
            comp = kmstr.compare(i, mmlen, mmstr);   // all substrings should be larger or equal.
            result = (comp >= 0);

            //std::cout << comp << " " ;
            if (comp == 0) {
                pos = i;
            } else if (comp < 0) {
                lower_pos = i;
            }
        }
        //std::cout << std::endl;

        return std::make_tuple(result, pos, lower_pos);
    }

    void minimize_kmers() {

       ::bliss::kmer::hash::minimizer::kmer_minimizer op;

      bool pass = true;
      T revcomp;
      MINI_MER mini_mer;
      int64_t match_pos = -1;
      int64_t lower_pos = -1;

      for (size_t i = 0; i < this->iterations; ++i) {
        mini_mer.getDataRef()[0] = op(this->kmers[i]);

        std::tie(pass, match_pos, lower_pos) = test_minimizer(this->kmers[i], mini_mer);

        if (!pass || (match_pos < 0) || (lower_pos >= 0)) {
            std::cout << "FAILED: kmer: i " << i << ", kmer: " << this->kmers[i] << ", minimizer: " << mini_mer << " at " << match_pos << " lower at " << lower_pos << std::endl; 
        }
        ASSERT_TRUE(pass);
      }

    }

    void minimize_kmols() {

       ::bliss::kmer::hash::minimizer::kmolecule_minimizer op;

      bool pass = true;
      int64_t match_pos = -1;
      int64_t lower_pos = -1;
      T revcomp;
      MINI_MER mini_mer;
      bool pass2 = true;
      int64_t match_pos2 = -1;
      int64_t lower_pos2 = -1;

      for (size_t i = 0; i < this->iterations; ++i) {
        mini_mer.getDataRef()[0] = op(this->kmers[i]);
        this->kmers[i].reverse_complement(revcomp);

        std::tie(pass, match_pos, lower_pos) = test_minimizer(this->kmers[i], mini_mer);
        std::tie(pass2, match_pos2, lower_pos2) = test_minimizer(revcomp, mini_mer);

        pass &= pass2;
        if (!pass || ((match_pos < 0) && (match_pos2 < 0)) || (lower_pos >= 0) || (lower_pos2 >= 0)) 
            std::cout << "FAILED: kmol revcomp: i " << i << ", kmer: " << this->kmers[i] << ", revcomp: " << revcomp << ", minimizer: " << mini_mer << " at str1 " << match_pos << " str2 " << match_pos2 << " lower1 at " << lower_pos << " lower2 at " << lower_pos2 << std::endl; 

        ASSERT_TRUE(pass);
      }

    }
};

template <typename T>
constexpr size_t KmerMinimizerTest<T>::iterations;

template <typename T>
std::vector<T> KmerMinimizerTest<T>::kmers;



// indicate this is a typed test
TYPED_TEST_CASE_P(KmerMinimizerTest);

TYPED_TEST_P(KmerMinimizerTest, minimize)
{
	this->minimize_kmers();
}

TYPED_TEST_P(KmerMinimizerTest, bi_minimize)
{
	this->minimize_kmols();
}





REGISTER_TYPED_TEST_CASE_P(KmerMinimizerTest, minimize, bi_minimize);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint64_t>,  // kmer_minimizer has 2 impl separated by a threshold
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>, 
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>, 
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>, 
    ::bliss::common::Kmer< 9, bliss::common::DNA,   uint32_t>,   
    ::bliss::common::Kmer< 8, bliss::common::DNA,   uint16_t>,   // below minimizer threshold
    ::bliss::common::Kmer< 8, bliss::common::DNA,   uint32_t>,   // below minimizer threshold
    ::bliss::common::Kmer< 5, bliss::common::DNA16, uint32_t>,  
    ::bliss::common::Kmer< 4, bliss::common::DNA16, uint16_t>,  
    ::bliss::common::Kmer< 4, bliss::common::DNA16, uint32_t>,  
    ::bliss::common::Kmer< 63, bliss::common::DNA,   uint64_t>,   // large
    ::bliss::common::Kmer< 31, bliss::common::DNA16, uint64_t>,  
    ::bliss::common::Kmer< 127, bliss::common::DNA,   uint64_t>,   // very large
    ::bliss::common::Kmer< 63, bliss::common::DNA16, uint64_t>,  
    ::bliss::common::Kmer< 255, bliss::common::DNA,   uint64_t>,   // super large
    ::bliss::common::Kmer< 127, bliss::common::DNA16, uint64_t>  
> KmerMinimizerTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bruno, KmerMinimizerTest, KmerMinimizerTestTypes);
