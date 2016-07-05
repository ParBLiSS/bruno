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

/**
 * range_test.cpp
 * Test range class
 *  Created on: Feb 18, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */

#include "io/data_block.hpp"
#include "partition/range.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <vector>
#include <list>
#include <limits>
#include <algorithm>

using namespace bliss::io;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DataBlockTest : public ::testing::Test
{
  protected:
    std::vector<T> src;
    size_t slen;
    size_t dlen;

    virtual void SetUp()
    {
      slen = (std::numeric_limits<T>::max() < 1000000) ? std::numeric_limits<T>::max() : 1000000;

      src.reserve(slen);
      for (size_t i = 0; i < slen; ++i) {
        src.push_back(i);
      }

      dlen = (slen / 2 < 1000) ? slen / 2 : 1000;
    }

    virtual void TearDown() {
    	std::vector<T>().swap(src);
    }

};

// indicate this is a typed test
TYPED_TEST_CASE_P(DataBlockTest);


// testing the equal function
TYPED_TEST_P(DataBlockTest, NoBuffer){
  typedef bliss::io::UnbufferedDataBlock<typename std::vector<TypeParam>::iterator, bliss::partition::range<size_t> >    DB_Type;

  bliss::partition::range<size_t> r(0, this->slen);
  DB_Type db;
  db.assign(this->src.begin(), this->src.end(), r);
  EXPECT_EQ(this->src.begin(), db.begin());
  EXPECT_EQ(this->src.end(),   db.end()  );


  bliss::partition::range<size_t> r2(10, this->slen / 2);
  DB_Type db2;
  db2.assign(this->src.begin() + r2.start, this->src.begin() + r2.end, r2);
  EXPECT_EQ(this->src.begin()+ r2.start, db2.begin());
  EXPECT_EQ(this->src.begin()+ r2.end  , db2.end()  );


}

// testing the equal function
TYPED_TEST_P(DataBlockTest, BufferList){
  typedef bliss::io::DataBlock<typename std::vector<TypeParam>::iterator, bliss::partition::range<size_t>, std::list<TypeParam> > DB_Type;

  bliss::partition::range<size_t> r(10, 10 + this->dlen);
  DB_Type db;
  db.assign(this->src.begin() + r.start, this->src.begin() + r.end, r);

  auto i = r.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

  std::vector<TypeParam> result1(db.begin(), db.end());



  bliss::partition::range<size_t> r2(15, 15 + this->dlen);
  db.assign(this->src.begin() + r2.start, this->src.begin() + r2.end, r2);


  i = r2.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Shifted Vectors x and y differ at index " << i;
  }

  i = r2.start;
  auto it2 =  result1.begin();
  for (auto it = db.begin(); it != db.end() && i < r2.end; ++it, ++i, ++it2) {
    EXPECT_NE(*it2, *it) << "Vectors x and y same at index " << i;
  }

}


// testing the equal function
TYPED_TEST_P(DataBlockTest, BufferVector){
  typedef bliss::io::BufferedDataBlock<typename std::vector<TypeParam>::iterator, bliss::partition::range<size_t> > DB_Type;

  bliss::partition::range<size_t> r(10, 10 + this->dlen);
  DB_Type db;
  db.assign(this->src.begin() + r.start, this->src.begin() + r.end, r);

  auto i = r.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

  std::vector<TypeParam> result1(db.begin(), db.end());



  bliss::partition::range<size_t> r2(15, 15 + this->dlen);
  db.assign(this->src.begin() + r2.start, this->src.begin() + r2.end, r2);


  i = r2.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Shifted Vectors x and y differ at index " << i;
  }

  i = r2.start;
  auto it2 = result1.begin();
  for (auto it = db.begin(); it != db.end() && i < r2.end; ++it, ++i, ++it2) {
    EXPECT_NE(*it2, *it) << "Vectors x and y same at index " << i;
  }

}





// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DataBlockTest, NoBuffer, BufferList, BufferVector);


/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DataBlockPtrTest : public ::testing::Test
{
  protected:
    T* src;
    size_t slen;
    size_t dlen;

    virtual void SetUp()
    {
      slen = (std::numeric_limits<T>::max() < 1000000) ? std::numeric_limits<T>::max() : 1000000;

      src = new T[slen];
      for (size_t i = 0; i < slen; ++i) {
        src[i] = i;
      }

      dlen = (slen / 2 < 1000) ? slen / 2 : 1000;
    }

    virtual void TearDown()
    {
      delete [] src;
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DataBlockPtrTest);


// testing the equal function
TYPED_TEST_P(DataBlockPtrTest, NoBuffer){
  typedef bliss::io::UnbufferedDataBlock<TypeParam*, bliss::partition::range<size_t> > DB_Type;

  bliss::partition::range<size_t> r(0, this->slen);
  DB_Type db;
  db.assign(this->src, this->src + this->slen, r);
  EXPECT_EQ(this->src,                db.begin());
  EXPECT_EQ((this->src + this->slen),   db.end()  );

  int i = r.start;
  for (typename DB_Type::iterator it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

  bliss::partition::range<size_t> r2(10, this->slen / 2);
  DB_Type db2;
  db2.assign(this->src + r2.start, this->src + r2.end, r2);
  EXPECT_EQ(this->src + r2.start, db2.begin());
  EXPECT_EQ(this->src + r2.end  , db2.end()  );

  i = r2.start;
  for (typename DB_Type::iterator it = db2.begin(); it != db2.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

}



// testing the equal function
TYPED_TEST_P(DataBlockPtrTest, BufferVector){
  typedef bliss::io::BufferedDataBlock<TypeParam*, bliss::partition::range<size_t> > DB_Type;

  bliss::partition::range<size_t> r(10, 10 + this->dlen);
  DB_Type db;
  db.assign(this->src + r.start, this->src + r.end, r);

  auto i = r.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }


  std::vector<TypeParam> result1(db.begin(), db.end());

  bliss::partition::range<size_t> r2(15, 15 + this->dlen);
  db.assign(this->src + r2.start, this->src + r2.end, r2);

  i = r2.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "shifted Vectors x and y differ at index " << i;
  }


  i = r2.start;
//  BL_INFO( "length = " << this->dlen );
  auto it2 = result1.begin();
  for (auto it = db.begin(); it != db.end() && i < r2.end; ++it, ++i, ++it2) {
//    BL_INFO( (int)(*it2) << " vs " << (int)(this->src[i]) );
    EXPECT_NE(*it2, *it) << "Vectors x and y same at index " << i;
  }


}

// testing the equal function
TYPED_TEST_P(DataBlockPtrTest, BufferList){
  typedef bliss::io::BufferedDataBlock<TypeParam*, bliss::partition::range<size_t> > DB_Type;

  bliss::partition::range<size_t> r(10, 10 + this->dlen);
  DB_Type db;
  db.assign(this->src + r.start, this->src + r.end, r);

  auto i = r.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }


  std::vector<TypeParam> result1(db.begin(), db.end());

  bliss::partition::range<size_t> r2(15, 15 + this->dlen);
  db.assign(this->src + r2.start, this->src + r2.end, r2);

  i = r2.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "shifted Vectors x and y differ at index " << i;
  }


  i = r2.start;
//  BL_INFO( "length = " << this->dlen );
  auto it2 = result1.begin();
  for (auto it = db.begin(); it != db.end() && i < r2.end; ++it, ++i, ++it2) {
//    BL_INFO( (int)(*it2) << " vs " << (int)(this->src[i]) );
    EXPECT_NE(*it2, *it) << "Vectors x and y same at index " << i;
  }


}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DataBlockPtrTest, NoBuffer, BufferVector, BufferList);



//////////////////////////
////  DEATH TESTS - test class named BlahDeathTest so gtest will not run these in a threaded context. and will run first
//
//// typedef RangeTest DataBlockDeathTest
//template<typename T>
//class DataBlockDeathTest : public ::testing::Test
//{
//  protected:
//    std::vector<T> src;
//    std::vector<T> dest;
//    size_t slen;
//    size_t dlen;
//
//    virtual void SetUp()
//    {
//      slen = (std::numeric_limits<T>::max() < 1000000) ? std::numeric_limits<T>::max() : 1000000;
//
//      src.reserve(slen);
//      for (T i = 0; i < slen; ++i) {
//        src.emplace_back(i);
//      }
//
//      dlen = (std::numeric_limits<T>::max() < 1000) ? std::numeric_limits<T>::max() : 1000;
//      dlen /= 2;
//      dest.reserve(dlen);
//    }
//};
//

//// annotate that DataBlockDeathTest is typed
//TYPED_TEST_CASE_P(DataBlockDeathTest);
//
//// failed construction due to asserts
//TYPED_TEST_P(DataBlockDeathTest, badRange){
//
//  typedef bliss::io::DataBlock<typename std::vector<TypeParam>::iterator, bliss::partition::range<size_t>> DB_Type;
//
//  bliss::partition::range<size_t> r2(0, this->slen);
//  std::string err_regex = ".data_block.hpp.* Assertion .* failed.*";
//
//  // basically, if start is larger than end.
//  EXPECT_EXIT(DB_Type(this->src.begin(), this->src.begin() + this->dlen, r2, this->dest.begin(), this->dlen), ::testing::KilledBySignal(SIGABRT), err_regex);
//}
//
//
//
//// register the death test cases
//REGISTER_TYPED_TEST_CASE_P(DataBlockDeathTest, badRange);


//////////////////// RUN the tests with different types.

//typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
//    int64_t, uint64_t, size_t> DataBlockTestTypes;
typedef ::testing::Types<int8_t, uint8_t> DataBlockTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DataBlockTest, DataBlockTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DataBlockPtrTest, DataBlockTestTypes);
//INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DataBlockDeathTest, DataBlockTestTypes);
