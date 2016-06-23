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
 * @file    mxx_support.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_IO_MXX_SUPPORT_HPP_
#define SRC_IO_MXX_SUPPORT_HPP_

#include <mxx/datatypes.hpp>

#include "common/kmer.hpp"

//#include "utils/system_utils.hpp"

#include <algorithm>  // for std::min

#include <unistd.h>   // for gethostname

#include "partition/range.hpp"
#include "common/sequence.hpp"
#include "debruijn/de_bruijn_node_trait.hpp"

namespace mxx {

  template<unsigned int size, typename A, typename WT>
    struct datatype_builder<typename bliss::common::Kmer<size, A, WT> > : 
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType , 
           bliss::common::Kmer<size, A, WT>::nWords> {

      typedef datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<unsigned int size, typename A, typename WT>
    struct datatype_builder<const typename bliss::common::Kmer<size, A, WT> > : 
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType , 
           bliss::common::Kmer<size, A, WT>::nWords> {

      typedef datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };

  
  template<>
    struct datatype_builder<bliss::common::LongSequenceKmerId> : 
    public datatype_builder<decltype(bliss::common::LongSequenceKmerId::id)> {

      typedef datatype_builder<decltype(bliss::common::LongSequenceKmerId::id)> baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<>
    struct datatype_builder<const bliss::common::LongSequenceKmerId> : 
    public datatype_builder<decltype(bliss::common::LongSequenceKmerId::id)> {

      typedef datatype_builder<decltype(bliss::common::LongSequenceKmerId::id)> baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<>
    struct datatype_builder<bliss::common::ShortSequenceKmerId> : 
    public datatype_builder<decltype(bliss::common::ShortSequenceKmerId::id)> {

      typedef datatype_builder<decltype(bliss::common::ShortSequenceKmerId::id)> baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<>
    struct datatype_builder<const bliss::common::ShortSequenceKmerId> : 
    public datatype_builder<decltype(bliss::common::ShortSequenceKmerId::id)> {

      typedef datatype_builder<decltype(bliss::common::ShortSequenceKmerId::id)> baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename T>
    struct datatype_builder<bliss::partition::range<T> > : 
    public datatype_contiguous<T , 
           sizeof(bliss::partition::range<T>)/sizeof(T) > {

      typedef datatype_contiguous<T,
      sizeof(bliss::partition::range<T>)/sizeof(T) > baseType;

      static MPI_Datatype get_type(){ 
        //== loop over the reads
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename T>
    struct datatype_builder<const bliss::partition::range<T> > : 
    public datatype_contiguous<T , 
           sizeof(bliss::partition::range<T>)/sizeof(T) > {

      typedef datatype_contiguous<T,
      sizeof(bliss::partition::range<T>)/sizeof(T) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename A, typename T>
    struct datatype_builder<::bliss::debruijn::graph::compact_edge<A, T> > :
    public datatype_builder<decltype(::bliss::debruijn::graph::compact_edge<A, T>::counts) > {

      typedef datatype_builder<decltype(::bliss::debruijn::graph::compact_edge<A, T>::counts) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename A, typename T>
    struct datatype_builder<const ::bliss::debruijn::graph::compact_edge<A, T> > :
    public datatype_builder<decltype(::bliss::debruijn::graph::compact_edge<A, T>::counts) > {

      typedef datatype_builder<decltype(::bliss::debruijn::graph::compact_edge<A, T>::counts) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };

//
//  template<typename A>
//    struct datatype_builder<bliss::de_bruijn::node::edge_exists<A> > :
//    public datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > {
//
//      typedef datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > baseType;
//
//      static MPI_Datatype get_type(){
//        return baseType::get_type();
//      }
//
//      static size_t num_basic_elements() {
//        return baseType::num_basic_elements();
//      }
//    };
//
//
//  template<typename A>
//    struct datatype_builder<const bliss::de_bruijn::node::edge_exists<A> > :
//    public datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > {
//
//      typedef datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > baseType;
//
//      static MPI_Datatype get_type(){
//        return baseType::get_type();
//      }
//
//      static size_t num_basic_elements() {
//        return baseType::num_basic_elements();
//      }
//    };
}

//std::ostream &operator<<(std::ostream &os, uint8_t const &t) {
//  return os << static_cast<uint32_t>(t);
//}
//std::ostream &operator<<(std::ostream &os, int8_t const &t) {
//  return os << static_cast<int32_t>(t);
//}
//
//template <typename T1, typename T2>
//std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> const &t) {
//  return os << t.first << ":" << t.second;
//}

#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
