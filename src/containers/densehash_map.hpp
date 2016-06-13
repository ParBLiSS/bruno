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
 * @file    densehash_map.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_WIP_DENSEHASH_VECMAP_HPP_
#define SRC_WIP_DENSEHASH_VECMAP_HPP_

#include <sparsehash/dense_hash_map>
#include <vector>
#include <functional>  // hash, equal_to, etc
#include <tuple>   // pair
#include <scoped_allocator>
#include <algorithm>
#include <cmath>   // ceil
#include <memory>  // allocator
#include <iostream>

#include "iterators/concatenating_iterator.hpp"

#include "containers/fsc_container_utils.hpp"

#include "utils/logging.h"

namespace fsc {  // fast standard container

  /**
   * @brief my version of hashed map.  because std::unordered_map does chaining for collision using a LINKED LIST.
   * @details std::unordered_map's use of linked list is fine for low collision.  for high collision, it's great for insertion but terrible for search or delete
   *          it is also not suitable for copying to a vector and sort/search.  in addition, the search (equal_range, count, and probably erase) are implemented
   *          using linear search.
   *
   *          This class attempts to address this shortcoming by maintaining the unordered map interface while replacing the LINKED_LIST with std::vector or std::multimap
   *
   *          Also note: google dense hash has strong requirements for the hash object to be move constructable/assignable (in our case, TransformedHash, and the underlying hash functions.
   *            because of the special constructors requirements, copy constructor and default constructors also need to be defined.
   *          Performance:  with high repeat input set (e.g. test.fastq, with 254K repeats per kmer, 40 unique kmers) std::unordered_map search time is approximately 1.5 sec for 1%, count is 0.56 sec.  build is 1.6 s.  pos+qual map)
   *            hashmap:  build 1.6s, find 1.5s, count 0.56s
   *            sorted vector:  build = 2 s, find 0.27s, count 0.01s
   *          with google dense_hash_map:  can't use it - google dense_hash_map requires uniqueness, and also requiers a NULL element.
   *
   *          tested using std::multimap and map.  VERY slow on insertion, especially for real data.  about 33% faster for find and count compared to unordered map for synthetic, high repeat datasets.  still about 4 times slower than sort based.
   *          for that dataset.  build = 3.6s, find 0.98s, count 0.73s
   *
   *          WE HAVE TO BUILD OUR OWN VERSION.
   *            prefer using std::vector as internal - self growing,
   *
   *          Note that this class has an incomplete implementation of the multimap interface - only those interfaces that I needed are implemented here.
   *
   *          DESIGN:
   *          instead of requiring sorting, since unordered_map (without multimap) has reasonable performance, we can
   *          build a multimap using unordered map by setting the mapped_type to be ::std::vector<::std::pair<Key, T> >
   *
   *          instead of building this at the distributed map level, we build it at the local unordered map level,
   *          which is easier to debug and benchmark.
   *
   *      vector holds pair<Key, T> because then we can directly copy (fast)
   *      a version with pair<T> would be more space efficient but slower when large amount of data is present.
   *
   *          this class internally uses a hashmap of vectors.
   *
   *
   *      memory usage:
   *        supercontainer = map, stores ::std::pair<K, std::vector<> > in linked list for chaining.
   *          each "bucket" hash size_t + pointer to head of linked list. - 16 bytes, HU hash unique elements.
   *          each link list node has ::std::pair<K, std::vector> as payload, and at least a next ptr. - 24 bytes,  U unique elements
   *          each vector stores  std::pair<K, T>, so 16 or 24 bytes.  N elements
   *
   *          total: (16N or 24N) + 24U + 16 HU.  assume good hash function, HU = U.
   *
   *          bottomline - large amount of memory is needed.
   *
   *
   *
   *
   *  this file contains the GOOGLE DENSE HASH MAP version of map and multimap that are compatible with kmer indexing.
   */



// key values span entire key space.
template <typename Key,
typename T,
typename LowerKeySpaceSelector,
typename Hash = ::std::hash<Key>,
typename Equal = ::std::equal_to<Key>,
typename Allocator = ::std::allocator<::std::pair<const Key, T> > >
class densehash_map {

  protected:

    using container_type =
        ::google::dense_hash_map<Key, T,
                       Hash, Equal, Allocator >;

    container_type lower_map;
    container_type upper_map;

    using container_iterator = typename container_type::iterator;
    using container_const_iterator = typename container_type::const_iterator;
    using container_range = ::std::pair<container_iterator, container_iterator>;
    using container_const_range = ::std::pair<container_const_iterator, container_const_iterator>;


    LowerKeySpaceSelector splitter;


    template <typename InputIt>
    InputIt partition_input(InputIt first, InputIt last) {
    	return ::std::stable_partition(first, last, splitter);
    }

  public:
    using key_type              = Key;
    using mapped_type           = T;
    using value_type            = ::std::pair<const Key, T>;
    using hasher                = Hash;
    using key_equal             = Equal;
    using allocator_type        = Allocator;
    using reference             = value_type&;
    using const_reference       = const value_type&;
    using pointer               = typename std::allocator_traits<Allocator>::pointer;
    using const_pointer         = typename std::allocator_traits<Allocator>::const_pointer;
    using iterator              = ::bliss::iterator::ConcatenatingIterator<container_iterator >;
    using const_iterator        = ::bliss::iterator::ConcatenatingIterator<container_const_iterator >;
    using size_type             = size_t;
    using difference_type       = ptrdiff_t;

    densehash_map(size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         lower_map(bucket_count / 2, hash, equal, alloc),
						 upper_map(bucket_count / 2, hash, equal, alloc)
						 {};

    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
    densehash_map(Key const & empty_key, Key const & deleted_key,
    		 	 	   Key const & upper_empty_key, Key const & upper_deleted_key,
					   const LowerKeySpaceSelector & _splitter = LowerKeySpaceSelector(),
             size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         lower_map(bucket_count / 2, hash, equal, alloc),
						 upper_map(bucket_count / 2, hash, equal, alloc)
						 {
    	reserve_keys(empty_key, deleted_key);
    	reserve_upper_keys(upper_empty_key, upper_deleted_key, _splitter);

//    	lower_map.max_load_factor(0.8);
//    	lower_map.min_load_factor(0.0);
//      upper_map.max_load_factor(0.8);
//      upper_map.min_load_factor(0.0);
    };

    template<class InputIt>
    densehash_map(InputIt first, InputIt last,
    		Key const & empty_key, Key const & deleted_key,
    		    		 	 	   Key const & upper_empty_key, Key const & upper_deleted_key,
                       size_type bucket_count = 128,
					   const LowerKeySpaceSelector & _splitter = LowerKeySpaceSelector(),
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
					   densehash_map(empty_key, deleted_key, upper_empty_key, upper_deleted_key,
							   _splitter, bucket_count, hash, equal, alloc) {

    	this->insert(first, last);
    };

    virtual ~densehash_map() {};

    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
        lower_map.set_empty_key(empty_key);
        lower_map.set_deleted_key(deleted_key);
    }

    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key, const LowerKeySpaceSelector & _splitter = LowerKeySpaceSelector()) {
    	printf("reserve_upper_keys\n");
        upper_map.set_empty_key(empty_key);
        upper_map.set_deleted_key(deleted_key);
        splitter = _splitter;
    }



    iterator begin() {
    	return iterator(std::vector<container_range> { container_range{lower_map.begin(), lower_map.end()},
                                                    container_range{upper_map.begin(), upper_map.end()} });
    }
    const_iterator begin() const {
      return cbegin();
    }
    const_iterator cbegin() const {
      return const_iterator(std::vector<container_const_range> { container_const_range{lower_map.begin(), lower_map.end()},
                                                                 container_const_range{upper_map.begin(), upper_map.end()} });
    }



    iterator end() {
      return iterator( upper_map.end() );
    }
    const_iterator end() const {
      return cend();
    }
    const_iterator cend() const {
      return const_iterator( upper_map.end() );
    }


    std::vector<Key> keys() const  {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }

    void keys(std::vector<Key> & ks) const  {
      ks.clear();
      ks.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        ks.emplace_back(it->first);
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        ks.emplace_back(it->first);
      }

      return ks;
    }

    std::vector<std::pair<Key, T>> to_vector() const  {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }

    void to_vector(  std::vector<std::pair<Key, T>> & vs) const  {
      vs.clear();
      vs.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        vs.emplace_back(*it);
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        vs.emplace_back(*it);
      }

    }



    bool empty() const {
      return lower_map.empty() && upper_map.empty();
    }

    size_type size() const {
      return lower_map.size() + upper_map.size();
    }

    size_type unique_size() const {
      return lower_map.size() + upper_map.size();
    }

    void clear() {
      lower_map.clear_no_resize();
      upper_map.clear_no_resize();
    }

    void resize(size_t const n) {
      lower_map.resize(static_cast<float>(n) / 2.0 / lower_map.max_load_factor() );
      upper_map.resize(static_cast<float>(n) / 2.0 / upper_map.max_load_factor() );

    }

    /// rehash for new count number of BUCKETS.  iterators are invalidated.
    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return lower_map.bucket_count() + upper_map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float load_factor() {
      return  static_cast<float>(size()) / static_cast<float>(bucket_count());
    }



    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {

        InputIt middle = partition_input(first, last);

        lower_map.resize(static_cast<float>(lower_map.size() + ::std::distance(first, middle)) / lower_map.max_load_factor() ) ;
        lower_map.insert(first, middle);

        upper_map.resize(static_cast<float>(upper_map.size() + ::std::distance(middle, last)) / upper_map.max_load_factor()) ;
        upper_map.insert(middle, last);
    }

    /// inserting a vector
    void insert(::std::vector<::std::pair<Key, T> > & input) {
    	insert(input.begin(), input.end());
    }

    std::pair<typename container_type::iterator, bool> insert(::std::pair<Key, T> const & x) {
      if (splitter(x.first)) {
        return lower_map.insert(x);
      }
      else {
        return upper_map.insert(x);
      }
    }

    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");


      if (first == last) return 0;

      size_t count = 0;


      InputIt middle = partition_input(first, last);


      // mark for erasure
      for (; first != middle; ++first) {
        auto iter = lower_map.find(*(first));
        if (iter == lower_map.end()) continue;

        if (pred(*iter)) {
        	lower_map.erase(iter);
        	++count;
        }
      }

      for (; first != last; ++first) {
        auto iter = upper_map.find(*(first));
        if (iter == upper_map.end()) continue;

        if (pred(*iter)) {
        	upper_map.erase(iter);
        	++count;
        }
      }

      return count;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {

        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");

        if (first == last) return 0;

        size_t count = 0;


        InputIt middle = partition_input(first, last);


        // mark for erasure
        for (; first != middle; ++first) {
          auto iter = lower_map.find(*(first));
          if (iter == lower_map.end()) continue;

          lower_map.erase(iter);
          ++count;
        }

        for (; first != last; ++first) {
          auto iter = upper_map.find(*(first));
          if (iter == upper_map.end()) continue;

          	upper_map.erase(iter);
          	++count;
        }


        return count;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {
    	size_t before = size();

    	for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
    		if (pred(*it))
    				lower_map.erase(it);
    	}
    	for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
    		if (pred(*it))
    				upper_map.erase(it);
    	}

    	return before - size();
    }

    size_type count(Key const & key) const {
    	if (splitter(key))
    		return lower_map.count(key);
    	else
    		return upper_map.count(key);
    }


    container_range equal_range(Key const & key) {
    	if (splitter(key)) {
    		return lower_map.equal_range(key);
    	}
    	else {
    		return upper_map.equal_range(key);
    	}

    }
    container_const_range equal_range(Key const & key) const {
    	if (splitter(key)) {
    		return lower_map.equal_range(key);
    	}
    	else {
    		return upper_map.equal_range(key);
    	}
    }
    // NO bucket interfaces

};


// Key values does not span entire key space.
template <typename Key,
typename T,
typename Hash,
typename Equal,
typename Allocator>
class densehash_map<Key, T, ::fsc::TruePredicate, Hash, Equal, Allocator> {

  protected:

    using container_type =
        ::google::dense_hash_map<Key, T,
                       Hash, Equal, Allocator >;

    container_type map;



  public:
    using key_type              = Key;
    using mapped_type           = T;
    using value_type            = ::std::pair<const Key, T>;
    using hasher                = Hash;
    using key_equal             = Equal;
    using allocator_type        = Allocator;
    using reference             = value_type&;
    using const_reference       = const value_type&;
    using pointer               = typename std::allocator_traits<Allocator>::pointer;
    using const_pointer         = typename std::allocator_traits<Allocator>::const_pointer;
    using iterator              = typename container_type::iterator;
    using const_iterator        = typename container_type::const_iterator;
    using size_type             = size_t;
    using difference_type       = ptrdiff_t;


    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
    densehash_map(size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         map(bucket_count, hash, equal, alloc) {
    };

    densehash_map(Key empty_key, Key deleted_key,
                   size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         map(bucket_count, hash, equal, alloc) {
    	reserve_keys(empty_key, deleted_key);

//      map.max_load_factor(0.8);
//      map.min_load_factor(0.0);

    };

    template<class InputIt>
    densehash_map(InputIt first, InputIt last,
                      Key empty_key, Key deleted_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       densehash_map(empty_key, deleted_key, std::distance(first, last), hash, equal, alloc) {

        insert(first, last);
    };

    virtual ~densehash_map() {};


    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
        map.set_empty_key(empty_key);
        map.set_deleted_key(deleted_key);
    }
    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key, const ::fsc::TruePredicate & _splitter = fsc::TruePredicate()) {
    	//printf("reserve_upper_keys no-op\n");
    }



    iterator begin() {
      return map.begin();
    }
    const_iterator begin() const {
      return cbegin();
    }
    const_iterator cbegin() const {
      return map.begin();
    }

    iterator end() {
      return map.begin();
    }
    const_iterator end() const {
      return cend();
    }
    const_iterator cend() const {
      return map.begin();
    }


    std::vector<Key> keys() const {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const {
      ks.clear();
      ks.reserve(size());

      for (auto it = map.begin(); it != map.end(); ++it) {
        ks.emplace_back(it->first);
      }
    }

    std::vector<std::pair<Key, T> > to_vector() const {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T> > & vs) const {
      vs.clear();
      vs.reserve(size());

     for (auto it = map.begin(); it != map.end(); ++it) {
        vs.emplace_back(*it);
      }
    }


    bool empty() const {
      return map.empty();
    }

    size_type size() const {
      return map.size();
    }
    size_type unique_size() const {
      return map.size();
    }

    void clear() {
      map.clear_no_resize();
    }

    void resize(size_t const n) {
      map.resize(static_cast<float>(n) / map.max_load_factor());
    }

    /// rehash for new count number of BUCKETS.  iterators are invalidated.
    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float load_factor() {
      return  static_cast<float>(map.size()) / static_cast<float>(map.bucket_count());
    }


    std::pair<iterator, bool> insert(const value_type& x) {
      return this->map.insert(std::forward<value_type>(x));
    }


    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {
      this->resize(map.size() + std::distance(first, last));

      map.insert(first, last);
    }

    /// inserting sorted range
    void insert(::std::vector<::std::pair<Key, T> > & input) {
      insert(input.begin(), input.end());
    }

    std::pair<iterator, bool> insert(::std::pair<Key, T> const & x) {
      return map.insert(x);
    }


    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;

      size_t count = 0;

      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*(first));
        if (iter == map.end()) continue;

        if (pred(*iter)) {
          map.erase(iter);
          ++count;
        }
      }
      return count;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");

        if (first == last) return 0;

        size_t count = 0;

        // mark for erasure
        for (; first != last; ++first) {
          auto iter = map.find(*(first));
          if (iter == map.end()) continue;

            map.erase(iter);
            ++count;
        }
        return count;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {
      size_t before = map.size();

      for (auto it = map.begin(); it != map.end(); ++it) {
        if (pred(*it))
            map.erase(it);
      }

      return before - map.size();
    }

    size_type count(Key const & key) const {
      return map.count(key);
    }


    ::std::pair<iterator, iterator> equal_range(Key const & key) {
      return map.equal_range(key);
    }

    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
      return map.equal_range(key);
    }
    // NO bucket interfaces

};







/**
 * @brief multimap implemented using densehash.
 * @details:  earlier attempt used a sorted vector.  while still faster than unordered map, the sorting step does not scale well.
 *            a second attempt tried to use key+counter as key, hash only on key, compare on full key+counter, and query via bucket interface and compare only key portion.
 *              this did not work as google dense hash map bucket points to array element, and does not represent a collection of entries with same hash value like
 *              std::unordered_multimap.
 *
 *            current implementation uses 1 vector for singletons, and 1 vector of vectors for multiple entries.  memory utilization is probably suboptimal (up to 2x)
 *            map stores index in these vectors, with sign bit signifying that it is a multiple entry.
 *
 *            this trades off sorting with more vector allocations and copying, but appears to be faster at least for the case when there are large number of singleton entries.
 *
 */
template <typename Key,
typename T,
typename LowerKeySpaceSelector,
typename Hash = ::std::hash<Key>,
typename Equal = ::std::equal_to<Key>,
typename Allocator = ::std::allocator<::std::pair<Key, T> > >
class densehash_multimap {

  protected:

    // data container
    using subcontainer_type = ::std::vector<::std::pair<Key, T>, Allocator >;
    using subiter_type = typename ::std::vector<::std::pair<Key, T>, Allocator >::iterator;
    using const_subiter_type = typename ::std::vector<::std::pair<Key, T>, Allocator >::const_iterator;

    // index in the vector - this is so we don't have to worry about pointer or iterators being invalidated when vector resizes
    // non-negative values indicate singletons.  negative values indicate multiple entries.
    // multiple entries is stored in a vector of vector.  the index is internal_val_type with sign bit removed.
    using internal_val_type = int64_t;


    using supercontainer_type =
        ::google::dense_hash_map<Key, internal_val_type,
                       Hash, Equal, Allocator >;

    subcontainer_type vec1;
    std::vector<subcontainer_type, Allocator> vecX;
    supercontainer_type lower_map;
    supercontainer_type upper_map;
    size_t s;

    // TODO: provide iterator implementation for  begin/end.

    LowerKeySpaceSelector splitter;

    template <typename InputIt>
    InputIt partition_input(InputIt first, InputIt last) {
      return ::std::stable_partition(first, last, splitter);
    }



    /**
     * @class    bliss::iterator::ConcatenatingIterator
     * @brief    this class presents a single/sequential view of a series of underlying iterator ranges
     * @details  random access iterator is not supported. other iterator categories are okay.
     *
     */
//    template<typename V>
//    class concat_iter :
//      public ::std::iterator<
//        typename ::std::forward_iterator_tag,
//        V
//      >
//    {
//      protected:
//        using superiterator_type = typename ::std::conditional<::std::is_const<V>::value,
//            typename supercontainer_type::const_iterator, typename supercontainer_type::iterator>::type;
//        using type = concat_iter<V>;
//
//      protected:
//        /// the current position in the ranges list
//
//        supercontainer_type & lmap;
//        supercontainer_type & umap;
//
//        superiterator_type curr_iter;
//
//        bool at_max;
//
//        /// enforce that iterator is at a dereferenceable position
//        void ensure_dereferenceable() {
//          if (at_max) return;
//
//          if (curr_iter == lmap.end()) {
//            curr_iter = umap.begin();
//          }
//
//          if (curr_iter == umap.end()) {
//            at_max = true;
//            return;
//          }
//
//          return;
//
//        }
//
//      public:
//        using difference_type = typename ::std::iterator_traits<superiterator_type>::difference_type;
//
//
//        /// constructor for end concat iterator.  _end refers to end of supercontainer.
//        concat_iter(supercontainer_type & _lmap, supercontainer_type & _umap) : lmap(_lmap), umap(_umap), curr_iter(_umap.end()), at_max(true) {};
//
//        /// constructor for end concat iterator.  _end refers to end of supercontainer.
//        concat_iter(supercontainer_type & _lmap, supercontainer_type & _umap, superiterator_type _iter) :
//          lmap(_lmap), umap(_umap), curr_iter(_iter) {
//          ensure_dereferenceable();
//        };
//
//        // note that explicit keyword cannot be on copy and move constructors else the constructors are not defined/found.
//
//        // copy constructor, assignment operator, move constructor, assignment operator should
//        // should be default since the member vars are simple.
//
//        bool is_at_max() const {
//          at_max = (curr_iter == umap.end());
//          return at_max;
//        }
//
//        /**
//         * @brief increment:  move to the next position in the concatenating iterator, which may cross range boundaries.
//         * @note  side effect: set at_end variable.
//         * @return
//         */
//        type& operator++() {
//          // if at end, return
//          if (!at_max) {
//            // now increment.  since we are careful to leave iterator at a dereferenceable state, curr_pos is not at a subcontainer's end.
//            // so just increment.
//            ++curr_iter;
//
//            // now make sure we don't end up at a subcontainer's end.
//            ensure_dereferenceable();
//          }
//          return *this;
//        }
//
//
//        /**
//         * post increment.  make a copy then increment that.
//         */
//        type operator++(int)
//        {
//          type output(*this);
//          this->operator++();
//          return output;
//        }
//
//        //=== input iterator specific
//
//        /// comparison operator
//        bool operator==(const type& rhs) const
//          {
//          if ((lmap != rhs.lmap) || (umap != rhs.umap)) throw std::logic_error("the iterators being compared do not have the same internal end iterators so they are not comparable.");
//
//          if (at_max && rhs.at_max) return true;
//          if (at_max || rhs.at_max) return false;
//
//          return (curr_iter == rhs.curr_iter);
//          }
//
//        /// comparison operator
//        bool operator!=(const type& rhs) const
//          {
//          if ((lmap != rhs.lmap) || (umap != rhs.umap)) throw std::logic_error("the iterators being compared do not have the same internal end iterators so they are not comparable.");
//
//          if (at_max && rhs.at_max) return false;
//          if (at_max || rhs.at_max) return true;
//
//          return (curr_iter != rhs.curr_iter);
//          }
//
//
//        inline V operator*() const {
//          return *curr_iter;
//        }
//
//        /*=== NOT output iterator.  this is a map, does not make sense to change the  */
//        /* content via iterator.                                                      */
//
//
//        //=== NOT full forward iterator - no default constructor
//
//        //=== NOT bidirectional iterator - no decrement because map produces forward iterator only.
//
//        //=== NOT full random access iterator - only have +, += but not -. -=.  no comparison operators. have offset dereference operator [].
//
//    };


  public:
    using key_type              = Key;
    using mapped_type           = T;
    using value_type            = ::std::pair<const Key, T>;
    using hasher                = Hash;
    using key_equal             = Equal;
    using allocator_type        = Allocator;
    using reference             = value_type&;
    using const_reference       = const value_type&;
    using pointer               = typename std::allocator_traits<Allocator>::pointer;
    using const_pointer         = typename std::allocator_traits<Allocator>::const_pointer;
    using iterator              = typename subcontainer_type::iterator;
    using const_iterator        = typename subcontainer_type::const_iterator;
    using size_type             = typename subcontainer_type::size_type;
    using difference_type       = typename subcontainer_type::difference_type;

  protected:
    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert_impl(InputIt first, InputIt last, supercontainer_type & map) {


        // get previous sizes so we know where to start from
        int64_t idx1 = vec1.size();
        int64_t idxX = vecX.size();

        // reserve but do not yet copy in.  we will have random access to check if already exist,
        // so don't copying in does not save a whole lot.
        vec1.reserve(vec1.size() + std::distance(first, last));

        // iterator over all and insert into map.
        std::pair<typename supercontainer_type::iterator, bool> insert_result;
        Key k;
        int64_t idx;
        for (InputIt it = first, max = last; it != max; ++it) {
          k = it->first;

          // try inserting
          insert_result = map.insert(std::make_pair(k, idx1));

          if (insert_result.second) {
            // successful map insertion, so now insert into vec1.
            vec1.emplace_back(*it);
            // new map entry already inserted. just need to increment.
            ++idx1;
          } else {
            // entry already there.  get the iterator and check the idx
            idx = insert_result.first->second;
            if (idx < 0) {
              // previously inserted multiple entries, so reuse the vector, and insert
              vecX[idx & ::std::numeric_limits<int64_t>::max()].emplace_back(*it);
              // no map update is needed.
            } else {
              // previously inserted 1 entry. so
              // update map with new vecX value
              insert_result.first->second = idxX | ~(::std::numeric_limits<int64_t>::max());
              // create a new entry in vecX
              vecX.emplace_back(subcontainer_type());
              // get previous value from vec1 and insert into vecX
              vecX[idxX].emplace_back(std::move(vec1[idx]));
              // insert new value into vecX
              vecX[idxX].emplace_back(*it);
              // update new vecX idx.
              ++idxX;
            }
          }
        }

        s += std::distance(first, last);

        // TODO: compact
    }



    template <typename InputIt, typename Pred>
    void erase_impl(InputIt first, InputIt last, Pred const & pred, supercontainer_type & map) {

      size_t dist = 0;
      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }

      }

    }



    template <typename InputIt>
    void erase_impl(InputIt first, InputIt last, supercontainer_type & map) {

      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {
          // multi.
          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          s -= vec.size();
          subcontainer_type().swap(vec);
        } else {
          // else single.  do nothing.
          --s;
        }

        // now clear the map entry.
        map.erase(iter);
      }

    }


    template <typename Pred>
    void erase_impl(Pred const & pred, supercontainer_type & map) {

      size_t dist = 0;
      int64_t idx;


      // mark for erasure
      auto max = map.end();
      for (auto iter = map.begin(); iter != max; ++iter) {
        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }
      }
    }


    size_type count_impl(Key const & key, supercontainer_type const & map) const {
      auto iter = map.find(key);
      if (iter == map.end()) return 0;

      if (iter->second < 0) {
        // multiple entries
        return vecX[iter->second &  ::std::numeric_limits<int64_t>::max()].size();
      } else {
        return 1;
      }
    }


    ::std::pair<iterator, iterator> equal_range_impl(Key const & key, supercontainer_type & map) {


      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.end(), vec1.end());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.begin() + iter->second, vec1.begin() + iter->second + 1);

      // found, has multiple values
      subcontainer_type & vec = vecX[iter->second & ::std::numeric_limits<int64_t>::max()];

      return std::make_pair(vec.begin(), vec.end());

    }
    ::std::pair<const_iterator, const_iterator> equal_range_impl(Key const & key, supercontainer_type const & map) const {
      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.cend(), vec1.cend());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.cbegin() + iter->second, vec1.cbegin() + iter->second + 1);

      // found, has multiple values
      subcontainer_type & vec = vecX[iter->second & ::std::numeric_limits<int64_t>::max()];

      return std::make_pair(vec.cbegin(), vec.cend());

    }

  public:


    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
    densehash_multimap(size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         lower_map(bucket_count / 2, hash, equal, alloc),
                         upper_map(bucket_count / 2, hash, equal, alloc), s(0UL) {
//      lower_map.max_load_factor(0.8);
//      lower_map.min_load_factor(0.0);
//      upper_map.max_load_factor(0.8);
//      upper_map.min_load_factor(0.0);

    };

    densehash_multimap(Key empty_key, Key deleted_key,
                             Key upper_empty_key, Key upper_deleted_key,
                             const LowerKeySpaceSelector & _splitter = LowerKeySpaceSelector(),
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         densehash_multimap(bucket_count / 2, hash, equal, alloc) {
    	reserve_keys(empty_key, deleted_key);
    	reserve_upper_keys(upper_empty_key, upper_deleted_key, _splitter);
    };

    template<class InputIt>
    densehash_multimap(InputIt first, InputIt last,
                      Key empty_key, Key deleted_key,
                      Key upper_empty_key, Key upper_deleted_key,
                       size_type bucket_count = 128,
                       const LowerKeySpaceSelector & _splitter = LowerKeySpaceSelector(),
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       densehash_multimap(empty_key, deleted_key, upper_empty_key, upper_deleted_key, _splitter, bucket_count, hash, equal, alloc) {
        this->insert(first, last);
    };

    virtual ~densehash_multimap() {
      std::cout << " singles: " << vec1.size() << " multiples: " << vecX.size() << std::endl;
    };


    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
        lower_map.set_empty_key(empty_key);
        lower_map.set_deleted_key(deleted_key);
    }

    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key, const LowerKeySpaceSelector & _splitter = LowerKeySpaceSelector()) {
    	printf("reserve_upper_keys\n");
        upper_map.set_empty_key(empty_key);
        upper_map.set_deleted_key(deleted_key);
        splitter = _splitter;
    }


    // TODO: begin and end iterator accessors.


    std::vector<Key> keys() const  {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const  {
      ks.clear();
      ks.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        ks.emplace_back(it->first);
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        ks.emplace_back(it->first);
      }

      return ks;
    }

    std::vector<std::pair<Key, T>> to_vector() const  {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T>> & vs) const  {
      vs.clear();
      vs.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        if (it->second < 0) {
          auto it2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].begin();
          auto max2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].end();
          for (; it2 != max2; ++it2) {
            vs.emplace_back(*it2);
          }
        } else {  // singleton
          vs.emplace_back(vec1[it->second]);
        }
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        if (it->second < 0) {
          auto it2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].begin();
          auto max2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].end();
          for (; it2 != max2; ++it2) {
            vs.emplace_back(*it2);
          }
        } else {  // singleton
          vs.emplace_back(vec1[it->second]);
        }
      }
    }


    bool empty() const {
      return lower_map.empty() && upper_map.empty();
    }

    size_type size() const {
      return s;
    }


    size_type unique_size() const {
      return lower_map.size() + upper_map.size();
    }

    void clear() {
      vec1.clear();
      vecX.clear();
      lower_map.clear_no_resize();
      upper_map.clear_no_resize();
      s = 0UL;
    }

    void resize(size_t const n) {
      lower_map.resize(static_cast<float>(n) / 2.0 / lower_map.max_load_factor());
      upper_map.resize(static_cast<float>(n) / 2.0 / lower_map.max_load_factor());
    }

    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return lower_map.bucket_count() + upper_map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float max_load_factor() {
      return static_cast<float>(size()) / static_cast<float>(bucket_count());
    }





    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<std::pair<Key, T>,
                      typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for insert cannot be converted to std::pair<Key, T> type");

        if (first == last) return;

        auto middle = partition_input(first, last);

        lower_map.resize(static_cast<float>(lower_map.size() + std::distance(first, middle)) / lower_map.max_load_factor());
        insert_impl(first, middle, lower_map);
        upper_map.resize(static_cast<float>(upper_map.size() + std::distance(first, middle)) / upper_map.max_load_factor());
        insert_impl(middle, last, upper_map);

        // TODO: compact
    }

    /// inserting sorted range
    void insert(::std::vector<::std::pair<Key, T> > & input) {
      // TODO: more memory efficient version of this.

      insert(input.begin(), input.end());
    }

    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;


      size_t before = s;

      auto middle = partition_input(first, last);
      erase_impl(first, middle, pred, lower_map);
      erase_impl(middle, last, pred, upper_map);

      return before - s;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;


      size_t before = s;

      auto middle = partition_input(first, last);
      erase_impl(first, middle, lower_map);
      erase_impl(middle, last, upper_map);

      return before - s;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {

      if (this->size() == 0) return 0;

      size_t before = s;

      erase_impl(pred, lower_map);
      erase_impl(pred, upper_map);

      return before - s;
    }

    size_type count(Key const & key) const {

      if (splitter(key)) {
        return count_impl(key, lower_map);
      } else {
        return count_impl(key, upper_map);
      }

    }



    ::std::pair<iterator, iterator> equal_range(Key const & key) {
      if (splitter(key)) {
        return equal_range_impl(key, lower_map);
      } else {
        return equal_range_impl(key, upper_map);
      }
    }
    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
      if (splitter(key)) {
        return equal_range_impl(key, lower_map);
      } else {
        return equal_range_impl(key, upper_map);
      }
    }
    // NO bucket interfaces

};




/**
 * specialization for densehash that does not require key space partitioning
 */
template <typename Key,
typename T,
typename Hash,
typename Equal,
typename Allocator>
class densehash_multimap<Key, T, ::fsc::TruePredicate, Hash, Equal, Allocator> {

  protected:

    // data container
    using tuple_allocator_type = typename Allocator::template rebind< ::std::pair<Key, T> >::other;
    using subcontainer_type = ::std::vector<::std::pair<Key, T>, tuple_allocator_type >;
    using subiter_type = typename subcontainer_type::iterator;
    using const_subiter_type = typename subcontainer_type::const_iterator;

    // index in the vector - this is so we don't have to worry about pointer or iterators being invalidated when vector resizes
    // non-negative values indicate singletons.  negative values indicate multiple entries.
    // multiple entries is stored in a vector of vector.  the index is internal_val_type with sign bit removed.
    using internal_val_type = int64_t;

    using super_allocator_type = typename Allocator::template rebind<std::pair< const Key, internal_val_type> >::other;
    using supercontainer_type =
        ::google::dense_hash_map<Key, internal_val_type,
                       Hash, Equal, super_allocator_type >;

    subcontainer_type vec1;
    using vector_allocator_type = typename Allocator::template rebind< subcontainer_type >::other;
    std::vector<subcontainer_type, vector_allocator_type> vecX;
    supercontainer_type map;
    size_t s;

    // TODO: provide iterator implementation for  begin/end.


  public:
    using key_type              = Key;
    using mapped_type           = T;
    using value_type            = ::std::pair<const Key, T>;
    using hasher                = Hash;
    using key_equal             = Equal;
    using allocator_type        = Allocator;
    using reference             = value_type&;
    using const_reference       = const value_type&;
    using pointer               = typename std::allocator_traits<Allocator>::pointer;
    using const_pointer         = typename std::allocator_traits<Allocator>::const_pointer;
    using iterator              = typename subcontainer_type::iterator;
    using const_iterator        = typename subcontainer_type::const_iterator;
    using size_type             = typename subcontainer_type::size_type;
    using difference_type       = typename subcontainer_type::difference_type;


    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
    densehash_multimap(size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         map(bucket_count, hash, equal, alloc), s(0UL) {
//      map.max_load_factor(0.8);
//      map.min_load_factor(0.0);
    };

    densehash_multimap(Key empty_key, Key deleted_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         densehash_multimap(bucket_count, hash, equal, alloc) {
    	reserve_keys(empty_key, deleted_key);
    };

    template<class InputIt>
    densehash_multimap(InputIt first, InputIt last,
                      Key empty_key, Key deleted_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       densehash_multimap(empty_key, deleted_key, bucket_count, hash, equal, alloc) {
        this->insert(first, last);
    };

    virtual ~densehash_multimap() {
      std::cout << " singles: " << vec1.size() << " multiples: " << vecX.size() << std::endl;
    };



    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
        map.set_empty_key(empty_key);
        map.set_deleted_key(deleted_key);
    }

    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key, const ::fsc::TruePredicate & _splitter = fsc::TruePredicate()) {
    	printf("reserve_upper_keys no-op\n");
    }


    // TODO: begin and end iterator accessors.
    std::vector<Key> keys() const  {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const  {
      ks.clear();
      ks.reserve(size());

      for (auto it = map.begin(); it != map.end(); ++it) {
        ks.emplace_back(it->first);
      }
    }

    std::vector<std::pair<Key, T>> to_vector() const  {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T>> & vs) const  {
      vs.clear();
      vs.reserve(size());

      int64_t idx;
      for (auto it = map.begin(); it != map.end(); ++it) {
        if (it->second < 0) {
          idx = it->second & ::std::numeric_limits<int64_t>::max();
          auto max2 = vecX[idx].end();
          for (auto it2 = vecX[idx].begin(); it2 != max2; ++it2) {
            vs.emplace_back(*it2);
          }
        } else {  // singleton
          vs.emplace_back(vec1[it->second]);
        }
      }
    }



    bool empty() const {
      return map.empty();
    }

    size_type size() const {
      return s;
    }

    size_type unique_size() const {
      return map.size();
    }

    void clear() {
      vec1.clear();
      vecX.clear();
      map.clear_no_resize();
      s = 0UL;
    }

    void resize(size_t const n) {
      map.resize(static_cast<float>(n) / map.max_load_factor());
    }

    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float max_load_factor() {
      return static_cast<float>(map.size()) / static_cast<float>(map.bucket_count());
    }



    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<std::pair<Key, T>,
                      typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for insert cannot be converted to std::pair<Key, T> type");

        if (first == last) return;


        // get previous sizes so we know where to start from
        int64_t idx1 = vec1.size();
        int64_t idxX = vecX.size();

        // reserve but do not yet copy in.  we will have random access to check if already exist,
        // so don't copying in does not save a whole lot.
        vec1.reserve(vec1.size() + std::distance(first, last));
        this->resize(vec1.size() + std::distance(first, last));

        // iterator over all and insert into map.
        std::pair<typename supercontainer_type::iterator, bool> insert_result;
        Key k;
        int64_t idx;
        for (InputIt it = first, max = last; it != max; ++it) {
          k = it->first;

          // try inserting
          insert_result = map.insert(std::make_pair(k, idx1));

          if (insert_result.second) {
            // successful map insertion, so now insert into vec1.
            vec1.emplace_back(*it);
            // new map entry already inserted. just need to increment.
            ++idx1;
          } else {
            // entry already there.  get the iterator and check the idx
            idx = insert_result.first->second;
            if (idx < 0) {
              // previously inserted multiple entries, so reuse the vector, and insert
              vecX[idx & ::std::numeric_limits<int64_t>::max()].emplace_back(*it);
              // no map update is needed.
            } else {
              // previously inserted 1 entry. so
              // update map with new vecX value
              insert_result.first->second = idxX | ~(::std::numeric_limits<int64_t>::max());
              // create a new entry in vecX
              vecX.emplace_back(subcontainer_type());
              // get previous value from vec1 and insert into vecX
              vecX[idxX].emplace_back(std::move(vec1[idx]));
              // insert new value into vecX
              vecX[idxX].emplace_back(*it);
              // update new vecX idx.
              ++idxX;
            }
          }
        }

        s += std::distance(first, last);

        // TODO: compact
    }

    /// inserting sorted range
    void insert(::std::vector<::std::pair<Key, T> > & input) {

      if (input.size() == 0) return;
//
//      if (this->size() > 0) {
        // there is data in there already.  so do normal insert
        this->insert(input.begin(), input.end());
//        return;
//      }
//
//      // else empty, so we can improve the memory utilization a little.
//
//      // get previous sizes so we know where to start from
//      input.swap(vec1);
//      int64_t idx1 = 0;
//
//      // iterator over all and insert into map.
//      std::pair<typename supercontainer_type::iterator, bool> insert_result;
//      Key k;
//      int64_t idx;
//      int64_t max = vec1.size();
//      for (int64_t i = 0; i < max; ++i) {
//        k = vec1[i].first;
//
//        // try inserting
//        insert_result = map.insert(std::make_pair(k, idx1));
//
//        if (insert_result.second) {
//          // new map entry already inserted. just need to increment.
//
//          // and also move the entry to the right place
//          if (i != idx1) ::std::swap(vec1[idx1], vec1[i]);
//
//          ++idx1;
//        } else {
//          // entry already there.  get the iterator and check the idx
//          idx = insert_result.first->second;
//          if (idx < 0) {
//            // previously inserted multiple entries, so reuse the vector, and insert
//            vecX[idx & ::std::numeric_limits<int64_t>::max()].emplace_back(vec1[i]);
//            // no map update is needed.
//          } else {
//            // previously inserted 1 entry. so
//            // update map with new vecX value
//            insert_result.first->second = static_cast<int64_t>(vecX.size()) | ~(::std::numeric_limits<int64_t>::max());
//            // create a new entry in vecX
//            vecX.emplace_back(subcontainer_type());
//            // get previous value from vec1 and insert into vecX
//            vecX.back().emplace_back(std::move(vec1[idx]));
//            // insert new value into vecX
//            vecX.back().emplace_back(vec1[i]);
//          }
//        }
//      }
//
//      s += vec1.size();
//
//      vec1.erase(vec1.begin() + idx1, vec1.end());
//      vec1.resize(idx1);
    }

    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;

      size_t before = s;
      size_t dist = 0;
      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }

      }

      return before - s;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;

      size_t before = s;

      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {
          // multi.
          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          s -= vec.size();
          subcontainer_type().swap(vec);
        } else {
          // else single.  do nothing.
          --s;
        }

        // now clear the map entry.
        map.erase(iter);
      }

      return before - s;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {

      if (this->size() == 0) return 0;

      size_t before = s;
      size_t dist = 0;
      int64_t idx;


      // mark for erasure
      auto  max = map.end();
      for (auto iter = map.begin(); iter != max; ++iter) {
        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }
      }

      return before - s;
    }

    size_type count(Key const & key) const {
      auto iter = map.find(key);
      if (iter == map.end()) return 0;

      if (iter->second < 0) {
        // multiple entries
        return vecX[iter->second &  ::std::numeric_limits<int64_t>::max()].size();
      } else {
        return 1;
      }
    }



    ::std::pair<iterator, iterator> equal_range(Key const & key) {
      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.end(), vec1.end());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.begin() + iter->second, vec1.begin() + iter->second + 1);

      // found, has multiple values
      subcontainer_type & vec = vecX[iter->second & ::std::numeric_limits<int64_t>::max()];

      return std::make_pair(vec.begin(), vec.end());

    }
    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.cend(), vec1.cend());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.cbegin() + iter->second, vec1.cbegin() + iter->second + 1);

      // found, has multiple values
      int64_t idx = iter->second & ::std::numeric_limits<int64_t>::max();

      return std::make_pair(vecX[idx].cbegin(), vecX[idx].cend());


    }
    // NO bucket interfaces

};







  /**
   * uncompacted version of the vecmap.  DEPRECATED AND NOT USED....
   *
   * internally has a single vector, sorted by (hash(key) % buckets) then by hash(key) then by key
   * insert into back of vector.   when insert, when map load_factor may trigger a rehash, we manually resort the vector and rebuild the map.
   *
   * note that as we are attempting to adhere to unordered_multimap template interface, we can't really add a comparacter parameter without causing
   * problems elsewhere, namely in distributed unordered hashvec map.  the transform that is inherently in the hash needs to be applied, else we might
   * get into a situation where std::less results in 2 separate keys while hash lumps the keys together.
   *
   * we also cannot sort by hash value, because of possibility of collision.
   */
  template <typename Key,
  typename T,
  typename Hash = ::std::hash<Key>,
  typename Comparator = ::std::less<Key>,
  typename Equal = ::std::equal_to<Key>,
  typename Allocator = ::std::allocator<::std::pair<Key, T> > >
  class densehash_vecmap {

    protected:
      struct Less {
        Comparator l;

        inline bool operator()(Key const &x, Key const &y ) {
          return l(x, y);
        }

        template <typename V>
        inline bool operator()(::std::pair<Key, V> const & x, ::std::pair<Key, V> const & y) {
          return l(x.first, y.first);
        }
        template <typename V>
        inline bool operator()(::std::pair<const Key, V> const & x, ::std::pair<const Key, V> const & y) {
          return l(x.first, y.first);
        }
      };


      // data container
      using subcontainer_type = ::std::vector<::std::pair<Key, T>, Allocator >;
      using subiter_type = typename ::std::vector<::std::pair<Key, T>, Allocator >::iterator;
      using const_subiter_type = typename ::std::vector<::std::pair<Key, T>, Allocator >::const_iterator;

      // index "pointers"
      using value_range_type = ::std::pair<typename subcontainer_type::iterator,
           typename subcontainer_type::iterator>;

      using superallocator_type = ::std::allocator<::std::pair<const Key, value_range_type > >;
      using supercontainer_type =
          ::google::dense_hash_map<Key, value_range_type,
                         Hash, Equal, superallocator_type >;

      subcontainer_type vec;
      supercontainer_type map;
      size_t s;


      // use of vector - increased cost during construction but query will be fast.
      // group all entries with the same key with vector - list chaing but already sorted.
      // group all entries with same hash with vector - list chaining but with randomly accessible container.

      inline size_t dist(value_range_type const & iters) const {
        return std::distance(iters.first, iters.second);
      }
      inline size_t dist(Key const & k) const {
        return dist(map.find(k));
      }
      inline size_t dist(typename supercontainer_type::iterator map_iter) {
        return (map_iter == map.end()) ? 0 : dist(map_iter->second);
      }
      inline size_t dist(typename supercontainer_type::const_iterator map_iter) const {
        return (map_iter == map.end()) ? 0 : dist(map_iter->second);
      }

      /**
       * @class    bliss::iterator::ConcatenatingIterator
       * @brief    this class presents a single/sequential view of a series of underlying iterator ranges
       * @details  random access iterator is not supported. other iterator categories are okay.
       *
       */
      template<typename V>
      class concat_iter :
        public ::std::iterator<
          typename ::std::random_access_iterator_tag,
          V
        >
      {
        protected:
          using subiterator_type = typename ::std::conditional<::std::is_const<V>::value,
              typename subcontainer_type::const_iterator, typename subcontainer_type::iterator>::type;
          using superiterator_type = typename ::std::conditional<::std::is_const<V>::value,
              typename supercontainer_type::const_iterator, typename supercontainer_type::iterator>::type;
          using type = concat_iter<V>;

          using inner_value_type = typename ::std::iterator_traits<subiterator_type>::value_type;

        public:
          template <typename KK, typename TT, typename HH, typename EE, typename AA, typename OutputIterator>
          OutputIterator
          copy(typename ::fsc::densehash_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> first,
               typename ::fsc::densehash_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> last, OutputIterator result);


        protected:
          /// the current position in the ranges list
          superiterator_type curr_iter;

          /// the current iterator position in the range of interest.
          subiterator_type curr_pos;

          superiterator_type max_iter;

          bool at_max;

          /// enforce that iterator is at a dereferenceable position
          void ensure_dereferenceable() {
            if (at_max) return;

            if (curr_iter == max_iter) {
              at_max = true;
              return;
            }

            // check to see if we are at end of subcontainer.  if so, move to next dereferenceable position.
            // end of a subcontainer is treated as same position as beginning of next subcontainer.
            while (curr_pos == curr_iter->second.second) {
              ++curr_iter;
              if (curr_iter == max_iter) {
                at_max = true;
                break; // reached the very end, can't deref max_iter.
              }
              else curr_pos = curr_iter->second.first;
            }
          }

        public:
          using difference_type = typename ::std::iterator_traits<subiterator_type>::difference_type;


          /// constructor for end concat iterator.  _end refers to end of supercontainer.
          concat_iter(superiterator_type _end) : curr_iter(_end), max_iter(_end), at_max(true) {};


          /// constructor for start concatenating iterator.  general version
          concat_iter(superiterator_type _iter, superiterator_type _end, subiterator_type _pos) :
            curr_iter(_iter), curr_pos(_pos), max_iter(_end), at_max(_iter == _end) {
            ensure_dereferenceable();
          };

          /// constructor for start concatenating iterator.  general version with checking that _pos belongs to _iter subcontainer via distance check.
          concat_iter(superiterator_type _iter, superiterator_type _end, subiterator_type _pos, difference_type distance_check) :
            concat_iter<V>(_iter, _end, _pos) {
            if (!at_max && (distance_check != ::std::distance(_iter->second.first, _pos) ) )
              throw std::logic_error("unordered_compact_vecmap constructor failing distance check, suggesting that _pos is not from same subcontainer as what _iter points to");
            ensure_dereferenceable();
          };


          // note that explicit keyword cannot be on copy and move constructors else the constructors are not defined/found.

          // copy constructor, assignment operator, move constructor, assignment operator should
          // should be default since the member vars are simple.

          bool is_at_max() const {
            at_max = (curr_iter == max_iter);
            return at_max;
          }

          /**
           * @brief increment:  move to the next position in the concatenating iterator, which may cross range boundaries.
           * @note  side effect: set at_end variable.
           * @return
           */
          type& operator++() {
            // if at end, return
            if (!at_max) {
              // now increment.  since we are careful to leave iterator at a dereferenceable state, curr_pos is not at a subcontainer's end.
              // so just increment.
              ++curr_pos;

              // now make sure we don't end up at a subcontainer's end.
              ensure_dereferenceable();
            }
            return *this;
          }


          /**
           * post increment.  make a copy then increment that.
           */
          type operator++(int)
          {
            type output(*this);
            this->operator++();
            return output;
          }

          //=== input iterator specific

          /// comparison operator
          bool operator==(const type& rhs) const
            {
            if (max_iter != rhs.max_iter) throw std::logic_error("the iterators being compared do not have the same internal end iterators so they are not comparable.");

            if (at_max && rhs.at_max) return true;
            if (at_max || rhs.at_max) return false;

            return ((curr_iter == rhs.curr_iter) && (curr_pos == rhs.curr_pos));
            }

          /// comparison operator
          bool operator!=(const type& rhs) const
            {
            if (max_iter != rhs.max_iter) throw std::logic_error("the iterators being compared do not have the same internal end iterators so they are not comparable.");

            if (at_max && rhs.at_max) return false;
            if (at_max || rhs.at_max) return true;

            return ((curr_iter != rhs.curr_iter) || (curr_pos != rhs.curr_pos));
            }


          template <
              typename VV = V,
              typename IV = inner_value_type,
              typename = typename ::std::enable_if<::std::is_constructible<VV, IV>::value>::type>
          inline V operator*() const {
            return *curr_pos;
          }

          /*=== NOT output iterator.  this is a map, does not make sense to change the  */
          /* content via iterator.                                                      */


          //=== NOT full forward iterator - no default constructor

          //=== NOT bidirectional iterator - no decrement because map produces forward iterator only.

          //=== NOT full random access iterator - only have +, += but not -. -=.  no comparison operators. have offset dereference operator [].


          /**
           * @brief     Advances this iterator by `n` positions.
           *            used by std::advance when randomaccess iterator.
           * @param n   The number of positions to advance.
           * @return    A reference to this after advancing.
           */
          type& operator+=(difference_type n)
          {
            // ::std::advance will use this or the ++ operator.
            if (n < 0) throw ::std::logic_error("::fsc::densehash_vecmap::iterator does not support decrement.");
            if (n == 0) return *this;  // nothing to add.
            if (at_max) return *this;  // iterator at the end.

            auto orig_iter = curr_iter;

            // dereferenceable right now
            auto curr_dist = ::std::distance(curr_pos, curr_iter->second.second);
            while (n >= curr_dist) {
              // not at end, and n is larger than curr dist, so go to next subcontainer.
              n -= curr_dist;  // consume some entries
              ++curr_iter;     // go to next container.
              at_max = (curr_iter == max_iter);
              if (at_max) return *this;  // if we are at end right now, then we can just return.
              else curr_dist = dist(curr_iter);    // see how much the next container has.

            }  // when exiting here, we are at a subcontainer that has more entries than n.  n could be 0.

            // now reset the curr_pos if curr_iter has been moved.
            if (curr_iter != orig_iter) curr_pos = curr_iter->second.first;
            ::std::advance(curr_pos, n);

            return *this;
          }


          /**
           * @brief     Advances a copy of this iterator by `n` positions.
           *
           * @param n   The number of positions to advance.
           * @return    The advanced iterator.
           */
          type operator+(difference_type n)
          {
            // reduced to += operator
            type output(*this);
            output += n;
            return output;
          }

          /**
           * @brief     Advances a copy of the `right` iterator by `n` positions.
           *
           * @param n   The number of positions to advance.
           * @return    The advanced iterator.
           */
          friend type operator+(difference_type n, const type& right)
          {
            // reduced to + operator
            return right + n;
          }

          /**
           * @brief     Returns the n'th element as seen from the current iterator
           *            position.
           *
           * @param n   The offset.
           *
           * @return    The element at offset `n` from the current position.
           */
          V operator[](difference_type n)
          {
            // reduce to the following:
            return *(*this + n);
          }

          /// difference between 2 iterators.  used by std::distance.
          friend difference_type operator-(const type& last, const type& first) {
            if (last == first) return 0;  // if both are at end, then we say dist is 0.
            if (first.at_max) return ::std::numeric_limits<difference_type>::lowest();  // first is at end, can't get to last.


            // now try to increment first until we get to last, or fail at getting to last.
            difference_type n = 0;
            auto cit = first.curr_iter;

            // init distance.  only meaningful here when first and last are not on same subcontainer.
            auto dist = std::distance(first.curr_pos, cit->second.second);

            // walk until either we are in same subcontainer, or at end of first iterator.
            while (cit != last.curr_iter) {
              n += dist;  //
              ++cit;
              if (cit == first.max_iter) break;
              else dist = dist(cit->second);
            }

            // at this point, we have cit == last.curr_iter, or cit == eit (cit == eit == last_curr_iter possible)
            if (cit == first.max_iter) {  // cit = eit
              // if at end of first, but not at curr of last, not reachable.
              // else if at end of first, and at curr of last (== end), then reached.  return n.
              return ((cit != last.curr_iter) ? ::std::numeric_limits<difference_type>::lowest() : n);
            }

            // else we have cit == last.curr_iter.  dist is either size of cit subcontainer, or from curr_pos to subcontainer's end, both not correct.
            // need to recalculate dist now.  first move cpos

            // recalc distance. if cit hasn't moved, use original.  else use beginning of current subcontainer.
            dist = std::distance(((cit == first.curr_iter) ? first.curr_pos : cit->second.first), last.curr_pos);
            // if dist is negative, then last.curr_pos is before cpos.  not reachable.
            // else add the distance to running total.
            return ((dist < 0) ? ::std::numeric_limits<difference_type>::lowest() : (n + dist));


            //            // iterate until both are in the same subcontainer, or first is at end.
            //            while (!it.at_end() && (it.curr_iter != last.curr_iter)) {
            //              n += ::std::distance(it.curr_pos, it.curr_iter->second.second);
            //              ++(it.curr_iter);
            //              if (it.curr_iter != it.max_iter) it.curr_pos = (it.curr_iter)->second.first;
            //            }
            //
            //            // both are at same place (including end)
            //            if (it == last) return n;
            //            // if first is at its end, and first != last, then last is not reachable.
            //            if (it.at_end) return ::std::numeric_limits<difference_type>::lowest();
            //
            //            // first is not at end, and first != last, then last is in same container as first.
            //            n += ::std::distance(it.curr_pos, (last.curr_iter == last.max_iter) ? last.curr_iter->second.first : last.curr_pos);
            //
            //            return n;
          }
      };



      /// compact the vector after disjoint entries are deleted.
      void compact() {
        if (vec.size() == 0 || map.size() == 0) {
          map.clear_no_resize();
          vec.clear();
          return;
        }

        // allocate a new vector
        subcontainer_type temp;
        temp.resize(vec.size());
        auto curr = temp.begin();
        auto prev = map.begin();
        // for each map entry, copy content over.  the range iterators in map are moved to temp
        auto max = map.end();
        for (auto it = map.begin(); it != max; ++it) {
          // copy over the content
          it->second.second = std::copy(it->second.first, it->second.second, curr);

          // update the iterators
          it->second.first = curr;
          curr = it->second.second;

          prev = it;
        }
        // erase everything after.  curr will be pointing to temp.end() now.
        temp.erase(curr, temp.end());

        // swap
        vec.swap(temp);

        // now last iter's second (was point to the last valid element, which then became the end of temp, and now needs to be come vec's end.
        prev->second.second = vec.end();
        // (all other iterators are not invalided.  just one that does not point to a real element (i.e. end iterator)

      }

      /// compact the vector after disjoint entries are deleted.  ASSUMPTION: all entries with same key are contiguous in memory, and map points to a subrange of each of these ranges.
      void inplace_compact() {
        if (vec.size() == 0 || map.size() == 0) {
            map.clear_no_resize();
            vec.clear();
            return;
          }

        Less less = Less();

        auto compacted_it = vec.begin();
        auto key = compacted_it->first;
        auto map_it = map.find(key);


        size_t i = 0;

        // go through all elements in vec, copy to beginning,
        auto max = vec.end();
        for (auto it = vec.begin(); it != max; ++i) {
          // get the current key
          key = it->first;
          // find the map entry
          map_it = map.find(key);


          // if in map, then compact
          if ((map_it != map.end()) && (dist(map_it) > 0)) {

              it = map_it->second.first;  // use original range to jump past stuff.
              std::advance(it, dist(map_it) - 1);
              it = std::adjacent_find(it, max, less);

              // in map.  we copy over the entries to keep.  advances the target it
              map_it->second.second = std::copy(map_it->second.first, map_it->second.second, compacted_it);
              // update the map
              map_it->second.first = compacted_it;
              compacted_it = map_it->second.second;

              //printf("iter %ld curr valid count = %ld, curr source count = %ld\n", i, std::distance(vec.begin(), compacted_it), std::distance(vec.begin(), it));

          } // else not in map.  don't update map, don't update current target it, but do update current it.

          // but advance it.  sorted, so != is same as <
          it = std::adjacent_find(it, max, less);
          // we are now at the last of the identical.  loop will increment to next.

          if (it != max) ++it;  // if nothing is found, max is returned, so need to check that.
        }
        // erase the extra
        vec.erase(compacted_it, vec.end());

        // the last map entry needs to have its second point to vec.end() now.  all other iterators are still valid.
        map_it->second.second = vec.end();

      }



      /// rehash to rebuild the hashmap index.
      void rebuild() {
        map.clear_no_resize();
        s = 0UL;

        if (vec.size() == 0) return;

        map.resize(map.size() + vec.size());
        Less less = Less();

        auto first = vec.begin();
        auto key = first->first;
        auto max = vec.end();
        for (auto it = vec.begin(); it != max;) {
          first = it;
          key = first->first;

          // find the last of the entries with same key
          it = std::adjacent_find(it, max, less);

          if (it != max) {
            // not last entry, so advance 1.
            ++it;
          }
          map.insert(::std::make_pair(std::move(key), std::move(std::make_pair(first, it) ) ) );
          s += std::distance(first, it);

        }

//        size_t bucket_max = 0;
//        for (size_t i = 0; i < map.bucket_count(); ++i) {
//          bucket_max = ::std::max(bucket_max, map.bucket_size(i));
//        }
//        printf("map size: %ld, map buckets %ld, max_bucket %ld, map loadfactor %f\n", map.size(), map.bucket_count(), bucket_max, map.load_factor());
      }




    public:
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<const Key, T>;
      using hasher                = Hash;
      using key_equal             = Equal;
      using allocator_type        = Allocator;
      using reference             = value_type&;
      using const_reference       = const value_type&;
      using pointer               = typename std::allocator_traits<Allocator>::pointer;
      using const_pointer         = typename std::allocator_traits<Allocator>::const_pointer;
      using iterator              = concat_iter<value_type>;
      using const_iterator        = concat_iter<const value_type>;
      using size_type             = typename subcontainer_type::size_type;
      using difference_type       = typename subcontainer_type::difference_type;


      //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
      densehash_vecmap(Key empty_key, Key deleted_key,
                   size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Less& less = Less(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                           map(bucket_count, hash, equal, alloc), s(0UL) {
        map.set_empty_key(empty_key);
        map.set_deleted_key(deleted_key);
      };

      template<class InputIt>
      densehash_vecmap(InputIt first, InputIt last,
                        Key empty_key, Key deleted_key,
                         size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Less& less = Less(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                         vec(first, last),
                         map(vec.size(), hash, equal, alloc), s(0UL) {
          map.set_empty_key(empty_key);
          map.set_deleted_key(deleted_key);

          std::sort(vec.begin(), vec.end(), less);

          this->rebuild();
      };

      virtual ~densehash_vecmap() {};



      iterator begin() {
        return iterator(map.begin(), map.end(), map.begin()->second.first, 0);
      }
      const_iterator begin() const {
        return cbegin();
      }
      const_iterator cbegin() const {
        return const_iterator(map.begin(), map.end(), map.begin()->second.first, 0);
      }



      iterator end() {
        return iterator(map.end());
      }
      const_iterator end() const {
        return cend();
      }
      const_iterator cend() const {
        return const_iterator(map.end());
      }

      bool empty() const {
        return map.empty();
      }

      size_type size() const {
        return s;
      }

      size_type unique_size() const {
        return map.size();
      }

      void clear() {
    	  vec.clear();
        map.clear_no_resize();
        s = 0UL;
      }

      void resize(size_t const n) {
        map.resize(n);
      }

      /// rehash for new count number of BUCKETS.  iterators are invalidated.
      void rehash(size_type count) {
        // only rehash if new bucket count is greater than old bucket count
        if (count > map.bucket_count())
          map.rehash(count);
      }

      /// bucket count.  same as underlying buckets
      size_type bucket_count() { return map.bucket_count(); }

      /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
      float max_load_factor() {
        return (map.size() == 0) ? map.max_load_factor() : map.max_load_factor() * (static_cast<float>(vec.size()) / static_cast<float>(map.size()));
      }



      // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
      template <class InputIt>
      void insert(InputIt first, InputIt last) {
          static_assert(::std::is_convertible<std::pair<Key, T>,
                        typename ::std::iterator_traits<InputIt>::value_type>::value,
                        "InputIt value type for insert cannot be converted to std::pair<Key, T> type");

          if (first == last) return;

          size_t prev_size = vec.size();

          // copy it in.  g++ 4.8.x has non-compliance to c++11 s.t. vec.insert returns void instead of iterator.
          // so compute the iterator directly.
          vec.reserve(prev_size + std::distance(first, last));
          vec.insert(vec.end(), first, last);

          Less less = Less();
          // sort the new part
          typename subcontainer_type::iterator middle = vec.begin() + prev_size;
          std::sort(middle, vec.end(), less);

          // merge with previous
          if (prev_size > 0) ::std::inplace_merge(vec.begin(), middle, vec.end(), less);

          // rebuild index
          this->rebuild();
      }

      /// inserting sorted range
      void insert(::std::vector<::std::pair<Key, T> > & input) {

        if (input.size() == 0) return;

        if (vec.empty()) {
          vec.swap(input);

          Less less = Less();
          // sort the new part
          std::sort(vec.begin(), vec.end(), less);

          this->rebuild();
        }
        else {
          this->insert(input.begin(), input.end());
        }
      }

      template <typename InputIt, typename Pred>
      size_t erase(InputIt first, InputIt last, Pred const & pred) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");

        if (first == last) return 0;

        size_t count = 0;
//        bool erased = false;

        // mark for erasure
        auto middle = map.begin()->second.first;
        for (; first != last; ++first) {
          auto iter = map.find(*(first));
          if (iter == map.end()) continue;

          middle = ::std::partition(iter->second.first, iter->second.second, pred);

          count += std::distance(iter->second.first, middle);
          iter->second.first = middle;

          if (dist(iter) == 0) map.erase(iter);
//          erased = true;
        }
//        if (erased) this->inplace_compact();
        s -= count;
        return count;
      }

      template <typename InputIt>
      size_t erase(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");


        if (first == last) return 0;

        size_t count = 0;

//        bool erased = false;

        // mark for erasure
        for (; first != last; ++first) {
          auto iter = map.find(*first);
          if (iter == map.end()) continue;

          count += dist(iter);
          //iter->second.first = iter->second.second;
          map.erase(iter);

//          erased = true;
        }

//        if (erased) this->inplace_compact();
        s -= count;
        return count;
      }

      template <typename Pred>
      size_t erase(Pred const & pred) {

        if (this->size() == 0) return 0;

        size_t before = this->size();

        auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);
        vec.erase(new_end, vec.end());

        this->rebuild();

        return before - this->size();
      }

      size_type count(Key const & key) const {
        if (map.find(key) == map.end()) return 0;
        else return dist(key);
      }



      void report() {
          BL_INFOF("vecmap bucket count: %lu\n", map.bucket_count());
          BL_INFOF("vecmap load factor: %f\n", map.load_factor());
          BL_INFOF("vecmap unique entries: %lu\n", map.size());
          BL_INFOF("vecmap total size: %lu\n", s);
      }


      size_type get_max_multiplicity() const {
        size_type max_multiplicity = 0;
        auto max = map.end();
        for (auto it = map.begin(); it != max; ++it) {
          max_multiplicity = ::std::max(max_multiplicity, dist(it));
        }
        return max_multiplicity;
      }

      size_type get_min_multiplicity() const {
        size_type min_multiplicity = ::std::numeric_limits<size_type>::max();
        auto max = map.end();
        size_type ss = 0;
        for (auto it = map.begin(); it != max; ++it) {
          ss = dist(it);
          if (ss > 0)
            min_multiplicity = ::std::min(min_multiplicity, ss);
        }
        return min_multiplicity;
      }

      double get_mean_multiplicity() const {
        return static_cast<double>(vec.size()) / double(map.size());
      }
      double get_stdev_multiplicity() const {
        double stdev_multiplicity = 0;
        auto max = map.end();
        double key_s;
        for (auto it = map.begin(); it != max; ++it) {
          key_s = dist(it);
          stdev_multiplicity += (key_s * key_s);
        }
        return stdev_multiplicity / double(map.size()) - get_mean_multiplicity();
      }


      ::std::pair<subiter_type, subiter_type> equal_range_value_only(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(subiter_type(), subiter_type());

        return ::std::make_pair(iter->second.first, iter->second.second);
      }
      ::std::pair<const_subiter_type, const_subiter_type> equal_range_value_only(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(const_subiter_type(), const_subiter_type());

        return ::std::make_pair(iter->second.first, iter->second.second);

      }


      ::std::pair<iterator, iterator> equal_range(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(iterator(map.end()), iterator(map.end()));

        return ::std::make_pair(iterator(iter, map.end(), iter->second.first, 0),
                                iterator(iter, map.end(), iter->second.second, dist(iter->second)));
      }
      ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(const_iterator(map.end()), const_iterator(map.end()));

        return ::std::make_pair(const_iterator(iter, map.end(), iter->second.first, 0),
                                const_iterator(iter, map.end(), iter->second.second, dist(iter->second)));

      }
      // NO bucket interfaces

  };



} // end namespace fsc.

namespace std {

  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::densehash_vecmap<Key, T, Hash, Equal, Allocator>::iterator first,
       typename ::fsc::densehash_vecmap<Key, T, Hash, Equal, Allocator>::iterator last, OutputIterator result) {

    // can last be reach from first?
    if ((last - first) <= 0) return result;

    // reachable.  so now walk.  do not need to do as much checking.
    auto out_iter = result;

    // now try to increment first until we get to last, or fail at getting to last.
    auto cit = first.curr_iter;
    auto cpos = first.curr_pos;

    // walk until either we are in same subcontainer.
    // since last is reachable from first, we don't need to check first's end_iter.
    while (cit != last.curr_iter) {
      out_iter = ::std::copy(cpos, cit->second.second, out_iter);
      ++cit;
      cpos = cit->second.first;
    }

    // now we're in same container.  take care the rest.  cpos was updated if cit was moved.
    out_iter = ::std::copy(cpos, last.curr_pos, out_iter);

    return out_iter;
  }

  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::densehash_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator first,
       typename ::fsc::densehash_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator last, OutputIterator result) {

    // can last be reach from first?
    if ((last - first) <= 0) return result;

    // reachable.  so now walk.  do not need to do as much checking.
    auto out_iter = result;

    // now try to increment first until we get to last, or fail at getting to last.
    auto cit = first.curr_iter;
    auto cpos = first.curr_pos;

    // walk until either we are in same subcontainer.
    // since last is reachable from first, we don't need to check first's end_iter.
    while (cit != last.curr_iter) {
      out_iter = ::std::copy(cpos, cit->second.second, out_iter);
      ++cit;
      cpos = cit->second.first;
    }

    // now we're in same container.  take care the rest.  cpos was updated if cit was moved.
    out_iter = ::std::copy(cpos, last.curr_pos, out_iter);


    return out_iter;
  }

}




#endif /* SRC_WIP_DENSEHASH_VECMAP_HPP_ */
