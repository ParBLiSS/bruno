/**
 * @file		object_pool.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   this file defines in-memory Object Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef OBJECTPOOL_HPP_
#define OBJECTPOOL_HPP_

#include <cassert>

#include <atomic>
#include <mutex>
#include <deque>
#include <unordered_set>
#include <stdexcept>

#include "utils/logging.h"
#include "concurrent/concurrent.hpp"
#include "concurrent/lockfree_queue.hpp"

// TODO: change to an object pool, not just a buffer pool.

namespace bliss
{
  namespace io
  {
    // TODO: move constructor and assignment operators between ObjectPools of different thread safeties.
    // TODO: use threadsafe_queue where appropriate.

    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            Each buffer is a block of preallocated memory that can be appended into.
     *
     *            The caller can acquire a new buffer from the pool and release an old buffer back into the pool, by Id.
     *            Released buffers are marked as empty and available.
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      life cycle:  pool acquire() by a single thread.
     *              application:  uses object, potentially by multiple threads
     *              pool release():  returns object to pool.  this class supports releasing an object by multiple threads
     *                only first one succeeds.
     *
     *
     *           note that the use of "in_use" set serves 2 purposes:  ensure there is no memory leak, and prevent multiple releases for the same object.
     *
     * @note      Each object should be acquired by a single thread.
     *            object may be used by multiple threads
     *            the object, when using this class, MAY be released by multiple threads.
     *              this is NOT the case for CommunicationLayer and SendMessageBuffer classes
     *                these 2 releases an object via a single thread
     *                CommLayer will also guarantee that all in-use buffers are released back to buffer.
     *
     *            It is important to address race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data while thread 2 releases it;
     *              thread 2 and 3 release the same object, thread 1 successfully acquires before thread 3 releases object).
     *
     *  @note:  for now, rely on default constructor of T.  in the future, custom allocator may be better
     *
     *  @tparam LockType    The thread safety property for the pool
     *  @tparam T           The object instance type.  object pool stores pointers to T instances.
     */
    template<bliss::concurrent::LockType LockType, class T >
    class ObjectPool
    {
      public:
        /**
         * @brief     Index Type used to reference the Objects in the ObjectPool.
         */
        using ObjectType = T;
        using ObjectPtrType = T*;   // shared pointer allows atomic operations as well as check for expired pointers
                                    // however, GCC does not support it.

        static const bliss::concurrent::LockType poolLT = LockType;


      protected:
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Object Ids.
         */
        bliss::concurrent::ThreadSafeQueue<ObjectPtrType>  available;
        std::unordered_set<ObjectPtrType>                  in_use;

        std::atomic<int64_t>  size_in_use;

        /**
         * @brief     mutex to control access.
         */
        mutable std::mutex mutex;
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;

        /**
         * NOT thread safe, so need to be wrapped in synchronized calls.
         */
        void clear_storage() {
          ObjectPtrType ptr;
          auto entry = available.tryPop();
          while (entry.first) {
             delete entry.second;
             entry = available.tryPop();
          }


          std::unique_lock<std::mutex> lock(mutex);
          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();
          lock.unlock();

          size_in_use.store(0, std::memory_order_release);
        }



//        /**
//         * @brief     Private move constructor with mutex lock.
//         * @param other   Source ObjectPool object to move
//         * @param l       Mutex lock on the Source ObjectPool object
//         */
//        ObjectPool(ObjectPool<LockType, T>&& other, const std::lock_guard<std::mutex>&) :
//          capacity(other.capacity)
//
//        {
//          other.capacity = 0;
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//        };


        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const ObjectPool<LockType, T>& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        ObjectPool<LockType, T>& operator=(const ObjectPool<LockType, T>& other) = delete;


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(_pool_capacity), in_use(), size_in_use(0)
          {};


//        /**
//         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
//         * @param other    source ObjectPool object to move.
//         */
//        explicit ObjectPool(ObjectPool<LockType, T>&& other) :
//          ObjectPool<LockType, T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
//
//        /**
//         * @brief     move assignment operator.
//         * @param other   source ObjectPool object to move from.
//         * @return        self, with member variables moved from other.
//         */
//        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) {
//          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
//                                        otherlock(other.mutex, std::defer_lock);
//          std::lock(mylock, otherlock);
//
//          capacity = other.capacity; other.capacity = 0;
//
//          clear_storage();
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//
//          return *this;
//        };
        explicit ObjectPool(ObjectPool<LockType, T>&& other) = delete;
        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) = delete;



        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          capacity = 0;

          clear_storage();
        };

        /**
         * @brief     Current size of the ObjectPool.  For debugging only.
         * @return  size, type IdType (aka int).
         */
        int64_t getAvailableCount() const  {

          return (isUnlimited() ? capacity : (capacity - size_in_use.load(std::memory_order_relaxed)));
        }


        /**
         * @brief     Current capacity of the ObjectPool
         * @return    capacity, type IdTyype (aka int)
         */
        const int64_t getCapacity() const {
          return capacity;
        }

        inline const bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }


        /**
         * @brief     Resets the entire ObjectPool: all Objects in pool are  marked as released and available.
         *
         * @note      This is not entirely thread safe.  the available set is cleared by a single thread, but other
         *            threads may be acquiring objects while they are being released.
         *
         *            It is envisioned that this function should be called from a single thread.
         *            however, care must be taken to ensure that all other threads have completed their work, else data loss is likely.
         */
        void reset() {
          std::lock_guard<std::mutex> lock(mutex);

          // move all from in_use to available
          for (auto iter = in_use.begin(); iter != in_use.end(); ++iter) {
            available.tryPush(*iter);
          }

          // TODO: clear the objects somehow?
          in_use.clear();
          size_in_use.store(0, std::memory_order_release);
        }


        /**
         * @brief     Get the next available Object by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the ObjectId if successful.
         */
        template<bliss::concurrent::LockType LT = LockType>
                typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, ObjectPtrType>::type tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          int64_t prev_size = size_in_use.fetch_add(1, std::memory_order_relaxed);
          if (prev_size >= capacity) {
            size_in_use.fetch_sub(1, std::memory_order_relaxed);
            // leave ptr as nullptr.
            if (this->isUnlimited()) {
              ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", prev_size);
            } else {
              WARNINGF("WARNING: pool is full. prev size %lu.", prev_size);
            }
          } else {

            // now get or create
            if (available.isEmpty()) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();

            } else {
              // has available for reuse.
              sptr = available.tryPop().second;
              // available may return null because it's concurrent queue.
              // cannot wait for available to return non-null, because available could be empty at this point.
            }

            if (sptr) {

              std::unique_lock<std::mutex> lock(mutex);

              // get the object ready for use.
              in_use.insert(sptr);  // store the shared pointer.

              lock.unlock();

            }
          }

          return sptr;

        }

        /**
         * @brief     Get the next available Object by id.  if none available, return false in the first argument of the pair.
         * @return    object pointer. if null, acquire was not successful.
         */
        template<bliss::concurrent::LockType LT = LockType>
                typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, ObjectPtrType>::type tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          int64_t prev_size = size_in_use.fetch_add(1, std::memory_order_relaxed);
          if (prev_size >= capacity) {

            size_in_use.fetch_sub(1, std::memory_order_relaxed);

            if (this->isUnlimited()) {
              ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", prev_size);
            } else {
              WARNINGF("WARNING: pool is full. prev size %lu.", prev_size);
            }

            // leave ptr as nullptr.
          } else {

            // now get or create
            if (available.isEmpty()) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();

            } else {
              // has available for reuse.
              sptr = available.tryPop().second;

              // available may return null because it's concurrent queue.
              // cannot wait for available to return non-null, because available could be empty at this point.
            }

            if (sptr) {
              while (spinlock.test_and_set());

              // get the object ready for use.
              in_use.insert(sptr);  // store the shared pointer.

              spinlock.clear();
            }
          }

          return sptr;
        }


        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD SAFE
         * @param ptr weak_ptr to object.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type releaseObject(ObjectPtrType ptr) {

          if (!ptr)
          {
            ERRORF("ERROR pool releasing a nullptr.");
            throw std::logic_error("ERROR pool releasing a nullptr.");
          }

          std::unique_lock<std::mutex> lock(mutex);
          int count = in_use.erase(ptr);
          lock.unlock();

          bool res = false;            // nullptr would not be in in_use.

          if (count > 0) { // object is in in_use.  if not in use, release could double count the entry.
            // now make object available.  make sure push_back is done one thread at a time.
            res = available.tryPush(ptr);
            size_in_use.fetch_sub(count, std::memory_order_release);
          } // else return false.

          return res;
        }


        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD SAFE
         * @param objectId    The id of the Object to be released.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, bool>::type releaseObject(ObjectPtrType ptr) {

          if (!ptr)
          {
            ERRORF("ERROR pool releasing a nullptr.");
            throw std::logic_error("ERROR pool releasing a nullptr.");
          }

          while (spinlock.test_and_set());
          int count = in_use.erase(ptr);
          spinlock.clear();

          bool res = false;            // nullptr would not be in in_use.

          if (count > 0) { // object is in in_use.  if not in use, release could double count the entry.
            // now make object available.  make sure push_back is done one thread at a time.
            res = available.tryPush(ptr);
            size_in_use.fetch_sub(count, std::memory_order_release);

          } // else return false.

          return res;

        }
    };



    template<class T >
    class ObjectPool<bliss::concurrent::LockType::NONE, T>
    {
      public:
        /**
         * @brief     Index Type used to reference the Objects in the ObjectPool.
         */
        using ObjectType = T;
        using ObjectPtrType = T*;   // shared pointer allows atomic operations as well as check for expired pointers
                                    // however, GCC does not support it.

        static const bliss::concurrent::LockType poolLT = bliss::concurrent::LockType::NONE;


      protected:


        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Object Ids.
         */
        std::deque<ObjectPtrType>                          available;
        std::unordered_set<ObjectPtrType>                  in_use;

        /**
         * NOT thread safe, so need to be wrapped in synchronized calls.
         */
        void clear_storage() {
          ObjectPtrType ptr;
          while (!available.empty()) {
           ptr = available.front();
           available.pop_front();
           delete ptr;
          }
          available.clear();

          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();

        }



//        /**
//         * @brief     Private move constructor with mutex lock.
//         * @param other   Source ObjectPool object to move
//         * @param l       Mutex lock on the Source ObjectPool object
//         */
//        ObjectPool(ObjectPool<lt, T>&& other, const std::lock_guard<std::mutex>&) :
//          capacity(other.capacity)
//
//        {
//          other.capacity = 0;
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//        };


        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const ObjectPool<poolLT, T>& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        ObjectPool<poolLT, T>& operator=(const ObjectPool<poolLT, T>& other) = delete;


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(), in_use()
          {};


//        /**
//         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
//         * @param other    source ObjectPool object to move.
//         */
//        explicit ObjectPool(ObjectPool<LockType, T>&& other) :
//          ObjectPool<LockType, T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
//
//        /**
//         * @brief     move assignment operator.
//         * @param other   source ObjectPool object to move from.
//         * @return        self, with member variables moved from other.
//         */
//        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) {
//          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
//                                        otherlock(other.mutex, std::defer_lock);
//          std::lock(mylock, otherlock);
//
//          capacity = other.capacity; other.capacity = 0;
//
//          clear_storage();
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//
//          return *this;
//        };
        explicit ObjectPool(ObjectPool<poolLT, T>&& other) = delete;
        ObjectPool<poolLT, T>& operator=(ObjectPool<poolLT, T>&& other) = delete;

        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          capacity = 0;

          clear_storage();
        };


        /**
         * @brief     Current size of the ObjectPool
         * @return  size, type IdType (aka int).
         */
        int64_t getAvailableCount() const  {
          return (isUnlimited() ? capacity : (capacity - in_use.size()));
        }

        /**
         * @brief     Current capacity of the ObjectPool
         * @return    capacity, type IdTyype (aka int)
         */
        const int64_t getCapacity() const {
          return capacity;
        }

        inline const bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }


        /**
         * @brief     Resets the entire ObjectPool: all Objects in pool are  marked as released and available.
         *
         * @note      This is not entirely thread safe.  the available set is cleared by a single thread, but other
         *            threads may be acquiring objects while they are being released.
         *
         *            It is envisioned that this function should be called from a single thread.
         *            however, care must be taken to ensure that all other threads have completed their work, else data loss is likely.
         */
        void reset() {

          // move all from in_use to available
          available.insert(available.end(), in_use.begin(), in_use.end());
          // TODO: clear the objects somehow?
          in_use.clear();

        }


        ObjectPtrType tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          size_t size = 0;

            // now get or create
            if (available.empty()) {
              if ((size = getAvailableCount()) > 0) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();
              }  else {
                // else already nullptr.
                if (this->isUnlimited()) {
                  ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", size);
                } else {
                  WARNINGF("WARNING: pool is full. prev size %lu.", size);
                }
              }
            } else {
              // has available for reuse.
              sptr = available.front();
              available.pop_front();


              // available may return null because it's concurrent queue.
              // cannot wait for available to return non-null, because available could be empty at this point.
            }

            // get the object ready for use.
            if (sptr) in_use.insert(sptr);  // store the shared pointer.
          return sptr;
        }

        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD UNSAFE
         * @param objectId    The id of the Object to be released.
         */
        bool releaseObject(ObjectPtrType ptr) {

          if (!ptr)
          {
            ERRORF("ERROR pool releasing a nullptr.");
            throw std::logic_error("ERROR pool releasing a nullptr.");
          }

          bool res = false;            // nullptr would not be in in_use.

          if (in_use.erase(ptr) > 0) { // object is in in_use.  if not in use, release could double count the entry.
            // now make object available.  make sure push_back is done one thread at a time.
            available.push_back(ptr);
            res = true;
          } // else return false.

          return res;
        }
    };

    template<bliss::concurrent::LockType LockType, class T>
    const bliss::concurrent::LockType ObjectPool<LockType, T>::poolLT;


  } /* namespace io */
} /* namespace bliss */

#endif /* OBJECTPOOL_HPP_ */
