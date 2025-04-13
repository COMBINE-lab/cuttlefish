/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#ifndef H_DATA_POOL
#define H_DATA_POOL

#include "Globals.h"

#include <vector>
#include <iostream>
#include <algorithm>

#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
namespace th = boost;
#else
#include <thread>
#include <mutex>
#include <condition_variable>

namespace th = std;
#endif

namespace rabbit {

namespace core {

/**
 * @brief DataPool class 
 * This class provide an data pool for reusing memory space
 */	
template <class _TDataType>
class TDataPool {
  typedef _TDataType DataType;
  typedef std::vector<DataType *> part_pool;

  const uint32 maxPartNum;
  const uint32 bufferPartSize;

  part_pool availablePartsPool;
  part_pool allocatedPartsPool;

  th::mutex mutex;
  th::condition_variable partsAvailableCondition;

 public:
  volatile uint32 partNum;
  static const uint32 DefaultMaxPartNum = 32;
  static const uint32 DefaultBufferPartSize = 1 << 23;

	/**
	 * @brief Constructor
	 * @param maxPartNum_ the maximum number of part contained in DataPool.
	 * @param bufferPartsize_ Bytes of each part in DataPool (eg. 1<<22 means 4MB each part)
	 */
  TDataPool(uint32 maxPartNum_ = DefaultMaxPartNum, uint32 bufferPartSize_ = DefaultBufferPartSize)
      : maxPartNum(maxPartNum_), bufferPartSize(bufferPartSize_), partNum(0) {
    if (bufferPartSize_ < DefaultBufferPartSize) std::cerr << "[warning]: your chunk buffer size maybe too small to hold one complete sequence!" << std::endl;
    availablePartsPool.resize(maxPartNum);
    allocatedPartsPool.reserve(maxPartNum);
  }

  ~TDataPool() {
    for (typename part_pool::iterator i = allocatedPartsPool.begin(); i != allocatedPartsPool.end(); ++i) {
      ASSERT(*i != NULL);
      delete *i;
    }
  }
	/**
	 * @brief Acquire an DataType data in DataPool and assign to part_	 
	 * @details Acquired data from availablePartsPool, if there is no available data in availablePartsPool, program will 
	 *       wait.
	 * @param part_ the pointer to acquireed space 
	 */
  void Acquire(DataType *&part_) {
    th::unique_lock<th::mutex> lock(mutex);

    while (partNum >= maxPartNum) partsAvailableCondition.wait(lock);

    ASSERT(availablePartsPool.size() > 0);

    DataType *&pp = availablePartsPool.back();
    availablePartsPool.pop_back();
    if (pp == NULL) {
      pp = new DataType(bufferPartSize);
      allocatedPartsPool.push_back(pp);
    } else {
      pp->Reset();
    }

    partNum++;
    part_ = pp;
  }
	/**
	 * @brief Realease data to DataPool
	 * @details Realease the data in part_ to AvailablePartsPool and notify
	 * @param part_ the pointer to be realesed to DataPool 
	 */
  void Release(const DataType *part_) {
    th::lock_guard<th::mutex> lock(mutex);

    ASSERT(part_ != NULL);
    ASSERT(partNum != 0 && partNum <= maxPartNum);
    ASSERT(std::find(allocatedPartsPool.begin(), allocatedPartsPool.end(), part_) != allocatedPartsPool.end());
    availablePartsPool.push_back((DataType *)part_);
    partNum--;

    partsAvailableCondition.notify_one();
  }
};

}  // namespace core

}  // namespace rabbit

#endif  // H_DATA_POOL
