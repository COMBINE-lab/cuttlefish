/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#ifndef H_DATAQUEUE
#define H_DATAQUEUE

#include <queue>

#include "Globals.h"

#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
namespace th = boost;
#else
#include <mutex>
#include <condition_variable>
namespace th = std;
#endif

namespace rabbit {

namespace core {

/*
 * @brief DataQueue class 
 * This class provide an data queue for synchronize data access
 */	
template <class _TDataType>
class TDataQueue {
  typedef _TDataType DataType;
  typedef std::queue<std::pair<int64, DataType *> > part_queue; 

  const uint32 threadNum;
  const uint32 maxPartNum;
  uint64 completedThreadMask;
  uint32 partNum;
  uint64 currentThreadMask;
  part_queue parts;

  th::mutex mutex;
  th::condition_variable queueFullCondition;
  th::condition_variable queueEmptyCondition;

 public:
  static const uint32 DefaultMaxPartNum = 64; //default max part number 
  static const uint32 DefaultMaxThreadtNum = 64; // default max thread number

	/*
	 * @brief Constructor 
	 * @param maxPartNum_ maximum parts in queue
	 * @param threadNum_ the number of threads push data into the queue
	 */
  TDataQueue(uint32 maxPartNum_ = DefaultMaxPartNum, uint32 threadNum_ = 1)
      : threadNum(threadNum_), maxPartNum(maxPartNum_), partNum(0), currentThreadMask(0) {
    ASSERT(maxPartNum_ > 0);
    ASSERT(threadNum_ >= 1);
    ASSERT(threadNum_ < 64);

    completedThreadMask = ((uint64)1 << threadNum) - 1;
  }

  ~TDataQueue() {}
  
	/*
	 * @brief return if the queue is empty 
	 */
  bool IsEmpty() { return parts.empty(); }

	/*
	 * @brief return if all task is complete
	 */
  bool IsCompleted() { return parts.empty() && currentThreadMask == completedThreadMask; }
	
  void SetCompleted() {
    th::lock_guard<th::mutex> lock(mutex);

    ASSERT(currentThreadMask != completedThreadMask);
    currentThreadMask = (currentThreadMask << 1) | 1;

    queueEmptyCondition.notify_all();
  }

  void Push(int64 partId_, const DataType *part_) {
    th::unique_lock<th::mutex> lock(mutex);

    while (partNum > maxPartNum) queueFullCondition.wait(lock);

    parts.push(std::make_pair(partId_, (DataType *)part_));
    partNum++;

    queueEmptyCondition.notify_one();
  }

  bool Pop(int64 &partId_, DataType *&part_) {
    th::unique_lock<th::mutex> lock(mutex);

    while ((parts.size() == 0) && currentThreadMask != completedThreadMask) queueEmptyCondition.wait(lock);

    if (parts.size() != 0) {
      partId_ = parts.front().first;
      part_ = parts.front().second;
      partNum--;
      parts.pop();
      queueFullCondition.notify_one();
      return true;
    }

    // assure this is impossible
    ASSERT(currentThreadMask == completedThreadMask);
    ASSERT(parts.size() == 0);
    return false;
  }

  bool TryPop(int64 &partId_, DataType *&part_) {
    th::unique_lock<th::mutex> lock(mutex);

    if (parts.size() == 0) {
      return false;
    } else {
      partId_ = parts.front().first;
      part_ = parts.front().second;
      partNum--;
      parts.pop();
      queueFullCondition.notify_one();
      return true;
    }

    // assure this is impossible
    ASSERT(currentThreadMask == completedThreadMask);
    ASSERT(parts.size() == 0);
    return false;
  }



  void Reset() {
    ASSERT(currentThreadMask == completedThreadMask);

    partNum = 0;
    currentThreadMask = 0;
  }
};

}  // namespace core

}  // namespace rabbit

#endif  // DATA_QUEUE_H
