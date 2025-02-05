#ifndef PARLAY_HASH_H_
#define PARLAY_HASH_H_

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <iterator>
#include <optional>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <utils/epoch.h>
#include "bigatomic.h"
#include "parallel.h"

constexpr bool PrintGrow = false;

namespace parlay {

template <typename Entries>
struct parlay_hash {
  using Entry = typename Entries::Entry;
  using K = typename Entry::Key;

  // *********************************************
  // Various parameters
  // *********************************************
  
  // set to grow by factor of 8 (2^3)
  static constexpr int log_grow_factor = 2;
  static constexpr int grow_factor = 1 << log_grow_factor;

  // groups of block_size buckets are copied over by a single thread
  // the block size typically grows with size, but starts here
  static constexpr long min_block_size = 4;

  // buffer_size is picked so state fits in a cache line (if it can)
  static constexpr long buffer_size = (sizeof(Entry) > 24) ? 1 : 48 / sizeof(Entry);

  // log_2 of the expected number of entries in a bucket (<= buffer_size)
  static constexpr long log_bucket_size = 
    (buffer_size == 1) ? 0 : ((buffer_size == 2) ? 1 : ((buffer_size <= 4) ? 2 : ((buffer_size <= 8) ? 3 : 3)));

  static long get_block_size(int num_bits) {
    return num_bits < 16 ? 16 : 256; }

  // The size of a bucket that causes the table to grow, i.e. if any
  // insert causes the bucket to reach the given size, then start
  // growing.
  // Technically this should be something like c log (n) / log(log n))
  // for a small constant c if each bucket is expected to hold 1
  // element, but.... each bucket can be expected to hold more than one.
  static long get_overflow_size(int num_bits) {
    if constexpr (log_bucket_size == 0) return num_bits < 18 ? 10 : 16;
    else if constexpr (log_bucket_size == 1) return num_bits < 18 ? 11 : 18;
    else if constexpr (log_bucket_size == 2) return num_bits < 18 ? 12 : 20;
    else if constexpr (log_bucket_size == 3) return num_bits < 18 ? 14 : 22;
    else return num_bits < 18 ? 20 : 24;
  }

  // clear_at_end will cause the scheduler and epoch-based collector
  // to clear their state on destruction
  static constexpr bool default_clear_at_end = true;
  bool clear_memory_and_scheduler_at_end;

  // a reference to the scheduler (null if not to be cleared)
  parlay::scheduler_type* sched_ref;

  // *********************************************
  // The state structure for each bucket
  // *********************************************
  
  // for overflow lists for each bucket
  struct link {
    Entry entry;
    link* next;
    link(const Entry& entry, link* next) : entry(entry), next(next) { }
  };

  // for delayed reclamation of links using an epoch-based collector
  epoch::memory_pool<link>* link_pool;

  link* new_link(const Entry& entry, link* l) {
    return link_pool->New(entry, l); }
  void retire_link(link* l) { link_pool->Retire(l);}

  // Each bucket contains a "state", which consists of a fixed size
  // buffer of entries (buffer_size) and an overflow list.  The first
  // buffer_size entries in the bucket are kept in the buffer, and any
  // overflow goes to the list.  The head stores both the pointer to
  // the overflow list (lower 56 bits) and the number of elements in
  // the buffer, or buffer_size+1 if overfull (top 8 bits).
  struct state {
  public:
    size_t list_head;
    Entry buffer[buffer_size];
    state() : list_head(0) {}
    state(const Entry& e) : list_head(1ul << 48) {
      buffer[0] = e;
    }
    static constexpr size_t forwarded_val = 1ul;
    
    size_t make_head(link* l, size_t bsize) {
      return (((size_t) l) | (bsize << 48)); }

    // update overflow list with new ptr (assumes buffer is full)
    state(const state& s, link* ptr)
      : list_head(make_head(ptr, buffer_size + (ptr != nullptr))) {
      for (int i=0; i < buffer_size; i++)
	buffer[i] = s.buffer[i];
    }

    // add entry to the bucket state (in buffer if fits, otherwise at head of overflow list)
    template <typename NL>
    state(const state& s, Entry e, const NL& new_link) {
      for (int i=0; i < std::min(s.buffer_cnt(), buffer_size); i++) 
	buffer[i] = s.buffer[i];
      if (s.buffer_cnt() < buffer_size) {
	buffer[s.buffer_cnt()] = e;
	list_head = make_head(nullptr, s.buffer_cnt() + 1);
      } else {
	link* l = new_link(e, s.overflow_list());
	list_head = make_head(l, buffer_size + 1);
      }
    }

    // add entry to buffer (assumes it fits) -- specialization of above
    state(const state& s, Entry e) : list_head(make_head(nullptr, s.buffer_cnt() + 1)) {
      for (int i=0; i < s.buffer_cnt(); i++) 
	buffer[i] = s.buffer[i];
      buffer[s.buffer_cnt()] = e;
    }

    // remove buffer entry j, replace with first from overflow list (assumes there is overflow)
    state(const state& s, link* ptr, int j)
      : list_head(make_head(ptr->next, buffer_size + (ptr->next != nullptr))) {
      for (int i=0; i < buffer_size; i++)
	buffer[i] = s.buffer[i];
      buffer[j] = Entry{ptr->entry};
    }

    // remove buffer entry j, replace with last entry in buffer (assumes no overflow)
    state(const state& s, int j) : list_head(make_head(nullptr, s.buffer_cnt() - 1)) {
      if (s.overflow_list() != nullptr) abort();
      for (int i=0; i < s.buffer_cnt(); i++)
	buffer[i] = s.buffer[i];
      buffer[j] = buffer[s.buffer_cnt() - 1];
    }

    state(bool x) : list_head(forwarded_val) {}
    
    bool is_forwarded() const {return list_head == forwarded_val ;}

    // number of entries in buffer, or buffer_size+1 if overflow
    long buffer_cnt() const {return (list_head >> 48) & 255ul ;}

    // number of entries in bucket (includes those in the overflow list)
    long size() const {
      if (buffer_cnt() <= buffer_size) return buffer_cnt();
      return buffer_size + list_length(overflow_list());
    }

    // get the overflow list
    link* overflow_list() const {
      return (link*) (list_head & ((1ul << 48) - 1));}
  };

  // returns std::optional(f(entry)) for entry with given key
  template <typename F>
  static auto find_in_list(const link* nxt, const K& k, const F& f) {
    using rtype = typename std::invoke_result<F,Entry>::type;
    long cnt = 0;
    while (nxt != nullptr && !nxt->entry.equal(k)) {
      nxt = nxt->next;
      cnt++;
    }
    if (nxt == nullptr)
      return std::pair(std::optional<rtype>(), cnt);
    else
      return std::pair(std::optional<rtype>(f(nxt->entry)), 0l);
  }

  // If k is found copies list elements up to k, and keeps the old
  // tail past k.  Returns the number of new nodes that will need to
  // be reclaimed, the head of the new list, and the link that is removed.
  // Returns [0, nullptr, nullptr] if k is not found
  std::tuple<int, link*, link*> remove_from_list(link* nxt, const K& k) {
    if (nxt == nullptr)
      return std::tuple(0, nullptr, nullptr);
    else if (nxt->entry.equal(k))
      return std::tuple(1, nxt->next, nxt);
    else {
      auto [len, ptr, removed] = remove_from_list(nxt->next, k);
      if (len == 0) return std::tuple(0, nullptr, nullptr);
      return std::tuple(len + 1, new_link(nxt->entry, ptr), removed);
    }
  }

  // update element with a given key in a list.  Uses path copying.
  // Returns a triple consisting of the position of the key in the list (1 based),
  // the head of the new list with the key updated, and the old link that is replaced.
  // If the key is not found, nothing is done, the last two results are nullptr, and
  // the first result is the length of the list.
  template <typename Constr>
  std::tuple<int, link*, link*> update_list(link* nxt, const K& k, const Constr& constr) {
    if (nxt == nullptr) 
      return std::tuple(0, nullptr, nullptr);
    else if (nxt->entry.equal(k))
      return std::tuple(1, link_pool->New(constr(std::optional(nxt->entry)), nxt->next), nxt);
    else {
      auto [len, ptr, updated] = update_list(nxt->next, k, constr);
      if (ptr == nullptr) return std::tuple(len + 1, nullptr, nullptr);
      return std::tuple(len + 1, link_pool->New(nxt->entry, ptr), updated);
    }
  }

  // retires first n elements of a list, but not the entries
  void retire_list_n(link* nxt, int n) {
    while (n > 0) {
      n--;
      link* tmp = nxt->next;
      retire_link(nxt);
      nxt = tmp;
    }
  }

  // Retires full list and their entries.  Used when destructing the
  // table.
  void retire_list_all(link* nxt) {
    while (nxt != nullptr) {
      link* tmp = nxt->next;
      entries_->retire_entry(nxt->entry);
      retire_link(nxt);
      nxt = tmp;
    }
  }

  // Retires full list, but not their entries. Used when copying to a
  // new list during expansion, i.e. the entries will be in the new
  // list and don't need to be retired.
  void retire_list(link* nxt) {
    while (nxt != nullptr) {
      link* tmp = nxt->next;
      retire_link(nxt);
      nxt = tmp;
    }
  }

  static long list_length(link* nxt) {
    long len = 0;
    while (nxt != nullptr) {
      len++;
      nxt = nxt->next;
    }
    return len;
  }

  // Find key if it is in the buffer. Return index.
  int find_in_buffer(const state& s, const K& k) {
    long len = s.buffer_cnt();
    for (long i = 0; i < std::min(len, buffer_size); i++)
      if (s.buffer[i].equal(k))
	return i;
    return -1;
  }

  // Apply f to all entries in the state.
  template <typename F>
  void for_each_in_state(const state& s, const F& f) {
    for (long i = 0; i < std::min(s.buffer_cnt(), buffer_size); i++)
      f(s.buffer[i]);
    link* l = s.overflow_list();
    while (l != nullptr) {
      f(l->entry);
      l = l->next;
    }
  }
    
  // Find entry with given key if in the bucket (state).  Return
  // optional of f applied to the entry if found, otherwise
  // std::nullopt.
  template <typename F>
  auto find_in_state(const state& s, const K& k, const F& f)
    -> std::optional<typename std::invoke_result<F,Entry>::type>
  {
    long len = s.buffer_cnt();
    for (long i = 0; i < std::min(len, buffer_size); i++)
      if (s.buffer[i].equal(k))
	return std::optional(f(s.buffer[i]));
    if (len <= buffer_size) return std::nullopt;
    return find_in_list(s.overflow_list(), k, f).first;
  }

  // A bucket is just an "atomic" state.
  // a big_atomic<x> is sort of like an std::atomic<x> but supports
  // load-linked, store-conditional, and is efficient when the x does
  // not fit in a machine word.
  using bckt = big_atomic<state>;

  // used for load-linked, store-conditionals
  using tag_type = typename big_atomic<state>::tag;

  // wrapper to ensure alignment
  struct alignas(64) bucket { bckt v; };

  // initialize an uninitialized bucket
  static void initialize(bucket& bck) {
    new (&bck.v) big_atomic<state>(state());
  }

  // *********************************************
  // The table structures
  // Each version increases in size, by grow_factor
  // *********************************************

  // status of a block of buckets, used when initializing and when copying to a new version
  enum status : char {Uninit, Initializing, Empty, Working, Done};

  // A single version of the table.
  // A version includes a sequence of "size" "buckets".
  // New versions are added as the hash table grows, and each holds a
  // pointer to the next larger version, if one exists.
  struct table_version {
    std::atomic<table_version*> next; // points to next version if created
    std::atomic<long> finished_block_count; //number of blocks finished copying
    long num_bits;  // log_2 of size
    size_t size; // number of buckets
    long block_size; // size of each block used for copying
    int overflow_size; // size of bucket to trigger next expansion
    bucket* buckets; // sequence of buckets
    //sequence<bucket> buckets; // sequence of buckets
    std::atomic<status>* block_status; // status of each block while copying

    // The index of a key is the highest num_bits of the lowest
    // 48-bits of the hash value.  Using the highest num_bits ensures
    // that when growing, a bucket will go to grow_factor contiguous
    // buckets in the next table.
    long get_index(const K& k) {
      size_t h = Entry::hash(k);
      return (h >> (48 - num_bits))  & (size-1u);}

    bckt* get_bucket(const K& k) {
      return &buckets[get_index(k)].v; }

    // initial table version, n indicating size
    table_version(long n) 
      : next(nullptr),
	finished_block_count(0),
	num_bits(std::max<long>((long) std::ceil(std::log2(min_block_size-1)),
				(long) std::ceil(std::log2(1.5*n)) - log_bucket_size)),
	size(1ul << num_bits),
	block_size(num_bits < 10 ? min_block_size : get_block_size(num_bits)),
	overflow_size(get_overflow_size(num_bits))
    {
      //if (PrintGrow) std::cout << "initial size: " << size << std::endl;
      buckets = (bucket*) malloc(sizeof(bucket)*size);
      block_status = (std::atomic<status>*) malloc(sizeof(std::atomic<status>) * size/block_size);
      parallel_for(0, size, [&] (long i) { initialize(buckets[i]);});
      parallel_for(0, size/block_size, [&] (long i) { block_status[i] = Empty;});
    }

    // expanded table version copied from smaller version t
    table_version(table_version* t)
      : next(nullptr),
	finished_block_count(0),
	num_bits(t->num_bits + log_grow_factor),
	size(t->size * grow_factor),
	block_size(get_block_size(num_bits)),
	overflow_size(get_overflow_size(num_bits))
    {
      buckets = (bucket*) malloc(sizeof(bucket)*size);
      block_status = (std::atomic<status>*) malloc(sizeof(std::atomic<status>) * size/block_size);
    }

    ~table_version() {
      free(buckets);
      free(block_status);
    }
  };

  // the current table version
  std::atomic<table_version*> current_table_version;

  // the initial table version, used for cleanup on destruction
  table_version* initial_table_version;

  // *********************************************
  // Functions for expanding the table
  // *********************************************

  // Called when table should be expanded (i.e. when some bucket is too large).
  // Allocates a new table version and links the old one to it.
  void expand_table(table_version* ht) {
    table_version* htt = current_table_version.load();
    if (htt->next == nullptr) {
      long n = ht->size;
      // if fail on lock, someone else is working on it, so skip
      get_locks().try_lock((long) ht, [&] {
	 if (ht->next == nullptr) {
	   ht->next = new table_version(ht);
	   //if (PrintGrow)
	   //  std::cout << "expand to: " << n * grow_factor << std::endl;
	 }
	 return true;});
    }
  }

  // Copies a bucket into grow_factor new buckets.
  void copy_bucket(table_version* t, table_version* next, long i) {
    long exp_start = i * grow_factor;
    // Clear grow_factor buckets in the next table version to put them in.
    for (int j = exp_start; j < exp_start + grow_factor; j++)
      initialize(next->buckets[j]); 
    // copy bucket to grow_factor new buckets in next table version
    while (true) {
      // the bucket to copy
      auto [s, tag] = t->buckets[i].v.ll();

      // insert into grow_factor buckets (states) for next larger table
      state hold[grow_factor];
      size_t mask = grow_factor-1;
      for_each_in_state(s, [&] (const Entry& entry) {
	size_t idx = next->get_index(entry.get_key()) & mask;
       	hold[idx] = state(hold[idx], entry,
			  [&] (const Entry& e, link* l) {return new_link(e,l);});
      });

      // now store the buckets into table
      for (int j = 0; j < grow_factor; j++)
	next->buckets[grow_factor * i + j].v.store_sequential(hold[j]);

      // try to replace original bucket with forwarded marker
      if (t->buckets[i].v.sc(tag, state(true))) {
	retire_list(s.overflow_list()); 
	break;
      }
      
      // If the attempt failed then someone updated bucket in the meantime so need to retry.
      // Before retrying need to clear out already added buckets.
      for (int j = exp_start; j < exp_start + grow_factor; j++) {
	state ss = next->buckets[j].v.load();
	retire_list(ss.overflow_list());
	next->buckets[j].v.store_sequential(state());
      }
    }
  }

  // If copying is ongoing (i.e., next is not null), and if the the
  // hash bucket given by hashid is not already copied, tries to copy
  // the block_size buckets that containing hashid to the next larger
  // table version.
  void copy_if_needed(table_version* t, long hashid) {
    table_version* next = t->next.load();
    if (next != nullptr) {
      long num_blocks = t->size/t->block_size;
      long block_num = hashid & (num_blocks -1);
      long start = block_num * t->block_size;
      status st = t->block_status[block_num];
      status old = Empty;
      if (st == Done) return;

      // if data is uninitialized, need to initialize
      // if (st == Uninit || st == Initializing) {
      // 	status x = Uninit;
      // 	if (t->block_status[block_num].compare_exchange_strong(x, Working)) {
      // 	  for (int i = start; i < start + t->block_size; i++) 
      // 	    initialize(t->buckets[i]);
      // 	  t->block_status[block_num] = Empty;
      // 	} else {
      // 	  while (t->block_status[block_num] == Initializing)
      // 	    for (volatile int i=0; i < 100; i++);
      // 	}
      // }
	
      // This is effectively a try lock on the block_num.
      // It blocks other updates on the buckets associated with the block.
      else if (st == Empty &&
	       t->block_status[block_num].compare_exchange_strong(old, Working)) {

	// initialize block_status for next grow round
	for (int i = 0; i < grow_factor; i++)
	  next->block_status[grow_factor*block_num + i] = Empty;
	
	// copy block_size buckets
	for (int i = start; i < start + t->block_size; i++) {
	  copy_bucket(t, next, i);
	}
	t->block_status[block_num] = Done;
	
	// If all blocks have been copied then can set current table
	// to next.  Note: this atomic fetch-and-add can be a
	// bottleneck and is the reason the block sizes are reasonably
	// large (e.g. 256).  A smarter combining tree could be used
	// if smaller block sizes are needed.
	if (++next->finished_block_count == num_blocks) {
	  //std::cout << "expand done" << std::endl;
	  current_table_version = next;
	}
      } else {
	// If another thread is working on the block, wait until Done
	while (t->block_status[block_num] == Working) {
	  for (volatile int i=0; i < 100; i++);
	}
      }
    }
  }
    
  // *********************************************
  // Construction and Destruction
  // *********************************************

  // Clear bucket, assuming it is not forwarded.
  void clear_bucket(bckt* b) {
    auto [s, tag] = b->ll();
    if (!s.is_forwarded() && b->sc(tag, state())) {
      for (int j=0; j < std::min(s.buffer_cnt(), buffer_size); j++) {
	entries_->retire_entry(s.buffer[j]);
      }
      retire_list_all(s.overflow_list());
    }
  }

  // Clears bucket or if the bucket is forwarded (during copying)
  // then clear the forwarded buckets.
  void clear_bucket_rec(table_version* t, long i) {
    bckt* b = &(t->buckets[i].v);
    state head = b->load();
    if (!head.is_forwarded())
      clear_bucket(b);
    else {
      table_version* next = t->next.load();
      for (int j = 0; j < grow_factor; j++)
	clear_bucket_rec(next, grow_factor * i + j);
    }
  }

  void clear_buckets() {
    table_version* ht = current_table_version.load();
    // clear buckets from current and future versions
    parallel_for(0, ht->size, [&] (size_t i) {
	clear_bucket_rec(ht, i);});
  }
  
  // Clear all memory.
  // Reinitialize to table of size 1 if specified, and by default.
  void clear(bool reinitialize = true) {
    clear_buckets();

    // now reclaim the arrays
    table_version* tv = initial_table_version;
    while (tv != nullptr) {
      table_version* tv_next = tv->next;
      delete tv;
      tv = tv_next;
    }
    // reinitialize
    if (reinitialize) {
      current_table_version = new table_version(1);
      initial_table_version = current_table_version;
    }
  }

  Entries* entries_;
  
  // Creates initial table version for the given size.  The
  // clear_at_end allows to free up the epoch-based collector's
  // memory, and the scheduler.
  parlay_hash(long n, Entries* entries, bool clear_at_end = default_clear_at_end)
    : entries_(entries),
      clear_memory_and_scheduler_at_end(clear_at_end),
      sched_ref(clear_at_end ?
		new parlay::scheduler_type(std::thread::hardware_concurrency()) :
		nullptr),
      link_pool(clear_at_end ?
		new epoch::memory_pool<link>() :
		&epoch::get_default_pool<link>()),
      current_table_version(new table_version(n)),
      initial_table_version(current_table_version.load())
  { }

  ~parlay_hash() {
    clear(false);
    if (clear_memory_and_scheduler_at_end) {
      delete sched_ref;
      delete link_pool;
    }
  }

  // *********************************************
  // Operations
  // *********************************************

  // Updates b, s, tag, and idx to the correct bucket, state, tag and
  // index if the the state s is forwarded.  Is called recursively,
  // but unlikely to go more than one level, and when not growing will
  // return immediately.
  void check_bucket_and_state(table_version* t, const K& k,
			      big_atomic<state>*& b, state& s, tag_type& tag, long& idx) {
    if (s.is_forwarded()) {
      table_version* nxt = t->next.load();
      idx = nxt->get_index(k);
      b = &(nxt->buckets[idx].v);
      std::tie(s, tag) = b->ll();
      check_bucket_and_state(nxt, k, b, s, tag, idx);
    }
  }

  // find in the bucket, or if forwarded (during copying) then follow
  // through to the next table, possibly reapeatedly, although
  // unlikely.
  template <typename F>
  auto find_in_bucket_rec(table_version* t, bckt* s, const K& k, const F& f)
    -> std::optional<typename std::invoke_result<F,Entry>::type>
  {
    state x = s->load();
    //if bucket is forwarded, go to next version
    if (x.is_forwarded()) {
      table_version* nxt = t->next.load();
      return find_in_bucket_rec(nxt, nxt->get_bucket(k), k, f);
    }
    return find_in_state(x, k, f);
  }

  // Finds the entry with the key
  // Returns an optional which is empty if the key is not in the table,
  // and contains f(e) otherwise, where e is the entry matching the key
  // NOTE: this is the most important function to opmitize for performance
  // Hence one hand inline and one prefetch (not used anywhere else in code).
  template <typename F>
  auto Find(const K& k, const F& f)
    -> std::optional<typename std::invoke_result<F,Entry>::type>
  {
    table_version* ht = current_table_version.load();
    long idx = ht->get_index(k);
    bckt* b = &(ht->buckets[idx].v);
    // if entries are direct, then safe to scan the buffer without epoch protection
    if constexpr (Entry::Direct) {
      auto [s, tag] = b->ll();
      if (s.is_forwarded()) 
	check_bucket_and_state(ht, k, b, s, tag, idx);
      for (long i = 0; i < std::min(s.buffer_cnt(), buffer_size); i++)
	if (s.buffer[i].equal(k))
	  return std::optional(f(s.buffer[i]));
      // if not found and not overfull, then done
      if (s.buffer_cnt() <= buffer_size) return std::nullopt;
      // otherwise need to search overflow, which requires protection
      return epoch::with_epoch([&, tag=tag, &s = s] {
        // if state has not changed, then just search list
	if (b->lv(tag)) return find_in_list(s.overflow_list(), k, f).first;
	return find_in_bucket_rec(ht, b, k, f);
      });
    } else { // if using indirection always use protection
      __builtin_prefetch(b); // allows read to be pipelined with epoch announcement
      return epoch::with_epoch([&] () -> std::optional<typename std::invoke_result<F,Entry>::type> {
	  return find_in_bucket_rec(ht, b, k, f);});


    }
  }

  // Inserts at key, and does nothing if key already in the table.
  // The constr function construct the entry to be inserted if needed.
  // Returns an optional, which is empty if sucessfully inserted or
  // contains f(e) if not, where e is the entry matching the key.
  template <typename Constr, typename F>
  auto Insert(const K& key, const Constr& constr, const F& f)
    -> std::optional<typename std::invoke_result<F,Entry>::type>
  {
    using rtype = std::optional<typename std::invoke_result<F,Entry>::type>;
    return epoch::with_epoch([&] () -> rtype {
			       auto [e, flag] = insert_(key, constr);
			       if (flag) return {};
			       return rtype(f(e));});
  }

  template <typename Constr>
  auto insert_(const K& key, const Constr& constr) -> std::pair<Entry, bool> {
    table_version* ht = current_table_version.load();
    long idx = ht->get_index(key);
    auto b = &(ht->buckets[idx].v);
    int delay = 200;
    while (true) {
      auto [s, tag] = b->ll();
      copy_if_needed(ht, idx);
      check_bucket_and_state(ht, key, b, s, tag, idx);
      long len = s.buffer_cnt();
      // if found in buffer then done
      for (long i = 0; i < std::min(len, buffer_size); i++)
	if (s.buffer[i].equal(key)) return std::pair(s.buffer[i], false);
      if (len < buffer_size) { // buffer has space, insert to end of buffer
	Entry new_e = constr();
	if (b->sc(tag, state(s, new_e))) return std::pair(new_e, true);
	entries_->retire_entry(new_e); // if failed need to ty again
      } else if (len == buffer_size) { // buffer full, insert new link
	Entry new_e = constr();
	link* new_head = new_link(new_e, nullptr);
	if (b->sc(tag, state(s, new_head))) 
	  return std::pair(new_e, true);
	entries_->retire_entry(new_head->entry); // if failed need to try again
	retire_link(new_head);
      } else { // buffer overfull, need to check if in list
	auto [x, list_len] = find_in_list(s.overflow_list(), key, identity);
	if (list_len + buffer_size > ht->overflow_size) expand_table(ht);
	if (x.has_value()) return std::pair(*x, false); // if in list, then done
	Entry new_e = constr();
	link* new_head = new_link(new_e, s.overflow_list());
	if (b->sc(tag, state(s, new_head))) // try to add to head of list
	  return std::pair(new_e, true);
	entries_->retire_entry(new_head->entry); // if failed need to ty again
	retire_link(new_head);
      }
      // delay before trying again, only marginally helps
      for (volatile int i=0; i < delay; i++);
      delay = std::min(2*delay, 5000); // 1000-10000 are about equally good
    }
  }

  template <typename Constr, typename G>
  auto Upsert(const K& key, const Constr& constr, G& g)
    -> std::optional<typename std::invoke_result<G,Entry>::type>
  {
    using rtype = std::optional<typename std::invoke_result<G,Entry>::type>;
    table_version* ht = current_table_version.load();
    long idx = ht->get_index(key);
    auto b = &(ht->buckets[idx].v);
    return epoch::with_epoch([&] () -> rtype {
      int delay = 200;
      while (true) {
	auto [s, tag] = b->ll();
	state out_s = s;
	copy_if_needed(ht, idx);
	check_bucket_and_state(ht, key, b, s, tag, idx);
	long len = s.buffer_cnt();
	bool cont = false;
	for (long i = 0; i < std::min(len, buffer_size); i++) {
	  if (s.buffer[i].equal(key)) {
	    Entry new_e = constr(std::optional(s.buffer[i]));
	    out_s.buffer[i] = new_e;
	    if (b->sc(tag, out_s)) return g(s.buffer[i]);
	    else {
	      entries_->retire_entry(new_e);
	      cont = true;
	      break;
	    }
	  }
	}
	if (cont) continue;
	if (len < buffer_size) { // buffer has space, insert to end of buffer
	  Entry new_e = constr(std::optional<Entry>());
	  if (b->sc(tag, state(s, new_e))) return std::nullopt;
	  entries_->retire_entry(new_e); // if failed need to ty again
	} else if (len == buffer_size) { // buffer just full, insert new link
	  link* new_head = new_link(constr(std::optional<Entry>()), nullptr);
	  if (b->sc(tag, state(s, new_head))) 
	    return std::nullopt;
	  entries_->retire_entry(new_head->entry); // if failed need to try again
	  retire_link(new_head);
	} else { // buffer overfull, need to check if in list
	  link* old_head = s.overflow_list();
	  auto [list_len, new_head, updated] = update_list(old_head, key, constr);
	  if (new_head != nullptr) {
	    if (b->sc(tag, state(s, new_head))) {// try to add to head of list
	      rtype r = std::optional(g(updated->entry));
	      retire_list_n(old_head, list_len); // retire old list
	      return r;
	    } else retire_list_n(new_head, list_len);
	  } else {
	    if (list_len + buffer_size > ht->overflow_size) expand_table(ht);
	    new_head = new_link(constr(std::optional<Entry>()), old_head);
	    if (b->sc(tag, state(s, new_head))) // try to add to head of list
	      return std::nullopt;
	    entries_->retire_entry(new_head->entry); // if failed need to ty again
	    retire_link(new_head);
	  }	    
	}
	// delay before trying again, only marginally helps
	for (volatile int i=0; i < delay; i++);
	delay = std::min(2*delay, 5000); // 1000-10000 are about equally good
      }
    });
  }

  // Removes entry with given key
  // Returns an optional which is empty if the key is not in the table,
  // and contains f(e) otherwise, where e is the entry that is removed.
  template <typename F>
  auto Remove(const K& key, const F& f)
    -> std::optional<typename std::invoke_result<F,Entry>::type>
  {
    using rtype = std::optional<typename std::invoke_result<F,Entry>::type>;
    table_version* ht = current_table_version.load();
    long idx = ht->get_index(key);
    auto b = &(ht->buckets[idx].v);
    // if entries are direct safe to scan the buffer without epoch protection
    if constexpr (Entry::Direct) {
      auto [s, tag] = b->ll();
      copy_if_needed(ht, idx);
      check_bucket_and_state(ht, key, b, s, tag, idx);
      if (s.buffer_cnt() <= buffer_size) {
	int i = find_in_buffer(s, key);
	if (i == -1) return std::nullopt;
	if (b->sc(tag, state(s, i))) {
	  rtype r = f(s.buffer[i]);
	  entries_->retire_entry(s.buffer[i]);
	  return r;
	} // if sc failed, will need to try again
      }
    }
    // if buffer overfull, or indirect, then need to protect
    return epoch::with_epoch([&] () -> rtype {
      int delay = 200;
      while (true) {
        auto [s, tag] = b->ll();
	copy_if_needed(ht, idx);
	check_bucket_and_state(ht, key, b, s, tag, idx);
	int i = find_in_buffer(s, key);
	if (i >= 0) { // found in buffer
	  if (s.buffer_cnt() > buffer_size) { // need to backfill from list
	    link* l = s.overflow_list();
	    if (b->sc(tag, state(s, l, i))) {
	      rtype r = f(s.buffer[i]);
	      entries_->retire_entry(s.buffer[i]);
	      retire_link(l);
	      return r;
	    } // if sc failed, will need to try again
	  } else { // buffer not overfull, can backfill within buffer
	    if (b->sc(tag, state(s, i))) {
	      rtype r = f(s.buffer[i]);
	      entries_->retire_entry(s.buffer[i]);
	      return r;
	    } // if sc failed, will need to try again
	  }
	} else { // not found in buffer
	  if (s.buffer_cnt() <= buffer_size) // if not overful, then done
	    return std::nullopt;
	  auto [cnt, new_list, removed] = remove_from_list(s.overflow_list(), key);
          if (cnt == 0) // if not found in list then done
	    return std::nullopt;
	  // if found, try to update with the new list that has the element removed
          if (b->sc(tag, state(s, new_list))) { 
	    rtype r = f(removed->entry); 
	    entries_->retire_entry(removed->entry);
            retire_list_n(s.overflow_list(), cnt); // retire old list
            return r;
          } // if sc failed, will need to try again
          retire_list_n(new_list, cnt - 1); // failed, retire new list
	}
	for (volatile int i=0; i < delay; i++);
	delay = std::min(2*delay, 5000); // 1000-10000 are about equally good
      }
    });
  }

  // Size of bucket, or if forwarded, then sum sizes of all forwarded
  // buckets, recursively.
  long bucket_size_rec(table_version* t, long i) {
    state head = t->buckets[i].v.load();
    if (!head.is_forwarded())
      return  head.size();
    else {
      long sum = 0;
      table_version* next = t->next.load();
      for (int j = 0; j < grow_factor; j++)
	sum += bucket_size_rec(next, grow_factor * i + j);
      return sum;
    }
  }

  long size() {
    table_version* ht = current_table_version.load();
    return epoch::with_epoch([&] {
       return parlay::tabulate_reduce(ht->size, [&] (size_t i) {
	   return bucket_size_rec(ht, i);});});
  }

  template <typename F>
  void for_each_bucket_rec(table_version* t, long i, const F& f) {
    state s = t->buckets[i].v.load();
    if (!s.is_forwarded())
      for_each_in_state(s, f);
    else {
      table_version* next = t->next.load();
      for (int j = 0; j < grow_factor; j++)
	for_each_bucket_rec(next, grow_factor * i + j, f);
    }
  }

  // Apply function f to all entries of the table.  Works while updates are going on, and guarantees that:
  //   any element whose delete linearizes before the invocation will not be included
  //   any element whose insert linearizes after the response will not be included
  //   any element that is present from invocation to response will be included
  // Elements that are inserted or deleted between the invocation and response might or might not appear.
  // template <typename F>
  // parlay::sequence<Entry> entries(const F& f) {
  //   table_version* ht = current_table_version.load();
  //   return epoch::with_epoch([&] {
  //     auto s = parlay::tabulate(ht->size, [&] (size_t i) {
  //       parlay::sequence<Entry> r;
  // 	for_each_in_bucket_rec(ht, i, [&] (const Entry& entry) {
  // 	  r.push_back(f(entry));});
  // 	return r;});
  //     return flatten(s);});
  // }

  // Applies f to all elments in table.
  // Same pseudo-linearizable guarantee as entries and size.
  template <typename F>
  void for_each(const F& f) {
    table_version* ht = current_table_version.load();
    return epoch::with_epoch([&] {
      for(long i = 0; i < ht->size; i++)
	for_each_bucket_rec(ht, i, f);});
  }

  // *********************************************
  // Iterator
  // *********************************************

  struct Iterator {
  public:
    using value_type        = typename Entries::Data;
    using iterator_category = std::forward_iterator_tag;
    using pointer           = value_type*;
    using reference         = value_type&;
    using difference_type   = long;

  private:
    std::vector<Entry> entries;
    Entry entry;
    table_version* t;
    int i;
    long bucket_num;
    bool single;
    bool end;
    void get_next_bucket() {
      while (entries.size() == 0 && ++bucket_num < t->size)
	bucket_entries(t, bucket_num, entries, [] (const Entry& e) {return e;});
      if (bucket_num == t->size) end = true;
    }

  public:
    Iterator(bool end) : i(0), bucket_num(-2l), single(false), end(true) {}
    Iterator(table_version* t) : t(t),
      i(0), bucket_num(-1l), single(false), end(false) {
      get_next_bucket();
    }
    Iterator(Entry entry) : entry(entry), single(true), end(false) {}
    Iterator& operator++() {
      if (single) end = true;
      else if (++i == entries.size()) {
	i = 0;
	entries.clear();
	get_next_bucket();
      }
      return *this;
    }
    Iterator& operator++(int) {
      Iterator tmp = *this;
      if (single) end = true;
      else if (++i == entries.size()) {
	i = 0;
	entries.clear();
	get_next_bucket();
      }
      return tmp;
    }
    template<bool D = Entry::Direct, std::enable_if_t<D, int> = 0>
    const value_type operator*() {
      if (single) return entry.get_entry();
      return entries[i].get_entry();}

    template<bool D = Entry::Direct, std::enable_if_t<!D, int> = 0>
    const value_type& operator*() { 
      if (single) return entry.get_entry();
      return entries[i].get_entry();}

    bool operator!=(const Iterator& iterator) {
      return !(end ? iterator.end : (bucket_num == iterator.bucket_num &&
				     i == iterator.i));
    }
    bool operator==(const Iterator& iterator) {
      return !(*this != iterator);}
  };

  Iterator begin() { return Iterator(current_table_version.load());}
  Iterator end() { return Iterator(true);}

  static constexpr auto identity = [] (const Entry& entry) {return entry;};
  static constexpr auto true_f = [] (const Entry& entry) {return true;};

  
  template <typename Constr>
  std::pair<Iterator,bool> insert(const K& key, const Constr& constr) {
    return epoch::with_epoch([&] {
      auto [e,flag] = insert_(key, constr);
      return std::pair(Iterator(e), flag);});
  }

  Iterator erase(Iterator pos) {
    Remove(*pos.first, true_f);
    return Iterator(true);
  }

  size_t erase(const K& key) {
    return Remove(key, true_f).has_value();
  }

  Iterator find(const K& k) {
    auto r = Find(k, identity);
    if (!r.has_value()) return Iterator(true);
    auto x = Iterator(*r);
    return x;
  }

};

  static constexpr bool default_clear_at_end = true;

  // conditionally rehash if type Hash::avalanching is not defined
  template<typename Hash, typename ignore = void>
  struct rehash {
    size_t operator()(size_t h) {
      size_t x = h * UINT64_C(0xbf58476d1ce4e5b9); // linear transform
      return (x ^ (x >> 31));  // non-linear transform
    }};

  template<typename Hash>
  struct rehash<Hash, typename Hash::is_avalanching> {
    size_t operator()(size_t i) {return i;}};

  // Definition where entries of the hash table are stored indirectly
  // through a pointer.  This means the entries themselves will never
  // move, but requires a level of indirection when accessing them.
  // Tags the high-bits of pointers with part of the hash function so
  // one can avoid the indirection if the tags do not match.
  // Currently used for all types that are not trivially copyable.
  template <typename EntryData>
  struct IndirectEntries {
    using DataS = EntryData;
    using Data = typename DataS::value_type;
    using Hash = typename DataS::Hash;
    using KeyEqual = typename DataS::KeyEqual;
      
    struct Entry {
      using K = typename DataS::K;
      using Key = std::pair<const K*,size_t>;
      static constexpr bool Direct = false;
      Data* ptr;
      static Data* tag_ptr(size_t hashv, Data* data) {
	return (Data*) (((hashv >> 48) << 48) | ((size_t) data));
      }
      Data* get_ptr() const {
	return (Data*) (((size_t) ptr) & ((1ul << 48) - 1)); }
      static unsigned long hash(const Key& k) {
	return k.second;}
      bool equal(const Key& k) const {
	return (((k.second >> 48) == (((size_t) ptr) >> 48)) &&
		KeyEqual{}(DataS::get_key(*get_ptr()), *k.first)); }
      Key get_key() const { return make_key(DataS::get_key(*get_ptr()));}
      Data& get_entry() const { return *get_ptr();}
      static Key make_key(const K& key) {
	return Key(&key, rehash<Hash>{}(Hash{}(key)));}
      Entry(Key k, Data* data) : ptr(tag_ptr(hash(k), data)) {}
      Entry() {}
    };

    bool clear_at_end;
    using Key = typename Entry::Key;

    // a memory pool for the entries
    epoch::memory_pool<Data>* data_pool;

    IndirectEntries(bool clear_at_end=false) 
      : clear_at_end(clear_at_end),
	data_pool(clear_at_end ?
		  new epoch::memory_pool<Data>() :
		  &epoch::get_default_pool<Data>()) {}
    ~IndirectEntries() {
      if (clear_at_end) { delete data_pool;}
    }

    // allocates memory for the entry
    Entry make_entry(const Key& k, const Data& data) {
      return Entry(k, data_pool->New(data)); }

    // retires the memory for the entry
    void retire_entry(Entry& e) {
      data_pool->Retire(e.get_ptr()); }
  };

  // Definition where entries of the hash table are stored directly.
  // This means the entries might be moved during updates, including
  // insersions, removals, and resizing.  Currently used for trivially
  // copyable types.
  template <typename EntryData>
  struct DirectEntries {
    using DataS = EntryData;
    using Data = typename DataS::value_type;
    using Hash = typename DataS::Hash;
    using KeyEqual = typename DataS::KeyEqual;
    using K = typename DataS::K;

    struct Entry {
      using K = typename DataS::K;
      using Key = K;
      static const bool Direct = true;
      Data data; 
      static unsigned long hash(const Key& k) {
	return rehash<Hash>{}(Hash{}(k));}
      bool equal(const Key& k) const { return KeyEqual{}(get_key(), k); }
      static Key make_key(const K& k) {return k;}
      const K& get_key() const {return DataS::get_key(data);}
      const Data& get_entry() const { return data;}
      Entry(const Data& data) : data(data) {}
      Entry() {}
    };

    DirectEntries(bool clear_at_end=false) {}
    Entry make_entry(const K& k, const Data& data) {
      return Entry(data); }

    // retiring is a noop since no memory has been allocated for entries
    void retire_entry(Entry& e) {}
  };

  // template <typename EntryData>
  // struct DirectEntriesX {
  //   using DataS = EntryData;
  //   using Data = typename DataS::value_type;
  //   using Hash = typename DataS::Hash;
  //   using KeyEqual = typename DataS::KeyEqual;
  //   using K = typename DataS::K;

  //   struct Entry {
  //     using K = typename DataS::K;
  //     using Key = K;
  //     static const bool Direct = true;
  //     std::array<long,1 + (sizeof(Data)-1)/8> data;
  //     static unsigned long hash(const Key& k) {
  // 	return rehash<Hash>{}(Hash{}(k));}
  //     bool equal(const Key& k) const { return KeyEqual{}(get_key(), k); }
  //     static Key make_key(const K& k) {return k;}
  //     const K& get_key() const { return DataS::get_key(*((Data*) &data));}
  //     const Data& get_entry() const { return *((Data*) &data);}
  //     Entry(const Data& d) { new (&data) Data(d); }
  //     Entry() {}
  //   };

  //   bool clear_at_end;

  //   // a memory pool for the entries
  //   epoch::retire_pool<Data>* data_pool;

  //   DirectEntriesX(bool clear_at_end=false) 
  //     : clear_at_end(clear_at_end),
  // 	data_pool(clear_at_end ?
  //   		  new epoch::retire_pool<Data>() :
  //   		  &epoch::get_default_retire_pool<Data>())
  //   {}
  //   ~DirectEntriesX() {
  //     if (clear_at_end) { delete data_pool;}
  //   }

  //   // allocates memory for the entry
  //   Entry make_entry(const K& k, const Data& data) {
  //     return Entry(data);}

  //   // retires the memory for the entry
  //   void retire_entry(Entry& e) {
  //     data_pool->Retire((Data*) &(e.data)); 
  //   }
  // };
  

}  // namespace parlay
#endif  // PARLAY_HASH_H_
