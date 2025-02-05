
#ifndef KEY_VALUE_COLLATOR_HPP
#define KEY_VALUE_COLLATOR_HPP



#include "Spin_Lock.hpp"
#include "Key_Value_Iterator.hpp"

#include <sys/types.h>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <atomic>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <sys/stat.h>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <thread>
#include <cassert>


// =============================================================================

namespace key_value_collator
{


template <typename T_buf_> class Buffer_Pool;

// A class to: collate a collection of key-value pairs, deposited from multiple
// producers; and to iterate over the collated key-value collection. Keys are of
// type `T_key_`, values are of type `T_val_`, and the keys are hashed to their
// corresponding partitions with `operator()(T_key_ key)` of class `T_hasher_`.
template <typename T_key_, typename T_val_, typename T_hasher_>
class Key_Value_Collator
{
public:

    typedef std::pair<T_key_, T_val_> key_val_pair_t;
    typedef std::vector<key_val_pair_t> buf_t;  // Type of the data buffers.

    typedef Key_Value_Iterator<T_key_, T_val_> iter_t;  // Type of the collation iterator.


private:

    const T_hasher_ hash;   // Hasher object to hash the keys to a numerical address-space.

    const std::string work_file_pref;   // Path to the temporary working files used by the collator.

    static constexpr char work_file_pref_default[] = ".";   // Default value for the temporary working files' prefixes.
    static constexpr char partition_file_ext[] = ".part";   // File extensions of the temporary partition files.

    static constexpr std::size_t partition_count = (1 << 7);    // Number of partitions for the keys.
    static constexpr std::size_t partition_buf_mem = (1LU * 1024 * 1024);   // Maximum memory for a partition buffer: 1MB.
    static constexpr std::size_t partition_buf_elem_th = partition_buf_mem / sizeof(key_val_pair_t);    // Maximum number of pairs to keep in a partition buffer.

    std::vector<buf_t> partition_buf;   // `partition_buf[i]` is the in-memory buffer for partition `i`.
    std::vector<std::ofstream> partition_file;  // `partition_file[i]` is the disk-storage file for partition `i`.

    Buffer_Pool<buf_t*> buf_pool;   // Managed buffer collection to copy-in and process incoming data from the producers.
    const std::size_t buf_count;    // Number of concurrent buffers for the producers.
    static constexpr std::size_t buf_count_default = 16;    // Default value for the concurrent buffer count.

    std::thread* mapper;    // The background thread mapping key-value pairs to corresponding partitions.
    std::atomic<bool> stream_incoming;  // Flag denoting whether the incoming key-value streams have ended or not.


    // Returns the disk-file path for the partition `p_id`.
    const std::string partition_file_path(std::size_t p_id) const;

    // Maps the key-value pairs from the producers to the partitions
    // corresponding to the keys.
    void map();

    // Maps the key-value pairs from the data buffer `buf` to the partitions
    // corresponding to the keys.
    void map_buffer(const buf_t& buf);

    // Returns the corresponding partition ID for the key `key`.
    std::size_t get_partition_id(const T_key_& key) const;

    // Flushes the buffer of the partition with ID `p_id` to disk and
    // clears the buffer.
    void flush(std::size_t p_id);


public:

    // Key_Value_Collator() = delete;
    Key_Value_Collator(const Key_Value_Collator&) = delete;
    Key_Value_Collator& operator=(const Key_Value_Collator&) = delete;

    // Constructs a key-value pair collection object to collate (and iterate
    // over) key-value pairs produced from multiple producers. Temporary disk-
    // files used throughout the process will be stored at the path-prefix
    // `work_file_pref`. `buf_count` concurrent buffers would be used to store
    // and process the deposited data. It should be set to at least the number
    // of producers to avoid throttling of the producers; and a good heuristic
    // choice for this is twice the number of producers.
    Key_Value_Collator(const std::string& work_file_pref = work_file_pref_default, std::size_t buf_count = buf_count_default);

    ~Key_Value_Collator();

    class Aggregate_Result;
    mutable Aggregate_Result agg_result;

    // Returns an available free buffer.
    buf_t& get_buffer();

    // Returns the buffer `buf` to the collator with deposited data.
    void return_buffer(buf_t& buf);

    // Closes the deposit stream incoming from the producers and flushes the
    // remaining in-memory content to disk. All deposit operations from the
    // producers must be made before invoking this.
    void close_deposit_stream();

    // Collates the deposited key-value pairs, using at most `thread_count`
    // processor-threads. Also generates an aggregate result of the keys
    // iff `aggregate = true`.
    void collate(uint32_t thread_count, bool aggregate = false) const;

    // Returns an iterator pointing at the beginning of the collection.
    iter_t begin() const;

    // Returns an iterator pointing at the end of the collection.
    iter_t end() const;

    std::size_t unique_key_count() const { return agg_result.unique_key_count; }    // Returns the number of unique keys.
    std::size_t pair_count() const { return agg_result.pair_count; }    // Returns the total number of key-value pairs.
    std::size_t mode_frequency() const { return agg_result.mode_count; }    // Returns the number of pairs with a most frequent key.
};


// A collection of objects of type `T_obj_`, ready to be used by multiple
// threads.
template <typename T_obj_>
class Object_Pool
{
private:

    std::vector<T_obj_> pool;   // The object collection.
    Spin_Lock lock_;    // Mutual-exclusion lock to avoid concurrent updates of the pool.
    std::atomic<std::size_t> size_; // Number of elements in the pool.


public:

    // Constructs an empty object-pool.
    Object_Pool(): size_(0)
    {}


    // Adds the object `obj` to the pool.
    void push(const T_obj_& obj)  // Should use move operand?
    {
        lock_.lock();
        pool.emplace_back(obj);
        size_.fetch_add(1, std::memory_order_acquire);  // size_++;
        lock_.unlock();
    }


    // Returns `true` iff the pool is empty.
    bool empty() const { return size_ == 0; }


    // Returns the size of the pool.
    std::size_t size() const { return size_; }


    // Tries to fetch an object to `obj` from the pool. Returns `true` iff
    // such an object is found.
    bool fetch(T_obj_& obj)
    {
        if(empty())
            return false;

        bool success = false;

        lock_.lock();

        if(!pool.empty())
        {
            size_.fetch_sub(1, std::memory_order_acquire);  // size_--;
            obj = pool.back();   // std::move(pool.back());
            pool.pop_back();

            success = true;
        }

        lock_.unlock();

        return success;
    }
};


// A managed collection of buffers of type `T_buf_`, ready to be used by
// multiple threads. Each buffer is either "free" or "full".
template <typename T_buf_>
class Buffer_Pool
{
private:

    Object_Pool<T_buf_> free_pool;  // Buffers available to be used.
    Object_Pool<T_buf_> full_pool;  // Buffers full of data.


public:

    Buffer_Pool() {}


    // Returns the number of available free buffers.
    std::size_t free_buf_count() const { return free_pool.size(); }


    // Returns the number of available data-full buffers.
    std::size_t full_buf_count() const { return full_pool.size(); }


    // Adds the buffer `buf` to the pool.
    void add_buf(const T_buf_& buf) { free_pool.push(buf); }


    // Tries to fetch a free buffer to `buf` from the pool. Returns `true` iff
    // such a buffer is found.
    bool fetch_free_buf(T_buf_& buf) { return free_pool.fetch(buf); }


    // Returns the buffer `buf` to the pool, presumably filled with data to be
    // used later.
    void return_full_buffer(const T_buf_& buf) { full_pool.push(buf); }


    // Tries to fetch a data-full buffer to `buf` from the pool. Returns `true`
    // iff such a buffer is found.
    bool fetch_full_buf(T_buf_& buf) { return full_pool.fetch(buf); }


    // Returns the buffer `buf` to the pool, to be reused later.
    void return_free_buf(const T_buf_& buf) { free_pool.push(buf); }
};


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Key_Value_Collator<T_key_, T_val_, T_hasher_>::Key_Value_Collator(const std::string& work_file_pref, const std::size_t buf_count):
    hash(),
    work_file_pref(work_file_pref),
    partition_buf(partition_count),
    partition_file(partition_count),
    buf_count(buf_count),
    mapper(nullptr),
    stream_incoming(true)
{
    static_assert(partition_buf_elem_th > 0, "Invalid configuration for partition buffer memory.");

    for(std::size_t p_id = 0; p_id < partition_count; ++p_id)
    {
        partition_buf[p_id].reserve(partition_buf_elem_th);
        partition_file[p_id].open(partition_file_path(p_id), std::ios::out | std::ios::binary);
    }

    for(std::size_t i = 0; i < buf_count; ++i)
        buf_pool.add_buf(new buf_t());

    mapper = new std::thread(&Key_Value_Collator<T_key_, T_val_, T_hasher_>::map, this);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Key_Value_Collator<T_key_, T_val_, T_hasher_>::~Key_Value_Collator()
{
    if(buf_pool.full_buf_count() > 0 || buf_pool.free_buf_count() != buf_count || mapper->joinable())
    {
        std::cerr << "Collator destructed while unprocessed buffers remained. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Release memory of the partition-buffers.
    buf_t* buf_p = nullptr;
    while(buf_pool.free_buf_count() > 0)
    {
        buf_pool.fetch_free_buf(buf_p);
        delete buf_p;
    }

    delete mapper;


    // Remove the partition-files.
    for(std::size_t p_id = 0; p_id < partition_count; ++p_id)
        if(std::remove(partition_file_path(p_id).c_str()))
        {
            std::cerr << "Error removing temporary files. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline typename Key_Value_Collator<T_key_, T_val_, T_hasher_>::buf_t& Key_Value_Collator<T_key_, T_val_, T_hasher_>::get_buffer()
{
    buf_t* buf_p;
    while(!buf_pool.fetch_free_buf(buf_p));

    return *buf_p;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Key_Value_Collator<T_key_, T_val_, T_hasher_>::return_buffer(buf_t& buf)
{
    buf_pool.return_full_buffer(&buf);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const std::string Key_Value_Collator<T_key_, T_val_, T_hasher_>::partition_file_path(const std::size_t p_id) const
{
    return work_file_pref + "." + std::to_string(p_id) + partition_file_ext;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Key_Value_Collator<T_key_, T_val_, T_hasher_>::map()
{
    buf_t* buf_p;

    while(stream_incoming || buf_pool.full_buf_count() > 0)
        if(buf_pool.fetch_full_buf(buf_p))
        {
            map_buffer(*buf_p);

            buf_p->clear();
            buf_pool.return_free_buf(buf_p);
        }
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Key_Value_Collator<T_key_, T_val_, T_hasher_>::map_buffer(const std::vector<Key_Value_Collator::key_val_pair_t>& buf)
{
    for(const auto& key_val_pair : buf)
    {
        const std::size_t p_id = get_partition_id(key_val_pair.first);
        auto& p_buf = partition_buf[p_id];
        p_buf.emplace_back(key_val_pair);

        assert(p_buf.size() <= partition_buf_elem_th);
        if(p_buf.size() == partition_buf_elem_th)
            flush(p_id);
    }
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Key_Value_Collator<T_key_, T_val_, T_hasher_>::flush(const std::size_t p_id)
{
    auto& buf  = partition_buf[p_id];
    auto& file = partition_file[p_id];
    if(!file.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(key_val_pair_t)))
    {
        std::cerr << "Error writing to partition file(s) of the collator. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    buf.clear();
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline std::size_t Key_Value_Collator<T_key_, T_val_, T_hasher_>::get_partition_id(const T_key_& key) const
{
    return hash(key) & (partition_count - 1);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Key_Value_Collator<T_key_, T_val_, T_hasher_>::close_deposit_stream()
{
    stream_incoming = false;
    if(!mapper->joinable())
    {
        std::cerr << "Early termination encountered for the key-mapper thread. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    mapper->join();


    // Flush the remaining in-memory partition contents, release their memory, and close the in-disk partitions.
    for(std::size_t p_id = 0; p_id < partition_count; ++p_id)
    {
        if(!partition_buf[p_id].empty())
            flush(p_id);

        buf_t().swap(partition_buf[p_id]);

        partition_file[p_id].close();
    }
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Key_Value_Collator<T_key_, T_val_, T_hasher_>::collate(const uint32_t thread_count, const bool aggregate) const
{
    std::vector<std::thread> worker;
    std::vector<Aggregate_Result> worker_aggregate(thread_count, Aggregate_Result());
    for(uint32_t t_id = 0; t_id < thread_count; ++t_id)
        worker.emplace_back(
            [this, aggregate](const uint32_t init_id, const uint32_t stride, Aggregate_Result& result)
            // Collates each partition with IDs starting from `init_id` and at stride lengths `stride`.
            {
                Aggregate_Result result_local;  // To avoid possible false-sharing.

                const auto file_size = [](const char* const file_name) -> off_t
                    {
                        struct stat st;
                        return stat(file_name, &st) == 0 ? st.st_size : 0;
                    };

                // Get the maximum partition size for this thread.
                off_t buf_sz = 0; // in bytes
                for(uint32_t p_id = init_id; p_id < partition_count; p_id += stride)
                    buf_sz = std::max(buf_sz, file_size(partition_file_path(p_id).c_str()));

                key_val_pair_t* const p_data = static_cast<key_val_pair_t*>(std::malloc(buf_sz));
                for(uint32_t p_id = init_id; p_id < partition_count; p_id += stride)
                {
                    // Read in the partition data to memory.

                    const std::string p_path = partition_file_path(p_id);
                    const auto p_bytes = file_size(p_path.c_str());

                    std::ifstream input(p_path.c_str(), std::ios::in | std::ios::binary);
                    input.read(reinterpret_cast<char*>(p_data), p_bytes);
                    assert(input.gcount() == p_bytes);
                    if(!input)
                    {
                        std::cerr << "Error reading the partition files. Aborting.\n";
                        std::exit(EXIT_FAILURE);
                    }

                    input.close();


                    // Sort the partition data and optionally get aggregate statistics.
                    const std::size_t elem_count = p_bytes / sizeof(key_val_pair_t);
                    std::sort(p_data, p_data + elem_count);

                    if(aggregate && elem_count > 0) // Aggregate results from this partition.
                    {
                        T_key_ curr_key = p_data[0].first;
                        std::size_t curr_key_freq = 1;
                        result_local.unique_key_count++;

                        if(result_local.mode_count < 1)
                            result_local.mode_count = 1;

                        for(std::size_t i = 1; i < elem_count; ++i)
                            if(p_data[i].first != curr_key)
                            {
                                curr_key = p_data[i].first;
                                curr_key_freq = 1;
                                result_local.unique_key_count++;
                            }
                            else
                                if(result_local.mode_count < ++curr_key_freq)
                                    result_local.mode_count = curr_key_freq;

                        result_local.pair_count += elem_count;
                    }


                    // Write the partition data back to disk.

                    std::remove(p_path.c_str());    // Remove the file, as ext4 fs driver close() waits before data
                                                    // are really written to the disk when done on an *existing* i-node.
                                                    // https://superuser.com/questions/865710/write-to-newfile-vs-overwriting-performance-issue
                    std::ofstream output(p_path.c_str(), std::ios::out | std::ios::binary);
                    if(!output.write(reinterpret_cast<const char*>(p_data), p_bytes))
                    {
                        std::cerr << "Error writing to the partition files. Aborting.\n";
                        std::exit(EXIT_FAILURE);
                    }

                    output.close();
                }

                std::free(p_data);

                result = result_local;
            },
            t_id, thread_count, std::ref(worker_aggregate[t_id])
        );


    for(uint32_t t_id = 0; t_id < thread_count; ++t_id)
    {
        if(!worker[t_id].joinable())
        {
            std::cerr << "Early termination encountered for a collator thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[t_id].join();
        agg_result.aggregate(worker_aggregate[t_id]);
    }
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Key_Value_Iterator<T_key_, T_val_> Key_Value_Collator<T_key_, T_val_, T_hasher_>::begin() const
{
    return iter_t(work_file_pref, partition_count);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Key_Value_Iterator<T_key_, T_val_> Key_Value_Collator<T_key_, T_val_, T_hasher_>::end() const
{
    return iter_t(work_file_pref, partition_count, true);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
const char Key_Value_Collator<T_key_, T_val_, T_hasher_>::partition_file_ext[];


// A class to pack aggregation results from `Key_Value_Collator`.
template <typename T_key_, typename T_val_, typename T_hasher_>
class Key_Value_Collator<T_key_, T_val_, T_hasher_>::Aggregate_Result
{
    friend class Key_Value_Collator<T_key_, T_val_, T_hasher_>;

private:

    std::size_t unique_key_count;  // Number of unique keys.
    std::size_t pair_count;        // Total number of key-value pairs.
    std::size_t mode_count;        // Number of pairs with a most frequent key.


    Aggregate_Result(): unique_key_count(0), pair_count(0), mode_count(0)
    {}


public:

    // Aggregates the result statistics with `other`'s ones.
    void aggregate(const Aggregate_Result& other)
    {
        unique_key_count += other.unique_key_count;
        pair_count += other.pair_count;
        if(mode_count < other.mode_count)
            mode_count = other.mode_count;
    }
};


template <typename T_key_>
class Identity_Functor
{
public:

    T_key_ operator()(const T_key_& key) const { return key; }
};

}



#endif
