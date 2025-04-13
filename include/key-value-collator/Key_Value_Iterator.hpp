
#ifndef KEY_VALUE_ITERATOR_HPP
#define KEY_VALUE_ITERATOR_HPP



#include "Spin_Lock.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include <iostream>


// =============================================================================

namespace key_value_collator
{


// A class to iterate over key-value pairs of type `(T_key_, T_val_)`, collated
// by the class `Key_Value_Collator`.
template <typename T_key_, typename T_val_>
class Key_Value_Iterator
{
    template <typename, typename, typename> friend class Key_Value_Collator;

    typedef std::pair<T_key_, T_val_> key_val_pair_t;

private:

    const std::string work_pref;    // Path to the working files used by the collator.
    const std::size_t partition_count;  // Number of partitions used by the collator.

    std::FILE* file_ptr;    // Pointer to the current (collated) partition file.
    std::size_t curr_p_id;  // ID of the current partition being iterated over.

    std::size_t pos;    // Absolute index (sequential ID of the key-block) into the collated collection.
    bool at_end;    // Whether the iterator is at the end of the collection.

    static constexpr std::size_t buf_sz = 5lu * 1024lu * 1024lu / sizeof(key_val_pair_t);   // Size of the buffer in elements: total 5MB.
    key_val_pair_t* buf;    // Buffer to read in chunks of key-value pairs.
    std::size_t buf_elem_count; // Number of pairs currently in the buffer.
    std::size_t buf_idx;    // Index of the next pair to process from the buffer.

    key_val_pair_t elem;    // Current pair to process.

    Spin_Lock lock; // Mutually-exclusive access lock for the iterator users.


    // Constructs an iterator for a key-value collator that has its collated
    // files at path prefix `work_pref` and used `partition_count` partitions.
    Key_Value_Iterator(const std::string& work_pref, std::size_t partition_count, bool at_end = false);

    // Returns the disk-file path for the partition `p_id`.
    // TODO: think management of codes repeated between collator and iterator.
    const std::string partition_file_path(std::size_t p_id) const { return work_pref + "." + std::to_string(p_id) + ".part"; }

    // Sets the file-handle to the file of the partition with ID `p_id`.
    void set_file_handle(std::size_t p_id);

    // Advances in the collection by one key-value pair. Sets the current file-
    // handle to null if the end of the collection has been reached.
    void advance();

    // Advances in the collection by one key-block, i.e. passes by all the
    // pairs that have the same key as the current one.
    void advance_key_block();


public:

    // Copy constructs an iterator from the iterator `other`. Only usable with
    // `other` iterators that are unmodified results of `begin()` and `end()`.
    Key_Value_Iterator(const Key_Value_Iterator& other);

    // Destructs the iterator.
    ~Key_Value_Iterator();

    // Returns the key of the current pair.
    T_key_ operator*();

    // Advances the iterator by one key-block.
    Key_Value_Iterator& operator++();

    // Returns `true` iff this iterator and `rhs` point to the same key-block of
    // the same collection.
    bool operator==(const Key_Value_Iterator& rhs) const;

    // Returns `true` iff this iterator and `rhs` point to different key-blocks
    // of some collection(s).
    bool operator!=(const Key_Value_Iterator& rhs) const;

    // Returns the absolute key-value pair index of the iterator's current
    // position.
    std::size_t pair_index() const { return pos; }

    // Tries to read in at most `count` key-value pairs into `buf`. Returns the
    // number of pairs read, which is 0 in case when the end of the collection
    // has been reached. It is thread-safe.
    std::size_t read(key_val_pair_t* buf, std::size_t count);

    // Dummy methods.
    void launch_production() {}
    bool launched() { return true; }
    bool value_at(size_t, cuttlefish::minimizer_t&) { return false; }
    bool tasks_expected(const size_t) const { return false; }
    void seize_production() {}
};


template <typename T_key_, typename T_val_>
inline Key_Value_Iterator<T_key_, T_val_>::Key_Value_Iterator(const std::string& work_pref, const std::size_t partition_count, const bool at_end):
    work_pref(work_pref),
    partition_count(partition_count),
    file_ptr(nullptr),
    curr_p_id(0),
    pos(0),
    at_end(at_end),
    buf(nullptr),
    buf_elem_count(0),
    buf_idx(0)
{}


template <typename T_key_, typename T_val_>
inline Key_Value_Iterator<T_key_, T_val_>::Key_Value_Iterator(const Key_Value_Iterator& other):
    work_pref(other.work_pref),
    partition_count(other.partition_count),
    file_ptr(other.file_ptr),
    curr_p_id(other.curr_p_id),
    pos(other.pos),
    at_end(other.at_end),
    buf(nullptr),
    buf_elem_count(other.buf_elem_count),
    buf_idx(other.buf_idx)
{
    if(other.file_ptr != nullptr || other.buf != nullptr)
    {
        std::cerr << "Cannot copy key-value iterator that is in use. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_key_, typename T_val_>
inline Key_Value_Iterator<T_key_, T_val_>::~Key_Value_Iterator()
{
    if(file_ptr != nullptr)
        std::fclose(file_ptr);

    std::free(buf);
}


template <typename T_key_, typename T_val_>
inline T_key_ Key_Value_Iterator<T_key_, T_val_>::operator*()
{
    if(buf == nullptr)
        advance();

    return elem.first;
}


template <typename T_key_, typename T_val_>
inline Key_Value_Iterator<T_key_, T_val_>& Key_Value_Iterator<T_key_, T_val_>::operator++()
{
    advance_key_block();

    return *this;
}



template <typename T_key_, typename T_val_>
inline bool Key_Value_Iterator<T_key_, T_val_>::operator==(const Key_Value_Iterator& rhs) const
{
    if(file_ptr == nullptr || rhs.file_ptr == nullptr)
        return file_ptr == nullptr && rhs.file_ptr == nullptr && at_end == rhs.at_end;

    return file_ptr == rhs.file_ptr && pos == rhs.pos;
}


template <typename T_key_, typename T_val_>
inline bool Key_Value_Iterator<T_key_, T_val_>::operator!=(const Key_Value_Iterator& rhs) const
{
    return !this->operator==(rhs);
}


template <typename T_key_, typename T_val_>
inline std::size_t Key_Value_Iterator<T_key_, T_val_>::read(key_val_pair_t* const buf, const std::size_t count)
{
    lock.lock();

    if(at_end)
    {
        lock.unlock();
        return 0;
    }

    if(file_ptr == nullptr && !at_end)
        set_file_handle(0);

    std::size_t buf_elem_count = 0;
    while(buf_elem_count == 0)
    {
        buf_elem_count = (file_ptr != nullptr ?
                            std::fread(static_cast<void*>(buf), sizeof(key_val_pair_t), count, file_ptr) : 0);

        if(buf_elem_count == 0)
        {
            if(++curr_p_id == partition_count)
            {
                std::fclose(file_ptr);

                file_ptr = nullptr;
                at_end = true;
                break;
            }

            set_file_handle(curr_p_id);
        }
    }

    pos += buf_elem_count;

    lock.unlock();

    return buf_elem_count;
}


template <typename T_key_, typename T_val_>
inline void Key_Value_Iterator<T_key_, T_val_>::advance()
{
    if(buf == nullptr)
    {
        set_file_handle(0);
        buf = static_cast<key_val_pair_t*>(std::malloc(buf_sz * sizeof(key_val_pair_t)));
    }

    while(buf_idx >= buf_elem_count)
    {
        buf_elem_count = std::fread(static_cast<void*>(buf), sizeof(key_val_pair_t), buf_sz, file_ptr);
        buf_idx = 0;

        if(buf_elem_count == 0)
        {
            if(++curr_p_id == partition_count)
            {
                std::fclose(file_ptr);

                file_ptr = nullptr;
                at_end = true;
                return;
            }

            set_file_handle(curr_p_id);
        }
    }

    elem = buf[buf_idx++];
    pos++;
}


template <typename T_key_, typename T_val_>
inline void Key_Value_Iterator<T_key_, T_val_>::advance_key_block()
{
    if(buf == nullptr)
        advance();

    const T_key_ key = elem.first;
    while(file_ptr != nullptr && elem.first == key)
        advance();
}


template <typename T_key_, typename T_val_>
inline void Key_Value_Iterator<T_key_, T_val_>::set_file_handle(const std::size_t p_id)
{
    if(file_ptr != nullptr)
        std::fclose(file_ptr);

    file_ptr = std::fopen(partition_file_path(p_id).c_str(), "rb");
    if(file_ptr == nullptr)
    {
        std::cerr << "Error opening partition files. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

}



#endif
