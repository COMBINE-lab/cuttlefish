
#ifndef MINIMIZER_INSTANCE_ITERATOR_HPP
#define MINIMIZER_INSTANCE_ITERATOR_HPP



#include "Minimizer_Instance.hpp"
#include "Spin_Lock.hpp"
#include "utility.hpp"

#include <cstddef>
#include <cstring>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdio>


// =============================================================================


// A class to iterate over a collection of type `T_container_` of minimizer
// instances.
template <typename T_container_>
class Minimizer_Instance_Iterator
{};


// The minimizer instances are stored in a `Minimizer_Instance*`-type container,
// having its size paired.
template <>
class Minimizer_Instance_Iterator<std::pair<Minimizer_Instance*, std::size_t>>
{
    typedef std::pair<Minimizer_Instance*, std::size_t> T_container_;

private:

    const Minimizer_Instance* const container;  // Pointer to the container.
    const Minimizer_Instance* const container_end;  // Pointer to the container-end.
    const Minimizer_Instance* ptr;  // Pointer to the current minimizer instance in a linear scan over the container.


public:

    Minimizer_Instance_Iterator(const T_container_& container):
        container(container.first),
        container_end(container.first + container.second),
        ptr(this->container)
    {}

    void operator++()
    {
        ptr++;
    }

    const Minimizer_Instance& operator*() const
    {
        return *ptr;
    }

    explicit operator bool() const
    {
        return ptr != container_end;
    }
};


// The minimizer instances are stored in a `FILE*`-type container.
template <>
class Minimizer_Instance_Iterator<std::FILE*>
{
private:

    std::FILE* file_ptr;    // Pointer to the container (file).

    std::size_t pos;    // Absolute index into the container.
    Minimizer_Instance* buffer; // Buffer to read in chunks of minimizer instances.
    std::size_t buf_elem_count; // Number of instances currently in the buffer.
    std::size_t buf_idx;    // Index of the next instance to process from the buffer.
    constexpr static std::size_t buf_sz = 5 * 1024 * 1024 / sizeof(Minimizer_Instance); // Size of the buffer (in elements): total 5 MB.

    Minimizer_Instance elem;    // Current instance to process.

    Spin_Lock lock; // Mutually-exclusive access lock for the consumers.


    // Sets the file pointer `file_ptr` appropriately at the `pos`'th indexed
    // minimizer instance in the file.
    void set_file_ptr();

    // Peeks into the file for a byte, without consuming it. Sets the file-
    // handle to null if the end-of-file has been reached.
    void peek();

    // Advances in the file by one minimizer instance. Sets the file-handle to
    // null if the end-of-file has been reached.
    void advance();

    // Advances in the file by one minimizer-block, i.e. passes by all the
    // instances that have the same minimizer as the current one.
    void advance_minimizer_block();


public:

    // Constructs an empty iterator, with a null file-handle.
    Minimizer_Instance_Iterator();

    // Constructs an iterator for the file container `in_file_ptr`.
    Minimizer_Instance_Iterator(std::FILE* file_ptr);

    // Copy constructs an iterator from another iterator `other`.
    Minimizer_Instance_Iterator(const Minimizer_Instance_Iterator& other);

    // Destructs the iterator.
    ~Minimizer_Instance_Iterator();

    // Returns the minimizer value of the current minimizer instance.
    cuttlefish::minimizer_t operator*();

    // Advances the iterator by one minimizer block in the container.
    Minimizer_Instance_Iterator& operator++();

    // Returns `true` iff this iterator and `rhs` point to the same position of
    // the same file.
    bool operator==(const Minimizer_Instance_Iterator& rhs) const;

    // Returns `true` iff this iterator and `rhs` point to different positions
    // of some file(s).
    bool operator!=(const Minimizer_Instance_Iterator& rhs) const;

    // Advances the iterator by one minimizer block in the container. Puts the
    // minimizer of the block to `min` and its satellite data (of type
    // `T_satellite_`) to `satellite`. Returns `true` iff there were blocks
    // remaining, i.e. the end-of-file had not been reached.
    template <typename T_satellite_>
    bool next(cuttlefish::minimizer_t& min, T_satellite_& count);

    // Tries to read in at most `buf_sz` minimizer instances into `buf`.
    // Returns the number of elements read, which is 0 in case when the EOF has
    // been reached. It is thread-safe.
    std::size_t next(Minimizer_Instance* buf, std::size_t buf_sz);

    // Dummy methods.
    void launch_production() {}
    bool launched() { return true; }
    bool value_at(const size_t consumer_id, cuttlefish::minimizer_t& elem) { (void)consumer_id; (void)elem; return false; }
    bool tasks_expected(const size_t consumer_id) const { (void)consumer_id; return false; }
    void seize_production() {}
};


inline cuttlefish::minimizer_t Minimizer_Instance_Iterator<std::FILE*>::operator*()
{
    if(buffer == nullptr)
        advance();

    return elem.minimizer();
}


inline Minimizer_Instance_Iterator<std::FILE*>& Minimizer_Instance_Iterator<std::FILE*>::operator++()
{
    advance_minimizer_block();
    return *this;
}


inline bool Minimizer_Instance_Iterator<std::FILE*>::operator==(const Minimizer_Instance_Iterator& rhs) const
{
    if(file_ptr == nullptr || rhs.file_ptr == nullptr)
        return file_ptr == nullptr && rhs.file_ptr == nullptr;

    return file_ptr == rhs.file_ptr && pos == rhs.pos;
}


inline bool Minimizer_Instance_Iterator<std::FILE*>::operator!=(const Minimizer_Instance_Iterator& rhs) const
{
    return !operator==(rhs);
}


template <>
inline bool Minimizer_Instance_Iterator<std::FILE*>::next<std::size_t>(cuttlefish::minimizer_t& min, std::size_t& count)
{
    if(file_ptr == nullptr)
        return false;

    if(buffer == nullptr)
        advance();

    count = 0;
    min = elem.minimizer();
    while(file_ptr != nullptr && elem.minimizer() == min)
    {
        advance();
        count++;
    }

    return true;
}


template <>
inline bool Minimizer_Instance_Iterator<std::FILE*>::next<std::vector<std::size_t>>(cuttlefish::minimizer_t& min, std::vector<std::size_t>& offsets)
{
    if(file_ptr == nullptr)
        return false;

    if(buffer == nullptr)
        advance();

    offsets.clear();
    min = elem.minimizer();
    while(file_ptr != nullptr && elem.minimizer() == min)
    {
        offsets.emplace_back(elem.offset());
        advance();
    }

    return true;
}


inline std::size_t Minimizer_Instance_Iterator<std::FILE*>::next(Minimizer_Instance* const buf, const std::size_t buf_sz)
{
    lock.lock();

    const std::size_t buf_elem_count = (file_ptr != nullptr ?
                                            std::fread(static_cast<void*>(buf), sizeof(Minimizer_Instance), buf_sz, file_ptr) :
                                            0);
    if(buf_elem_count == 0)
        file_ptr = nullptr;

    pos += buf_elem_count;

    lock.unlock();

    return buf_elem_count;
}


inline void Minimizer_Instance_Iterator<std::FILE*>::advance()
{
    if(buffer == nullptr)
        buffer = allocate<Minimizer_Instance>(buf_sz);

    if(buf_idx >= buf_elem_count)
    {
        buf_elem_count = std::fread(static_cast<void*>(buffer), sizeof(Minimizer_Instance), buf_sz, file_ptr);
        buf_idx = 0;

        if(buf_elem_count == 0)
        {
            file_ptr = nullptr;
            return;
        }
    }

    elem = buffer[buf_idx++];
    pos++;
}


inline void Minimizer_Instance_Iterator<std::FILE*>::advance_minimizer_block()
{
    if(buffer == nullptr)
        advance();

    const cuttlefish::minimizer_t min = elem.minimizer();
    while(file_ptr != nullptr && min == elem.minimizer())
        advance();
}




#endif
