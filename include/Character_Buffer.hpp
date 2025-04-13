
#ifndef CHARACTER_BUFFER_HPP
#define CHARACTER_BUFFER_HPP



#include "Spin_Lock.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "FASTA_Record.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>


// A buffer class to contain contiguous characters. It flushes to a sink of
// type `T_sink_` when it overflows or is destructed. Writing to the provided
// sink (in the constructor) is thread-safe.
template <typename T_sink_>
class Character_Buffer
{
private:

    // The buffer is to have a maximum capacity of `CAPACITY` (it is non-
    // binding when a string with length  larger than that is added).
    static constexpr std::size_t cap_ = 128 * 1024ULL;

    std::string buf;    // The character buffer.
    T_sink_& sink;  // Reference to the sink to flush the buffer content to.


    // Ensures that the buffer has enough space for additional `append_size`
    // number of bytes, using flush if required.
    void ensure_space(std::size_t append_size);

    // Flushes the buffer content to the sink, and clears the buffer.
    void flush();


public:

    // Constructs a character buffer object that would flush its content to `sink`.
    Character_Buffer(T_sink_& sink);

    Character_Buffer(const Character_Buffer& rhs) = default;

    Character_Buffer(Character_Buffer&& rhs) = default;

    // Appends the content of `str` to the buffer. Flushes are possible.
    template <typename T_container_>
    void operator+=(const T_container_& str);

    // Appends the content of the FASTA record `fasta_rec` to the buffer. Flushes
    // are possible.
    template <bool Colored_ = false>
    void operator+=(const FASTA_Record& fasta_rec);

    // Appends the content of the FASTA record `fasta_cycle` to the buffer, that is
    // supposed to be a cycle in a de Bruijn graph `G(·, k)`. The cyclic FASTA
    // sequence is rotated around its index `pivot` — the entire sequence is
    // right-rotated so that the `pivot`-index character is at index 0 finally. A
    // line-break is added at the end of the sequence, since the user might not be
    // able to provide it with the "to be rotated" sequence.
    template <uint16_t k>
    void rotate_append_cycle(const FASTA_Record& fasta_cycle, std::size_t pivot);

    // Returns the `len`-length suffix of the buffer.
    const char* suffix(std::size_t len) const;

    // Flushes the buffer if not empty.
    void close();

    // Destructs the buffer object, flushing it if content are present.
    ~Character_Buffer();
};


// Helper class to actually flush the content of the `Character_Buffer` class to its
// sink of type `T_sink`.
// It's used to circumvent the C++ constraint that partial specialization of a
// a member function is not possible without partially specializing the entire
// class. We need to specialize the actual flushing mechanism to support various
// types of sinks, e.g. `std::ofstream`, `spdlog::logger` etc.
// Since the sole purpose of the class is to support the `Character_Buffer` class
// circumvent some contraint, everything is encapsulated in its specializations
// as private, with `Character_Buffer` as friend.
// TODO: remove this entire class as multi-specializations are not used, which
// was its sole purpose.
template <typename T_sink_>
class Character_Buffer_Flusher
{};


template <>
class Character_Buffer_Flusher<std::ofstream>
{
    template <typename> friend class Character_Buffer;

private:

    // Mutual-exclusion lock to control multi-threaded access to otherwise not thread-
    // safe sinks (e.g. `std::ofstream`). Note that, the lock is per sink-type, not per
    // actual sink — which is a limitation.
    static Spin_Lock lock;


    // Writes the content of `buf` to the sink `sink`.
    static void write(const std::string& buf, std::ofstream& sink);
};


template <>
class Character_Buffer_Flusher<Async_Logger_Wrapper>
{
    template <typename> friend class Character_Buffer;

private:

    // Writes the content of `buf` to the sink `sink`. Note that
    // `buf` is modified in the process — a null-terminator (`\0`) is appended at the
    // end — which is expected to be not problematic under the assumption that the
    // buffer is cleared after the write (i.e. flush).
    static void write(std::string& buf, const Async_Logger_Wrapper& sink);
};


template <typename T_sink_>
inline Character_Buffer<T_sink_>::Character_Buffer(T_sink_& sink):
    sink(sink)
{
    buf.reserve(cap_);
}


template <typename T_sink_>
template <typename T_container_>
inline void Character_Buffer<T_sink_>::operator+=(const T_container_& str)
{
    ensure_space(str.size());
    buf.append(str.begin(), str.end());
}


template <typename T_sink_>
template <bool Colored_>
inline void Character_Buffer<T_sink_>::operator+=(const FASTA_Record& fasta_rec)
{
    if constexpr(!Colored_)
        ensure_space(fasta_rec.header_size() + 1 + fasta_rec.seq_size() + 1);   // Two extra bytes for the line-breaks.
    else
        ensure_space(fasta_rec.header_size() + (18 + 1) * fasta_rec.color_list_size() + 1 + fasta_rec.seq_size() + 1);

    fasta_rec.append_header(buf);   // Append the header.
    if constexpr(Colored_)
        fasta_rec.append_color_list(buf);
    buf.push_back('\n');    // Break line.
    fasta_rec.append_seq(buf);  // Append the sequence.
    buf.push_back('\n');    // Break line.
}


template <typename T_sink_>
template <uint16_t k>
inline void Character_Buffer<T_sink_>::rotate_append_cycle(const FASTA_Record& fasta_rec, const std::size_t pivot)
{
    ensure_space(fasta_rec.header_size() + 1 + fasta_rec.seq_size() + 1);   // Two extra bytes for two line-breaks.

    fasta_rec.append_header(buf);   // Append the header.
    buf.push_back('\n');    // Break line.
    fasta_rec.template append_rotated_cycle<k>(buf, pivot); // Append the sequence right-rotated around index `pivot`.
    buf.push_back('\n');    // End the sequence.
}


template <typename T_sink_>
inline void Character_Buffer<T_sink_>::ensure_space(const std::size_t append_size)
{
    if(buf.size() + append_size >= cap_)    // Using `>=` since for async logging, a `\0` is inserted at the end of `buffer`.
        flush();
}


template <typename T_sink_>
inline const char* Character_Buffer<T_sink_>::suffix(const std::size_t len) const
{
    return buf.data() + (buf.size() - len);
}


template <typename T_sink_>
inline void Character_Buffer<T_sink_>::flush()
{
    Character_Buffer_Flusher<T_sink_>::write(buf, sink);

    buf.clear();
}


template <typename T_sink_>
inline void Character_Buffer<T_sink_>::close()
{
    if(!buf.empty())
        flush();
}


template <typename T_sink_>
inline Character_Buffer<T_sink_>::~Character_Buffer()
{
    close();
}


inline void Character_Buffer_Flusher<std::ofstream>::write(const std::string& buf, std::ofstream& output)
{
    lock.lock();

    output.write(buf.data(), buf.size());
    if(!output)
    {
        std::cerr << "Error writing the output. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    lock.unlock();
}


inline void Character_Buffer_Flusher<Async_Logger_Wrapper>::write(std::string& buf, const Async_Logger_Wrapper& sink)
{
    buf.push_back('\0');

    sink.write(buf.data());
}


#endif
