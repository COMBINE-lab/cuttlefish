
#ifndef STRING_BUFFER_HPP
#define STRING_BUFFER_HPP



#include "Spin_Lock.hpp"

#include <vector>
#include <fstream>
#include <iostream>


// A buffer class to contain contiguous strings. The buffer is to have a maximum
// capacity of `CAPACITY` (although it is non-binding when a string with length 
// larger than that is added), and it flushes to a sink of type `T_sink_` when it
// overflows or is destructed. Writing to the provided sink (in the constructor)
// is thread-safe — w/ a limiting contention for access per sink-type, not per sink.
template <std::size_t CAPACITY, typename T_sink_>
class String_Buffer
{
private:

    std::vector<char> buffer;   // The string buffer.
    T_sink_& sink;  // Reference to the sink to flush the buffer content to.


    // Flushes the buffer content to the sink, and clears the buffer.
    void flush();


public:

    // Constructs a string buffer object that would flush its content to `sink`.
    String_Buffer(T_sink_& sink);

    // Appends the content of the string `str` to the buffer. Flushes are possible.
    void operator+=(const std::string& str);

    // Destructs the buffer object, flushing it if content are present.
    ~String_Buffer();
};


// Helper class to actually flush the content of the `String_Buffer` class to its
// sink of type `T_sink`.
// It's used to circumvent the C++ constraint that partial specialization of a
// a member function is not possible without partially specializing the entire
// class. We need to specialize the actual flushing mechanism to support various
// types of sinks, e.g. `std::ofstream`, `spdlog::logger` etc.
template <typename T_sink_>
class String_Buffer_Flusher
{
    // Since the sole purpose of the class is to support the `String_Buffer` class
    // circumvent some contraint, everything is encapsulated here as private, with
    //  `String_Buffer` as friend.
    template <std::size_t, typename> friend class String_Buffer;

private:

    // Mutual-exclusion lock to control multi-threaded access to otherwise not thread-
    // safe sinks (e.g. `std::ofstream`). Note that, the lock is per sink-type, not per
    // actual sink — which is a limitation.
    static Spin_Lock lock;


    // Writes `len` characters from the memory location `str_buf` to the sink `sink`.
    static void write(const char* str_buf, std::size_t len, T_sink_& sink);
};


template <std::size_t CAPACITY, typename T_sink_>
inline String_Buffer<CAPACITY, T_sink_>::String_Buffer(T_sink_& sink):
    sink(sink)
{
    buffer.reserve(CAPACITY);
}


template <std::size_t CAPACITY, typename T_sink_>
inline void String_Buffer<CAPACITY, T_sink_>::operator+=(const std::string& str)
{
    if(buffer.size() + str.length() >= CAPACITY)
    {
        flush();
        
        if(str.length() >= CAPACITY)
        {
            std::cerr <<    "A single output string overflows the string-buffer capacity.\n"
                            "Output string length: " << str.length() << ", string-buffer capacity: " << CAPACITY << ".\n"
                            "Please consider increasing the buffer capacity parameter in build for future use.\n";
            
            buffer.reserve(str.length());
        }
    }


    buffer.insert(buffer.end(), str.begin(), str.end());
}


template <std::size_t CAPACITY, typename T_sink_>
inline void String_Buffer<CAPACITY, T_sink_>::flush()
{
    String_Buffer_Flusher<T_sink_>::write(buffer.data(), buffer.size(), sink);

    buffer.clear();
}


template <std::size_t CAPACITY, typename T_sink_>
inline String_Buffer<CAPACITY, T_sink_>::~String_Buffer()
{
    if(!buffer.empty())
        flush();
}


template <>
inline void String_Buffer_Flusher<std::ofstream>::write(const char* const str_buf, const std::size_t len, std::ofstream& output)
{
    lock.lock();

    output.write(str_buf, len);

    if(output.fail())
    {
        std::cerr << "Error writing the output. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    lock.unlock();
}


template <typename T_sink_> Spin_Lock String_Buffer_Flusher<T_sink_>::lock; // Definition of the static lock of `String_Buffer_Flusher`.



#endif
