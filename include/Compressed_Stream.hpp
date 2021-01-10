
#ifndef COMPRESSED_STREAM_HPP
#define COMPRESSED_STREAM_HPP



// #include "zstr/zstr.hpp"
#include "lz4_stream/lz4_stream.h"
#include "FastPFor/codecfactory.h"

#include <fstream>


template <typename T_elem_>
class Compressed_Stream
{
private:

    // std::unique_ptr<std::ostream> output{nullptr};  //
    // std::unique_ptr<std::istream> input{nullptr};   //
    std::unique_ptr<std::ofstream> op{nullptr}; //
    std::unique_ptr<std::ifstream> ip{nullptr}; //
    std::unique_ptr<lz4_stream::ostream> output{nullptr};   //
    std::unique_ptr<lz4_stream::istream> input{nullptr};    //


public:
    
    // 
    Compressed_Stream(const std::string& file_path, bool is_compression = true);

    // 
    ~Compressed_Stream();

    // 
    std::streamsize read(void* buf, size_t elem_count) const;

    // 
    void write(void* buf, size_t elem_count) const;
};


template <typename T_elem_>
inline Compressed_Stream<T_elem_>::Compressed_Stream(const std::string& file_path, const bool is_compression)
{
    // if(is_compression)
    //     output = std::unique_ptr<std::ostream>(new zstr::ofstream(file_path));
    // else
    //     input = std::unique_ptr<std::istream>(new zstr::ifstream(file_path));

    if(is_compression)
    {
        op = std::unique_ptr<std::ofstream>(new std::ofstream(file_path));
        output = std::unique_ptr<lz4_stream::ostream>(new lz4_stream::ostream(*op));
    }
    else
    {
        ip = std::unique_ptr<std::ifstream>(new std::ifstream(file_path));
        input = std::unique_ptr<lz4_stream::istream>(new lz4_stream::istream(*ip));
    }
}


template <typename T_elem_>
inline Compressed_Stream<T_elem_>::~Compressed_Stream()
{
    output.reset(nullptr);
    op.reset(nullptr);

    ip.reset(nullptr);
    input.reset(nullptr);
}


template <typename T_elem_>
inline std::streamsize Compressed_Stream<T_elem_>::read(void* const buf, const size_t elem_count) const
{
    assert(!is_compression);

    char* const buffer = (char*)buf;
    input->read(buffer, elem_count * sizeof(T_elem_));
    const std::streamsize elems_read = input->gcount() / sizeof(T_elem_);
    
    return elems_read;
}


template <typename T_elem_>
inline void Compressed_Stream<T_elem_>::write(void* const buf, const size_t elem_count) const
{
    assert(is_compression);

    char* const buffer = (char*)buf;
    output->write(buffer, elem_count * sizeof(T_elem_));
}



#endif
