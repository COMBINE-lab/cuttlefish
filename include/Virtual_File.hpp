
#ifndef VIRTUAL_FILE_HPP
#define VIRTUAL_FILE_HPP



#include "utility.hpp"

#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <cstdlib>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// In-memory file accessor with bounded memory. Supports only non-decreasing
// indexing into the file.
template <typename T_>
class Virtual_File
{
private:

    static constexpr std::size_t buf_sz_default = 16 * 1024;    // 16KB.

    const std::size_t buf_sz;   // Maximum number of bytes from the file to keep in memory.
    const std::size_t buf_elem_count;   // Maximum number of elements from the file to keep in memory.
    const std::size_t file_elem_count;  // Number of elements in the file.

    T_* const buf;  // The in-memory file buffer.

    std::size_t chunk_start_idx;    // Index into the file where the chunk currently loaded into the buffer starts.
    std::size_t chunk_end_idx;  // Non-inclusive index into the file where the chunk currently loaded into the buffer ends.

    std::FILE* const fp;    // Handle to the file.

    std::size_t next_acc_idx;   // Next valid index to access into the file; used for error-checking.


    // Reads in a chunk of data from the file into the buffer and returns the
    // number of elements read.
    std::size_t read();


public:

    // Constructs a virtual-file for the file at path `file_path`. It uses
    // `buf_bytes` number of bytes.
    Virtual_File(const char* file_path, std::size_t buf_bytes = buf_sz_default);

    ~Virtual_File();

    // Returns the size of the file in elements.
    std::size_t size() const { return file_elem_count; }

    // Returns the data at index `idx` of the file.
    T_ operator[](std::size_t idx);


    // Invalidate move- and copy-constructors, and copy-assignment.
    Virtual_File(Virtual_File&&) = delete;
    Virtual_File(const Virtual_File&) = delete;
    Virtual_File& operator=(const Virtual_File&) = delete;
    Virtual_File& operator=(Virtual_File&) = delete;
};


template <typename T_>
inline Virtual_File<T_>::Virtual_File(const char* const file_path, const std::size_t buf_bytes):
      buf_sz(buf_bytes)
    , buf_elem_count(buf_sz / sizeof(T_))
    , file_elem_count(std::filesystem::file_size(file_path) / sizeof(T_))
    , buf(allocate<T_>(buf_sz))
    , chunk_start_idx(0)
    , chunk_end_idx(0)
    , fp(std::fopen(file_path, "rb"))
    , next_acc_idx(0)
{
    assert(buf_elem_count > 0);
    assert(std::filesystem::file_size(file_path) % sizeof(T_) == 0);
}


template <typename T_>
inline Virtual_File<T_>::~Virtual_File()
{
    deallocate(buf);

    std::fclose(fp);
}


template <typename T_>
inline std::size_t Virtual_File<T_>::read()
{
    const std::size_t elems_to_read = std::min(file_elem_count - chunk_end_idx, buf_elem_count);
    const std::size_t elems_read = std::fread(buf, sizeof(T_), elems_to_read, fp);
    assert(elems_read == elems_to_read);
    (void)elems_read;

    return elems_to_read;
}


template <typename T_>
inline T_ Virtual_File<T_>::operator[](const std::size_t idx)
{
    assert(idx >= next_acc_idx);
    next_acc_idx = idx;

    while(idx >= chunk_end_idx)
    {
        chunk_start_idx = chunk_end_idx;
        chunk_end_idx += read();
    }

    return buf[idx - chunk_start_idx];
}

}



#endif
