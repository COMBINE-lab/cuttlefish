
#ifndef CODEC_STREAM_HPP
#define CODEC_STREAM_HPP



#include "SpinLock/SpinLock.hpp"
#include "lz4/lz4.h"

#include <vector>
#include <cstddef>
#include <numeric>


namespace codec_stream
{

// Forward declaration of the class — required for the `Token` class.
template <typename T_elem_> class Codec_Stream;


// A token class to distinguish different user-threads of the same `Codec_Stream` object.
class Token
{
    template <typename T_elem_> friend class Codec_Stream;

private:

    size_t token_num_; // Token-number.


    // Constructs a Token object with the token-number `token`.
    Token(size_t token): token_num_{token}
    {}

    // Returns the token-number.
    size_t token_num() const
    {
        return token_num_;
    }
};


// A codec (compression-decompression) stream class.
template <typename T_elem_>
class Codec_Stream
{
    // Type of the codec block sizes (in bytes). It must be able to hold `max_elems_per_blk * sizeof(T_elem_)`.
    typedef uint32_t blk_sz_t;

private:

    const size_t user_count;    // Number of threads using the codec-stream.
    const size_t max_elems_per_blk;    // Maximum size (in elements) of each buffer (per-user).
    const size_t buf_sz;    // Size (in bytes) of each buffer (per-user).
    const std::string file_pref;    // Prefix of the path to the codec files.
    const std::string stream_file;  // Path to the codec-streaming file.
    const std::string blk_elem_file;    // Path to the metafile containing the number of elements in the codec blocks.
    const std::string blk_bytes_file;   // Path to the metafile containing the number of bytes in the compressed blocks.
    const bool is_compression;  // Whether the codec-stream is being used for compression or decompression.

    std::vector<std::vector<uint8_t>> buffers;  // Collection of buffers, one for each user of the stream.
    std::vector<Token> tokens;  // Collection of user-tokens.
    size_t tokens_used{0};  // Number of tokens produced.

    std::vector<blk_sz_t> elems_per_blk;    // Number of elements in the codec blocks.
    std::vector<blk_sz_t> bytes_per_blk;    // Number of bytes in the compressed blocks.
    size_t read_idx{0}; // Index of the next block to read from disk.

    SpinLock mutex_lock;    // Mutual exclusion lock for accessing disk files.
    std::FILE* stream{nullptr}; // File pointer to the compressed output.

    std::vector<uint64_t> elem_bytes;   // Total number of bytes in the elements provided (per-user).
    std::vector<uint64_t> compressed_bytes; // Total number of bytes in the compressed output (per-user).
    std::vector<uint64_t> decompressed_bytes;   // Total number of bytes in the decompressed output (per-user).

    // The larger the acceleration value, the faster the algorithm, but also the lesser the compression.
    // It's a trade-off. It can be fine tuned, with each successive value providing roughly +~3% to speed.
    // An acceleration value of "1" is the same as regular `LZ4_compress_default()`.
    static constexpr int LZ4_ACCELERATION_FACTOR = 1;   // "Acceleration" factor for the LZ4 compression algorithm.


    // Writes the content of the vector `vec` to the file at path `file_path`.
    template <typename T_vec_elem_> static void write(const std::string& file_path, std::vector<T_vec_elem_>& vec);

    // Reads the content of the file at path `file_path` to the vector `vec`.
    template <typename T_vec_elem_> static void read(const std::string& file_path, std::vector<T_vec_elem_>& vec);

    // Closes compression stream.
    void close_write();

    // Closes decompression stream.
    void close_read();


public:

    // Constructs a `Codec_Stream` that can be used concurrently by up-to `user_count` threads,
    // supporting codec operations of up-to `max_elems_per_blk` number of elements at a time.
    // The stream uses three files — each having the path prefix `file_pref`. The stream works
    // as a compression stream iff `is_compression` is `true`, and works as a decompression
    // stream otherwise.
    Codec_Stream(size_t user_count, size_t max_elems_per_blk, const std::string& file_pref, bool is_compression = true);

    // Destructs the codec-stream.
    ~Codec_Stream();

    // Returns a free and unique token that needs to be passed when using the codec-stream.
    const Token& get_token();

    // Writes (compressed) `elem_count` number of elements from the memory location `data` to
    // the codec-stream file. `token` is used to recognize the invoking user thread.
    void write(const Token& token, const T_elem_* data, size_t elem_count);

    // Writes (compressed) `elem_count` number of elements from the vector `data` to the
    // codec-stream file. `token` is used to recognize the invoking user thread.
    void write(const Token& token, const std::vector<T_elem_>& data, size_t elem_count);

    // Reads up-to `max_elems_per_blk` (decompressed) elements from the codec-stream file into
    // the memory location `data`, and returns the number of elements read. `token` is used to
    // recognize the invoking user thread. Returns 0 if end-of-file has been reached.
    size_t read(const Token& token, T_elem_* data);

    // Reads up-to `max_elems_per_blk` (decompressed) elements from the codec-stream file into
    // the vector `data`, and returns the number of elements read. `token` is used to recognize
    // the invoking user thread. Returns 0 if end-of-file has been reached.
    size_t read(const Token& token, std::vector<T_elem_>& data);

    // Closes the codec-stream.
    void close_stream();

    // Delete the codec-stream specific files.
    void remove_files() const;
};


template <typename T_elem_>
inline Codec_Stream<T_elem_>::Codec_Stream(const size_t user_count, const size_t max_elems_per_blk, const std::string& file_pref, const bool is_compression):
    user_count{user_count},
    max_elems_per_blk{max_elems_per_blk},
    buf_sz(LZ4_compressBound(max_elems_per_blk * sizeof(T_elem_))),
    file_pref(file_pref),
    stream_file(file_pref + ".stream"),
    blk_elem_file(file_pref + ".elems"),
    blk_bytes_file(file_pref + ".bytes"),
    is_compression{is_compression}
{
    tokens.reserve(user_count);

    buffers.resize(user_count);
    for(size_t idx = 0; idx < user_count; ++idx)
        buffers[idx].resize(buf_sz);
    
    elem_bytes.resize(user_count), compressed_bytes.resize(user_count), decompressed_bytes.resize(user_count);


    if(is_compression)
        stream = std::fopen(stream_file.c_str(), "wb");
    else
    {
        stream = std::fopen(stream_file.c_str(), "rb");

        read(blk_elem_file, elems_per_blk);
        read(blk_bytes_file, bytes_per_blk);
    }

    if(std::ferror(stream))
    {
        std::cerr << "Error opening codec-stream file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_elem_>
inline Codec_Stream<T_elem_>::~Codec_Stream()
{
    if(stream != nullptr)
        std::fclose(stream);
}


template <typename T_elem_>
inline const Token& Codec_Stream<T_elem_>::get_token()
{
    if(tokens_used == user_count)
    {
        std::cerr << "Number of tokens requested exceeded the consumer-count of the codec-stream. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
    
    mutex_lock.lock();

    tokens.push_back(tokens_used++);
    const Token& token = tokens.back();

    mutex_lock.unlock();

    return token;
}


template <typename T_elem_>
inline void Codec_Stream<T_elem_>::write(const Token& token, const T_elem_* const data, const size_t elem_count)
{
    assert(elem_count <= max_elems_per_blk);

    const size_t token_num = token.token_num(); // User thread-id.
    const char* const src = (const char*)data;  // Source buffer to compress.
    char* const dst = (char*)(buffers[token_num].data());   // Destination buffer for the compressed output.
    const size_t src_size = elem_count * sizeof(T_elem_);   // Number of bytes to compress.

    // Compress.
    // TODO: test using `LZ4_compress_fast_extState` too (maybe performance gains?).
    const int bytes_count = LZ4_compress_fast(src, dst, src_size, buf_sz, LZ4_ACCELERATION_FACTOR);
    if(bytes_count == 0)
    {
        std::cerr << "Compression failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Write the compression output to disk; also record the counts of the compression-input and -output elements.
    mutex_lock.lock();

    elems_per_blk.emplace_back(elem_count);
    bytes_per_blk.emplace_back(bytes_count);
    
    std::fwrite((const void*)dst, 1, bytes_count, stream);
    if(std::ferror(stream))
    {
        std::cerr << "Error writing compressed data. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
    
    mutex_lock.unlock();

    elem_bytes[token_num] += (elem_count * sizeof(T_elem_));
    compressed_bytes[token_num] += bytes_count;
}


template <typename T_elem_>
inline void Codec_Stream<T_elem_>::write(const Token& token, const std::vector<T_elem_>& data, const size_t elem_count)
{
    write(token, data.data(), elem_count);
}


template <typename T_elem_>
inline size_t Codec_Stream<T_elem_>::read(const Token& token, T_elem_* const data)
{
    const size_t token_num = token.token_num(); // User thread-id.
    char* const src = (char*)(buffers[token_num].data());   // Source buffer to decompress.

    // Read the compressed block into memory.
    mutex_lock.lock();

    if(read_idx == elems_per_blk.size())
    {
        mutex_lock.unlock();
        return 0;
    }

    const blk_sz_t blk_elems = elems_per_blk[read_idx];   // Number of elements compressed in the next block.
    const blk_sz_t blk_bytes = bytes_per_blk[read_idx];   // Number of bytes in the next compressed block.
    read_idx++;

    std::fread((void*)src, 1, blk_bytes, stream);
    if(std::ferror(stream))
    {
        std::cerr << "Error reading compressed data. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    mutex_lock.unlock();

    
    // Decompress.
    char* const dst = (char*)data;
    const int bytes_count = LZ4_decompress_safe(src, dst, blk_bytes, buf_sz);
    if(bytes_count <= 0)
    {
        std::cerr << "Decompression failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    const size_t decompressed_elems = bytes_count / sizeof(T_elem_);
    if(decompressed_elems != blk_elems)
    {
        std::cerr << "Decompressed element count does not match the expected count. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    compressed_bytes[token_num] += blk_bytes;
    decompressed_bytes[token_num] += bytes_count;

    return decompressed_elems;
}


template <typename T_elem_>
inline void Codec_Stream<T_elem_>::close_stream()
{
    if(is_compression)
        close_write();
    else
        close_read();
}


template <typename T_elem_>
inline void Codec_Stream<T_elem_>::close_write()
{
    if(stream == nullptr)
        return;

    
    mutex_lock.lock();

    std::fflush(stream);
    if(std::ferror(stream) || std::fclose(stream) != 0)
    {
        std::cerr << "Error with codec-stream file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    stream = nullptr;

    write<typename decltype(elems_per_blk)::value_type>(blk_elem_file, elems_per_blk);
    write<typename decltype(bytes_per_blk)::value_type>(blk_bytes_file, bytes_per_blk);

    mutex_lock.unlock();

    std::cout << "Total input bytes provided:  " << std::accumulate(elem_bytes.begin(), elem_bytes.end(), 0) << "\n";
    std::cout << "Total compressed byte count: " << std::accumulate(compressed_bytes.begin(), compressed_bytes.end(), 0) << "\n";
}


template <typename T_elem_>
inline void Codec_Stream<T_elem_>::close_read()
{
    if(stream == nullptr)
        return;

    mutex_lock.lock();
    
    if(std::fclose(stream) != 0)
    {
        std::cerr << "Error closing codec-stream file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    stream = nullptr;
    mutex_lock.unlock();

    std::cout << "Total compressed bytes read:   " << std::accumulate(compressed_bytes.begin(), compressed_bytes.end(), 0) << "\n";
    std::cout << "Total decompressed byte count: " << std::accumulate(decompressed_bytes.begin(), decompressed_bytes.end(), 0) << "\n";
}


template <typename T_elem_>
inline void Codec_Stream<T_elem_>::remove_files() const
{
    if(std::remove(stream_file.c_str()) != 0 || std::remove(blk_elem_file.c_str()) != 0 || std::remove(blk_bytes_file.c_str()) != 0)
    {
        std::cerr << "Error deleting codec-stream files. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_elem_>
template <typename T_vec_elem_>
inline void Codec_Stream<T_elem_>::write(const std::string& file_path, std::vector<T_vec_elem_>& vec)
{
    std::FILE* fp = std::fopen(file_path.c_str(), "wb");
    std::fwrite(vec.data(), sizeof(T_vec_elem_), vec.size(), fp);
    std::fflush(fp);

    if(std::ferror(fp))
    {
        std::cerr << "Error writing to the codec-stream metafiles. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    std::fclose(fp);
}


template <typename T_elem_>
template <typename T_vec_elem_>
inline void Codec_Stream<T_elem_>::read(const std::string& file_path, std::vector<T_vec_elem_>& vec)
{
    std::ifstream file(file_path.c_str(), std::ios::binary | std::ios::ate);
    std::streamsize file_sz = file.tellg();
    file.seekg(0, std::ios::beg);

    const size_t elem_count = file_sz / sizeof(T_vec_elem_);
    vec.resize(elem_count);

    file.read((char*)vec.data(), file_sz);
    if(file.bad())
    {
        std::cerr << "Error reading codec-stream metafile " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

}


#endif
