
#ifndef UTILITY_HPP
#define UTILITY_HPP



#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <type_traits>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <chrono>

// TODO: wrap everything here in some namespaces.
// =============================================================================

// Returns a random string of length `len`, using characters from `alphabet`.
std::string get_random_string(size_t len, const char* alphabet =    "0123456789"
                                                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                    "abcdefghijklmnopqrstuvwxyz");

// Returns `true` iff `pref` is a prefix of `s`.
bool is_prefix(const std::string& s, const std::string& pref);

// Returns `true` iff there exists a file in the file system with the path
// `file_path`.
bool file_exists(const std::string& file_path);

// Returns `true` iff these exists a directory in the file system with the
// path `dir_path`.
bool dir_exists(const std::string& dir_path);

// Returns the file size is bytes of the file at path `file_path`. Returns
// `0` in case the file does not exist.
std::size_t file_size(const std::string& file_path);

// Returns `true` iff there exists some file in the file system path
// `path` with its name being prefixed by `prefix`.
bool file_prefix_exists(const std::string& path, const std::string& prefix);

// Loads the binary file at path `file_path` and returns its size in bytes.
std::size_t load_file(const char* file_path, char* buf);

// Loads the binary file at path `file_path` and returns its size in bytes.
std::size_t load_file(const std::string& file_path, char* buf);

// Loads `sz` bytes from the binary file at path `file_path`.
void load_file(const std::string& file_path, const std::size_t sz, char* buf);

// Returns a string that is a copy of `s` but has all the whitespaces removed.
std::string remove_whitespaces(const char* s);

// Given the collection of strings `s`, returns the concatenated string
// `s0 : s1 : ... : s_m`, where successive strings are separated by `delimiter`.
const std::string concat_strings(const std::vector<std::string>& s, const std::string& delimiter = ", ");

// Removes the file at path `file_path` from disk. Returns `true` iff the
// removal is successful.
bool remove_file(const std::string& file_path);

// Clears the content of the file at path `file_path`.
void clear_file(const std::string& file_path);

// Returns the name of the file present at the path `file_path`.
const std::string filename(const std::string& file_path);

// Returns the directory of the file present at the path `file_path`.
const std::string dirname(const std::string& file_path);

// Moves the file present at path `from_path` to the path `to_path`.
void move_file(const std::string& from_path, const std::string& to_path);

// Returns the value corresponding to `metric` from the pseudo-filesystem
// `/proc`. Returns `0` in case of errors encountered.
std::size_t process_metric(const std::string& metric);

// Returns the maximum memory ("high-water-mark") used by the running
// process in bytes. Returns `0` in case of errors encountered.
std::size_t process_peak_memory();

// Returns the current memory ("resident-set-size") used by the running process
// in bytes. Returns `0` in case of errors encountered.
std::size_t process_cur_memory();

// Force-frees the memory allocated to the container `container`.
template <typename T_container_>
void force_free(T_container_& container)
{
    T_container_().swap(container);
}

// Returns pointer to a memory-allocation for `size` elements of type `T_`.
template <typename T_>
T_* allocate(const std::size_t size)
{
    return static_cast<T_*>(std::malloc(size * sizeof(T_)));
}

// Returns pointer to a memory-allocation for `size` elements of type `T_`,
// initialized to `0`s.
template <typename T_>
T_* allocate_zeroed(const std::size_t size)
{
    return static_cast<T_*>(std::calloc(size, sizeof(T_)));
}

// Returns pointer to a memory-allocation for `size` elements of type `T_`,
// whose alignment is specified is `alignment`.
template <typename T_>
T_* aligned_allocate(const std::size_t size, const std::size_t alignment = 8)
{
    return static_cast<T_*>(std::aligned_alloc(alignment, size * sizeof(T_)));
}

// Returns pointer to a memory-reallocation for `size` elements of type `T_`.
template <typename T_>
T_* reallocate(T_* const ptr, const std::size_t size)
{
    return static_cast<T_*>(std::realloc(ptr, size * sizeof(T_)));
}

// Deallocates the pointer `ptr`, allocated with `allocate`.
template <typename T_>
void deallocate(T_* const ptr)
{
    std::free(ptr);
}

// Resizes the type-`T_` container `container` geometrically with the growth
// factor `gf` such that it has a size of at least `sz`.
template <typename T_>
void resize_geometric(T_& container, const std::size_t sz, const double gf = 2.0)
{
    assert(gf > 1.0);

    std::size_t curr_sz = std::max(container.size(), 1lu);
    while(curr_sz < sz)
        curr_sz *= gf;

    if(container.size() < curr_sz)
        container.resize(curr_sz);
}

// Allocates the type-`T_` container `p` that currently has space for `cur_sz`
// elements geometrically with the growth factor `gf` such that it has enough
// space for at least `req_sz` elements, and returns the new size. If `keep_`
// is `true`, then the existing elements are kept.
template <typename T_, bool keep_ = false>
std::size_t reserve_geometric(T_*& p, const std::size_t curr_sz, const std::size_t req_sz, const double gf = 2.0)
{
    assert(gf > 1.0);

    if(curr_sz >= req_sz)
        return curr_sz;

    std::size_t new_sz = std::max(curr_sz, 1lu);
    while(new_sz < req_sz)
        new_sz *= gf;

    if constexpr(keep_)
        p = reallocate(p, new_sz);
    else
    {
        deallocate(p);
        p = allocate<T_>(new_sz);
    }

    return new_sz;
}

// Returns the corresponding integer value of type `T_to_` for a type-`T_`
// enum-value `enum_val`.
template <typename T_, typename T_to_ = std::size_t>
constexpr T_to_ as_int(const T_ enum_val)
{
    static_assert(std::is_enum_v<T_>);
    return static_cast<T_to_>(enum_val);
}

// Returns `true` iff `x` is a power of 2.
constexpr bool is_pow_2(const uint64_t x)
{
    return (x > 0 ? ((x & (x - 1)) == 0) : false);
}

// Returns the smallest power of 2 at least as large as `x`. `x` must be in
// `[1, 2^63]`.
constexpr uint64_t ceil_pow_2(uint64_t x)
{
    assert(x > 0 && x <= (1lu << 63));

    if((x & (x - 1)) == 0)
        return x;

    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);

    return x + 1;
}

// Returns the floor of 2-based logarithm of `x`.
constexpr uint64_t log_2(uint64_t x)
{
    assert(x > 0);
    return (x == 1 ? 0 : 1 + log_2(x >> 1));
}


namespace memory
{

// Returns the resident set size of `v`.
template <typename T_>
std::size_t RSS(const std::vector<T_>& v)
{
    static_assert(std::is_standard_layout_v<T_>);
    return sizeof(v) + v.capacity() * sizeof(T_);
}

}


// Wrapper class for a data-element of type `T_` to ensure that in a linear
// collection of `T_`'s, each element is aligned to a cache-line boundary.
template <typename T_>
class alignas(L1_CACHE_LINE_SIZE)
    Padded
{
private:

    T_ data_;


public:

    Padded()
    {}

    Padded(const T_& data):
      data_(data)
    {}

    Padded(T_&& data):
        data_(std::move(data))
    {}

    T_& unwrap() { return data_; }

    const T_& unwrap() const { return data_; }

    template <typename T_archive_> void serialize(T_archive_& archive) { archive(data_); }
};


// Wrapper class for a buffer of elements of type `T_`.
template <typename T_>
class Buffer
{
private:

    std::size_t cap_;   // Capacity of the buffer.
    T_* buf_;   // The raw buffer.


public:

    // Constructs an empty buffer.
    Buffer():
          cap_(0)
        , buf_(nullptr)
    {}

    // Constructs a buffer with capacity `cap`.
    Buffer(const std::size_t cap):
          cap_(cap)
        , buf_(allocate<T_>(cap))
    {}

    ~Buffer() { deallocate(buf_); }

    Buffer(Buffer&& rhs) { *this = std::move(rhs); }

    Buffer& operator=(Buffer&& rhs)
    {
        cap_ = rhs.cap_;
        buf_ = rhs.buf_;
        rhs.cap_ = 0;
        rhs.buf_ = nullptr;

        return *this;
    }

    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;

    // Returns the memory region of the buffer.
    T_* data() { return buf_; }

    // Returns the memory region of the buffer.
    const T_* data() const { return buf_; }

    // Returns reference to the `idx`'th element of the buffer.
    T_& operator[](const std::size_t idx) { return buf_[idx]; }

    // Returns the `idx`'th element of the buffer.
    const T_& operator[](const std::size_t idx) const { return buf_[idx]; }

    // Returns the capacity of the buffer.
    auto capacity() const { return cap_; }

    // Ensures that the buffer have space for at least `new_cap` elements. No
    // guarantees are made for the existing elements.
    void reserve_uninit(const std::size_t new_cap) { cap_ = reserve_geometric(buf_, cap_, new_cap); }

    // Ensures that the buffer have space for at least `new_cap` elements.
    void reserve(const std::size_t new_cap) { cap_ = reserve_geometric<T_, true>(buf_, cap_, new_cap); }

    // Resizes the buffer to have capacity `cap`. No guarantees are made for
    // the existing elements.
    void resize_uninit(const std::size_t cap) { deallocate(buf_); buf_ = allocate<T_>(cap); cap_ = cap; }

    // Resizes the buffer to have capacity `cap`, initialized with `0`.
    void resize_init(const std::size_t cap) { deallocate(buf_); buf_ = allocate_zeroed<T_>(cap); cap_ = cap; }

    // Frees the buffer's memory.
    void free() { deallocate(buf_); buf_ = nullptr; cap_ = 0; }

    // Returns the resident set size of the buffer.
    std::size_t RSS() const { return sizeof(buf_) + sizeof(cap_) + capacity() * sizeof(T_); }

    // Serializes the buffer to the `cereal` archive `archive`.
    template <typename T_archive_>
    void save(T_archive_& archive) const
    {
        archive(cap_);
        for(std::size_t i = 0; i < cap_; ++i)
            archive(buf_[i]);
    }

    // Deserializes the buffer from the `cereal` archive `archive`.
    template <typename T_archive_>
    void load(T_archive_& archive)
    {
        archive(cap_);
        resize_uninit(cap_);
        for(std::size_t i = 0; i < cap_; ++i)
            archive(buf_[i]);
    }
};


namespace timer
{
    inline auto now() { return std::chrono::high_resolution_clock::now(); }

    inline auto duration(const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };
}


namespace type
{
    template <typename T_> inline T_& mut_ref(const T_& v) { return const_cast<T_&>(v); }
}



#endif
