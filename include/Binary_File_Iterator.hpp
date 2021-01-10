
#ifndef BINARY_FILE_ITERATOR_HPP
#define BINARY_FILE_ITERATOR_HPP



#include "Compressed_Stream.hpp"


template <typename T_elem_>
class Binary_File_Iterator
{
    typedef Binary_File_Iterator iterator;

private:

    const std::string file_path;    //
    T_elem_ elem;   //
    std::unique_ptr<Compressed_Stream<T_elem_>> stream{nullptr}; //
    size_t offset{0};   // 
    T_elem_* buf{nullptr};   //
    std::streamsize elems_in_buf{0};  //
    std::streamsize elems_read{0};   //
    static constexpr size_t BUF_SZ = 100000;    //
    static constexpr size_t END_OFF = std::numeric_limits<size_t>::max();


    // 
    void peek();

    // 
    void advance();

    
public:

    // 
    Binary_File_Iterator()
    {}

    // Should be prohibited.
    Binary_File_Iterator(const iterator& other);

    // 
    Binary_File_Iterator(const std::string& file_path, const bool at_end = false);

    // 
    ~Binary_File_Iterator();

    // 
    const T_elem_& operator*();

    // 
    void operator++();

    // 
    bool operator==(const iterator& rhs) const;

    // 
    bool operator!=(const iterator& rhs) const;

    // 
    iterator& operator=(const iterator&) = delete;

    // Dummy methods.
    void launch_production() {}
    bool launched() { return true; }
    bool value_at(const size_t consumer_id, T_elem_& elem) { (void)consumer_id; (void)elem; return false; }
    bool tasks_expected(const size_t consumer_id) const { (void)consumer_id; return false; }
    void seize_production() {}
};


template <typename T_elem_>
inline Binary_File_Iterator<T_elem_>::Binary_File_Iterator(const iterator& other):
    file_path(other.file_path),
    elem(other.elem),
    offset{other.offset}
{}


template <typename T_elem_>
inline Binary_File_Iterator<T_elem_>::Binary_File_Iterator(const std::string& file_path, const bool at_end):
    file_path(file_path),
    offset{at_end ? END_OFF : 0}
{
    peek();
}


template <typename T_elem_>
inline Binary_File_Iterator<T_elem_>::~Binary_File_Iterator()
{
    if(buf != nullptr)
        delete[] buf;
}


template <typename T_elem_>
inline void Binary_File_Iterator<T_elem_>::peek()
{
    FILE* fp = std::fopen(file_path.c_str(), "rb");

    if(!fp)
    {
        std::cerr << "Error opening working file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    const int ch = std::fgetc(fp);
    std::ungetc(ch, fp);

    if(ch == EOF)
    {
        stream = nullptr;
        offset = END_OFF;
    }

    std::fclose(fp);
}


template <typename T_elem_>
inline void Binary_File_Iterator<T_elem_>::advance()
{
    if(stream == nullptr)
    {
        stream = std::unique_ptr<Compressed_Stream<T_elem_>>(new Compressed_Stream<T_elem_>(file_path, false));
        buf = new T_elem_[BUF_SZ];
    }


    offset++;

    if(elems_read == elems_in_buf)
    {
        elems_in_buf = stream->read(buf, BUF_SZ);
        elems_read = 0;

        if(elems_in_buf == 0)
        {
            stream = nullptr;
            offset = END_OFF;
            return;
        }
    }

    elem = buf[elems_read++];
}


template <typename T_elem_>
inline const T_elem_& Binary_File_Iterator<T_elem_>::operator*()
{
    if(stream == nullptr)
        advance();
    
    return elem;
}


template <typename T_elem_>
inline void Binary_File_Iterator<T_elem_>::operator++()
{
    advance();
}


template <typename T_elem_>
inline bool Binary_File_Iterator<T_elem_>::operator==(const iterator& rhs) const
    {
        return stream == rhs.stream && offset == rhs.offset;
    }

template <typename T_elem_>
inline bool Binary_File_Iterator<T_elem_>::operator!=(const iterator& rhs) const
{
    return !operator==(rhs);
}



#endif
