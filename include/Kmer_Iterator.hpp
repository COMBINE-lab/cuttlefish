
#ifndef KMER_ITERATOR_HPP
#define KMER_ITERATOR_HPP



#include "globals.hpp"
#include "kmc_api/kmc_file.h"
#include "Kmer_Container.hpp"


// Iterator class to iterate over KMC databases on disk.
class Kmer_Iterator
{
    friend class Kmer_Container;


public:

    typedef Kmer_Iterator iterator;
    
    // iterator traits
    typedef std::input_iterator_tag iterator_category;
    typedef cuttlefish::kmer_t value_type;
    typedef int difference_type;
    typedef cuttlefish::kmer_t* pointer;
    typedef cuttlefish::kmer_t& reference;

    typedef const cuttlefish::kmer_t* const_ptr_t;



private:

    const Kmer_Container* kmer_container; // The associated k-mer container on which to iterate on.
    CKMCFile kmer_database_input;   // The input reader object (from KMC databases).
    CKmerAPI kmer_object;   // Current KMC k-mer object that this iterator is holding.
    cuttlefish::kmer_t kmer;    // K-mer present inside the `kmer_object` api.
    bool at_begin;  // Whether this iterator points to the beginning of the KMC database or not.


    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position based on the value of `at_begin`.
    Kmer_Iterator(const Kmer_Container* kmer_container, bool at_begin = true);

    // Opens the KMC database (internally buffered) to read k-mers.
    void open_kmer_database();

    // Advances the iterator forward by offset one.
    void advance();



public:

    // Copy constructs an iterator from the another one `other`.
    Kmer_Iterator(const iterator& other);

    // Assigns the iterator `rhs` to this one, and returns the new iterator.
    const iterator& operator=(const iterator& rhs);

    // Returns a `Kmer` object corresponding to the iterator.
    value_type operator*() const;

    // Returns a const pointer to the `Kmer` object corresponding to the iterator.
    const_ptr_t operator->() const;

    // Advances the iterator by offset one, and returns the new iterator.
    const iterator& operator++();

    // Advances the iterator by offset one, and returns the old iterator.
    // TODO
    iterator operator++(int);

    // Returns true iff this and `rhs` -- both the iterators refer to the same
    // container and have the same KMC k-mer object.
    bool operator==(const iterator& rhs) const;

    // Returns true iff the iterators this and `rhs` -- either they refer to
    // different containers, or have different KMC k-mer objects.
    bool operator!=(const iterator& rhs) const;
};



inline Kmer_Iterator::Kmer_Iterator(const Kmer_Container* const kmer_container, const bool at_begin):
    kmer_container(kmer_container), kmer_object(), at_begin(at_begin)
{
    if(at_begin)
    {
        open_kmer_database();
        kmer_object = CKmerAPI(kmer_container->kmer_length());
        advance();
    }
}


inline void Kmer_Iterator::open_kmer_database()
{
    if(!kmer_database_input.OpenForListing(kmer_container->kmc_file_name))
    {
        std::cerr << "Error opening KMC database with prefix " << kmer_container->kmc_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline void Kmer_Iterator::advance()
{
    uint32_t count;
    if(!kmer_database_input.ReadNextKmer(kmer_object, count))
    {
        kmer_object = CKmerAPI();
        kmer_database_input.Close();
    }
    else
        kmer.from_CKmerAPI(kmer_object);// = cuttlefish::kmer_t(kmer_object);
}


inline Kmer_Iterator::Kmer_Iterator(const iterator& other):
    kmer_container(other.kmer_container), kmer_object(other.kmer_object), at_begin(other.at_begin)
{
    if(at_begin)
    {
        open_kmer_database();
        advance();
    }
}


inline const Kmer_Iterator& Kmer_Iterator::operator=(const iterator& rhs)
{
    kmer_container = rhs.kmer_container;
    kmer_object = rhs.kmer_object;
    at_begin = rhs.at_begin;

    if(at_begin)
    {
        open_kmer_database();
        advance();
    }

    return *this;
}


inline Kmer_Iterator::value_type Kmer_Iterator::operator*() const
{
    // return cuttlefish::kmer_t(kmer_object);
    return kmer;
}


inline Kmer_Iterator::const_ptr_t Kmer_Iterator::operator->() const
{
    return &kmer;
}


inline const Kmer_Iterator& Kmer_Iterator::operator++()
{
    advance();
    if(at_begin)
        at_begin = false;

    return *this;
}


inline bool Kmer_Iterator::operator==(const iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && kmer_object == rhs.kmer_object;
}


inline bool Kmer_Iterator::operator!=(const iterator& rhs) const
{
    return !(this->operator==(rhs));
}



#endif
