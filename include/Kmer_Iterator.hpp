
#ifndef KMER_ITERATOR_HPP
#define KMER_ITERATOR_HPP



#include "kmc_api/kmc_file.h"
#include "Kmer.hpp"
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



private:

    Kmer_Container* kmer_container; // The associated k-mer container on which to iterate on.
    CKmerAPI kmer_object;   // Current KMC k-mer object that this iterator is holding.


    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position based on the value of `at_begin`.
    Kmer_Iterator(Kmer_Container* const kmer_container, const bool at_begin = true);

    // Advances the iterator forward by offset one.
    void advance();



public:

    // Copy constructs an iterator from the another one `other`.
    Kmer_Iterator(const iterator& other);

    // Returns a `Kmer` object corresponding to the current iterator position.
    value_type operator*() const;

    // Advances the iterator by offset one, and returns the new iterator.
    const iterator& operator++();

    // Advances the iterator by offset one, and returns the old iterator.
    iterator operator++(int);

    // Assigns the iterator `rhs` to this one, and returns the new iterator.
    const iterator& operator=(const iterator& rhs);

    // Returns true iff this and `rhs` -- both the iterators refer to the same
    // container and have the same KMC k-mer object.
    bool operator==(const iterator& rhs) const;

    // Returns true iff the iterators this and `rhs` -- either they refer to
    // different containers, or have different KMC k-mer objects.
    bool operator!=(const iterator& rhs) const;
};



inline Kmer_Iterator::Kmer_Iterator(Kmer_Container* const kmer_container, const bool at_begin):
    kmer_container(kmer_container), kmer_object()
{
    if(at_begin)
    {
        kmer_object = CKmerAPI(kmer_container->kmer_length());
        advance();
    }
}


inline void Kmer_Iterator::advance()
{
    uint32_t dummy;
    if(!kmer_container->kmer_database.ReadNextKmer(kmer_object, dummy))
        kmer_object = CKmerAPI();
}


inline Kmer_Iterator::Kmer_Iterator(const iterator& other):
    kmer_container(other.kmer_container), kmer_object(other.kmer_object)
{}


inline Kmer_Iterator::value_type Kmer_Iterator::operator*() const
{
    return cuttlefish::kmer_t(kmer_object);
}


inline const Kmer_Iterator& Kmer_Iterator::operator++()
{
    advance();

    return *this;
}


inline Kmer_Iterator Kmer_Iterator::operator++(int)
{
    Kmer_Iterator clone(*this);
    advance();

    return clone;
}


inline const Kmer_Iterator& Kmer_Iterator::operator=(const iterator& rhs)
{
    kmer_container = rhs.kmer_container;
    kmer_object = rhs.kmer_object;

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
