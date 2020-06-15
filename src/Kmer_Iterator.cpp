
#include "Kmer_Iterator.hpp"


Kmer_Iterator::Kmer_Iterator(Kmer_Container* const kmer_container, const bool at_begin):
    kmer_container(kmer_container), kmer_object(), at_begin(at_begin)
{
    if(at_begin)
    {
        open_kmer_database();
        kmer_object = CKmerAPI(kmer_container->kmer_length());
        advance();
    }
}


void Kmer_Iterator::open_kmer_database()
{
    if(!kmer_database_input.OpenForListing(kmer_container->kmc_file_name))
    {
        std::cerr << "Error opening KMC database with prefix " << kmer_container->kmc_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


void Kmer_Iterator::advance()
{
    uint32_t dummy;
    if(!kmer_database_input.ReadNextKmer(kmer_object, dummy))
    {
        std::cerr << "\n\nHERE\n\n";
        kmer_object = CKmerAPI();
        kmer_database_input.Close();
    }
    else
        kmer.from_CKmerAPI(kmer_object);// = cuttlefish::kmer_t(kmer_object);
}


Kmer_Iterator::Kmer_Iterator(const iterator& other):
    kmer_container(other.kmer_container), kmer_object(other.kmer_object), at_begin(other.at_begin)
{
    if(at_begin)
    {
        open_kmer_database();
        advance();
    }
}


const Kmer_Iterator& Kmer_Iterator::operator=(const iterator& rhs)
{
    kmer_container = rhs.kmer_container;
    kmer_object = rhs.kmer_object;
    at_begin = rhs.at_begin;

    if(at_begin)
    {
        open_kmer_database();
        advance();
    }
    std::cerr << "\n\n\ncopied iterator\n\n\n"; 

    return *this;
}


Kmer_Iterator::value_type Kmer_Iterator::operator*() const
{
    return kmer;
    // return cuttlefish::kmer_t(kmer_object);
}


const Kmer_Iterator& Kmer_Iterator::operator++()
{
    advance();
    if(at_begin)
        at_begin = false;

    return *this;
}


bool Kmer_Iterator::operator==(const iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && kmer_object == rhs.kmer_object;
}


  bool Kmer_Iterator::operator!=(const iterator& rhs) const
{
    return !(this->operator==(rhs));
}