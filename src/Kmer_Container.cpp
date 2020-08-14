
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"


template <uint16_t k>
Kmer_Container<k>::Kmer_Container(const std::string& kmc_file_path):
    kmc_file_path(kmc_file_path)
{
    CKMCFile kmer_database;
    if(!kmer_database.OpenForListing(kmc_file_path))
    {
        std::cout << "Error opening KMC database files with prefix " << kmc_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    if(!kmer_database.Info(kmer_database_info))
    {
        std::cout << "Error reading from KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    kmer_database.Close();
}


template <uint16_t k>
const std::string& Kmer_Container<k>::container_location() const
{
    return kmc_file_path;
}


template <uint16_t k>
uint32_t Kmer_Container<k>::kmer_length() const
{
    return kmer_database_info.kmer_length;
}


template <uint16_t k>
uint64_t Kmer_Container<k>::size() const
{
   return kmer_database_info.total_kmers;
}


template <uint16_t k>
typename Kmer_Container<k>::iterator Kmer_Container<k>::begin() const
{
    return iterator(this);
}


template <uint16_t k>
typename Kmer_Container<k>::iterator Kmer_Container<k>::end() const
{
    return iterator(this, false);
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Container)
