
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"


Kmer_Container::Kmer_Container(const std::string& kmc_file_name)
{
    if(!kmer_database.OpenForListing(kmc_file_name))
    {
        std::cout << "Error opening KMC database files with prefix " << kmc_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


uint32_t Kmer_Container::kmer_length() const
{
    CKMCFileInfo database_info;
    if(!kmer_database.Info(database_info))
    {
        std::cout << "Error reading from KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return database_info.kmer_length;
}


uint64_t Kmer_Container::size() const
{
    CKMCFileInfo database_info;
    if(!kmer_database.Info(database_info))
    {
        std::cout << "Error reading from KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    return database_info.total_kmers;
}


Kmer_Iterator Kmer_Container::begin()
{
    return Kmer_Iterator(this);
}


Kmer_Iterator Kmer_Container::end()
{
    return Kmer_Iterator(this, false);
}


void Kmer_Container::close()
{
    if(!kmer_database.Close())
    {
        std::cerr << "Error closing the KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}
