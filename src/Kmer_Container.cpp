
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"


Kmer_Container::Kmer_Container(const std::string& kmc_file_name):
    kmc_file_name(kmc_file_name)
{
    CKMCFile kmer_database;
    if(!kmer_database.OpenForListing(kmc_file_name))
    {
        std::cout << "Error opening KMC database files with prefix " << kmc_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    if(!kmer_database.Info(kmer_database_info))
    {
        std::cout << "Error reading from KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    kmer_database.Close();
}


std::string Kmer_Container::container_location() const
{
    return kmc_file_name;
}


uint32_t Kmer_Container::kmer_length() const
{
    return kmer_database_info.kmer_length;
}


uint64_t Kmer_Container::size() const
{
   return kmer_database_info.total_kmers;
}


Kmer_Container::iterator Kmer_Container::begin() const
{
    return iterator(this);
}


Kmer_Container::iterator Kmer_Container::end() const
{
    return iterator(this, false);
}
