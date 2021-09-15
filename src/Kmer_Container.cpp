
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"

template <uint16_t k>
Kmer_Container<k>::Kmer_Container(const std::string& kmc_file_path):
    kmc_file_path(kmc_file_path)
{
    CKMC_DB kmer_database;
    if(!kmer_database.read_parameters(kmc_file_path))
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


    const uint16_t kmer_len = kmer_length();
    if(kmer_len != k)
    {
        std::cerr << "Expected k value " << k << ", but is provided with a " << kmer_len << "-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
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
uint64_t Kmer_Container<k>::size(const std::string& kmc_db_path)
{
    const Kmer_Container<k> kmer_container(kmc_db_path);
    return kmer_container.size();
}


template <uint16_t k>
void Kmer_Container<k>::remove(const std::string& kmc_db_path)
{
    const std::string kmc_pref_file(kmc_db_path + ".kmc_pre");
    const std::string kmc_suff_file(kmc_db_path + ".kmc_suf");

    if(!remove_file(kmc_pref_file) || !remove_file(kmc_suff_file))
    {
        std::cerr << "Error removing the KMC database file from path prefix " << kmc_db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


// template <uint16_t k>
// typename Kmer_Container<k>::iterator Kmer_Container<k>::end() const
// {
//     return iterator(this, false);
// }

// template <uint16_t k>
// typename Kmer_Container<k>::buf_iterator Kmer_Container<k>::buf_begin() const
// {
//     return buf_iterator(this, true, false);
// }


// template <uint16_t k>
// typename Kmer_Container<k>::buf_iterator Kmer_Container<k>::buf_end() const
// {
//     return buf_iterator(this, false, true);
// }


template <uint16_t k>
typename Kmer_Container<k>::spmc_iterator Kmer_Container<k>::spmc_begin(const size_t consumer_count) const
{
    return spmc_iterator(this, consumer_count);   
}


template <uint16_t k>
typename Kmer_Container<k>::spmc_iterator Kmer_Container<k>::spmc_end(const size_t consumer_count) const
{
    return spmc_iterator(this, consumer_count, false, true);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_ALL, Kmer_Container)
