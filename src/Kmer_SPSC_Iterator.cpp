
#include "Kmer_SPSC_Iterator.hpp"
#include "Kmer_Container.hpp"
#include "globals.hpp"


template <uint16_t k>
Kmer_SPSC_Iterator<k>::Kmer_SPSC_Iterator()
{}


template <uint16_t k>
Kmer_SPSC_Iterator<k>::Kmer_SPSC_Iterator(const std::string& kmer_db_path):
    kmer_db_path{kmer_db_path},
    kmer_count{Kmer_Container<k>(kmer_db_path).size()},
    kmers_read{0},
    suff_buf{nullptr},
    kmers_parsed_off_buf{0},
    kmers_available_in_buf{0},
    buf_stat{Buffer_Status::pending}
{}


template <uint16_t k>
Kmer_SPSC_Iterator<k>::Kmer_SPSC_Iterator(const Kmer_Container<k>* const kmer_container):
    kmer_db_path{kmer_container->container_location()},
    kmer_count{kmer_container->size()},
    kmers_read{0},
    suff_buf{nullptr},
    kmers_parsed_off_buf{0},
    kmers_available_in_buf{0},
    buf_stat{Buffer_Status::pending}
{}


template <uint16_t k>
Kmer_SPSC_Iterator<k>::~Kmer_SPSC_Iterator()
{
    if(suff_buf != nullptr)
        delete[] suff_buf;

    pref_buf.clear();
    pref_buf.shrink_to_fit();
}


template <uint16_t k>
void Kmer_SPSC_Iterator<k>::init(const std::string& kmer_db_path)
{
    this->kmer_db_path = kmer_db_path;
    kmer_count = Kmer_Container<k>(kmer_db_path).size();
    kmers_read = 0;
    suff_buf = nullptr;
    kmers_parsed_off_buf = 0;
    kmers_available_in_buf = 0;
    buf_stat = Buffer_Status::pending;
}


template <uint16_t k>
void Kmer_SPSC_Iterator<k>::open_kmer_database(const std::string& db_path)
{
    if(!kmer_database.open_for_cuttlefish_listing(db_path))
    {
        std::cerr << "\nError opening k-mer database with path prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
void Kmer_SPSC_Iterator<k>::close_kmer_database()
{
    if(!kmer_database.Close())
    {
        std::cerr << "\nError closing k-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
void Kmer_SPSC_Iterator<k>::launch()
{
    open_kmer_database(kmer_db_path);

    suff_buf = new uint8_t[suf_buf_size];
    pref_buf.clear();
    set_buffer_status(Buffer_Status::pending);
}


template <uint16_t k>
void Kmer_SPSC_Iterator<k>::close()
{
    close_kmer_database();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_SPSC_Iterator)
