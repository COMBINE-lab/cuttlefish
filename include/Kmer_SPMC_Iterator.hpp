
#ifndef KMER_SPMC_ITERATOR_HPP
#define KMER_SPMC_ITERATOR_HPP



#include "globals.hpp"
#include "Kmer.hpp"
#include "Kmer_Container.hpp"
#include "kmc_api/kmc_file.h"

#include <vector>
#include <thread>
#include <iostream> // TODO: remove after testing.


// An "iterator" class to iterate over a k-mer database on disk, where a single producer thread
// (sequentially) reads the raw binary representations of the k-mers from disk, and a number of
// different consumer threads fetch (and parse) the raw binary k-mers.
template <uint16_t k>
class Kmer_SPMC_Iterator
{
    friend class Kmer_Container<k>;


public:

    typedef Kmer_SPMC_Iterator iterator;

    // Iterator traits.
    // typedef std::input_iterator_tag iterator_category;
    // typedef Kmer<k> value_type;
    // typedef int difference_type;
    // typedef Kmer<k>* pointer;
    // typedef Kmer<k>& reference;

    // Will `->` be supported?
    // typedef const Kmer<k>* const_ptr_t;


private:

    const Kmer_Container<k>* const kmer_container;  // The associated k-mer container over which to iterate.
    CKMCFile kmer_database; // The k-mer database object.

    const uint64_t kmer_count;  // Number of k-mers present in the underlying database.
    const size_t consumer_count;  // Total number of consumers.

    uint64_t kmers_read;    // Number of raw k-mers read (off disk) by the iterator.

    std::unique_ptr<std::thread> reader{nullptr};   // The thread doing the actual disk-read (of the binary data).

    static constexpr size_t BUF_SZ_PER_CONSUMER = (1 << 24);   // Size of the consumer-specific buffers (in bytes).
    std::vector<uint8_t*> buffer; // Buffers for the raw binary k-mers, for the consumer threads.

    // Status of the tasks for each consumer thread.
    enum class Task_Status: uint8_t
    {
        pending,    // k-mers yet to be provided;
        available,  // k-mers are available and waiting to be parsed and processed;
        no_more,    // no k-mers will be provided anymore.
    };

    // TODO: replace the raw pointer with a vector maybe?
    volatile Task_Status* task_status{nullptr}; // Collection of the task statuses of the consumers.
    
    // Index of the prefix (into the in-memory KMC prefix buffer) to start parsing (and using) k-mers from, for a consumer.
    std::vector<uint64_t> pref_idx;
    
    // Index of the suffix (into the in-disk KMC suffix collection) to start parsing (and using) k-mers from, for a consumer.
    std::vector<uint64_t> suff_idx;

    // Number of k-mers read-off disk and made available for a consumer; i.e. number of raw suffixes present in the current buffer.
    std::vector<uint64_t> kmers_available;
    
    // Number of k-mers parsed by a consumer from its current buffer.
    std::vector<uint64_t> kmers_parsed;


    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position, based on `at_begin` and `at_end`. The
    // iterator is to support `consumer_count` number of different consumers.
    Kmer_SPMC_Iterator(const Kmer_Container<k>* kmer_container, size_t consumer_count, bool at_begin = true, bool at_end = false);

    // Opens the k-mer database with the path prefix `db_path` into `kmer_database`.
    static void open_kmer_database(const std::string& db_path, CKMCFile& kmer_database);

    // Closes the k-mer database `kmer_database`.
    static void close_kmer_database(CKMCFile& kmer_database);

    // Reads raw binary k-mer representations from the underlying k-mer database, and
    // makes those available for consumer threads. Reading continues until the database
    // has been depleted.
    void read_raw_kmers();

    // Returns the id (number) of an idle consumer thread.
    size_t get_idle_consumer() const;


public:

    // Copy constructs an iterator from another one `other`.
    // Note: probably this shouldn't exist, like the `operator=`. But the BBHash code
    // requires this to be implemented.
    Kmer_SPMC_Iterator(const iterator& other);

    // Destructs the iterator.
    ~Kmer_SPMC_Iterator();

    // Prohibits assignment-copying. This is a complex object with a background disk-reader
    // thread and a KMC database object in some arbitrary state. These should not be allowed
    // to be copied. Besides, there is a number of constant fields for the iterator.
    iterator& operator=(const iterator& rhs) = delete;

    // Tries to fetch and parse the next k-mer for the consumer with id `consumer_id` into `kmer`.
    // Returns `true` iff it's successful, i.e. k-mers were remaining for this consumer.
    bool value_at(size_t consumer_id, Kmer<k>& kmer);

    // Returns `true` iff this and `rhs` — both the iterators refer to the same container and
    // the same number of raw k-mers have been read (from disk) for both.
    bool operator==(const iterator& rhs) const;

    // Returns `true` iff the iterators, this and `rhs` — either they refer to different containers,
    // or a different number of raw k-mers have been read (from disk) for them.
    bool operator!=(const iterator& rhs) const;

    // Launches the background disk-read of raw binary k-mers.
    void launch_production();

    // Waits for the disk-reads of the raw k-mers to be completed, and then waits for the consumers
    // to finish their ongoing tasks; then signals them that no more data are to be provided, and
    // also closes the k-mer database. 
    void seize_production();

    // Returns `true` iff tasks might be provided to the consumer with id `consumer_id` in future.
    bool tasks_expected(size_t consumer_id) const;

    // Returns `true` iff a task is available for the consumer with id `consumer_id`.
    bool task_available(size_t consumer_id) const;
};


template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::Kmer_SPMC_Iterator(const Kmer_Container<k>* const kmer_container, const size_t consumer_count, const bool at_begin, const bool at_end):
    kmer_container(kmer_container),
    kmer_count{kmer_container->size()},
    consumer_count{consumer_count},
    kmers_read{at_end ? kmer_count : 0}
{
    if(!(at_begin ^ at_end))
    {
        std::cerr << "Invalid position provided for SPMC k-mer iterator construction. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::Kmer_SPMC_Iterator(const iterator& other):
    kmer_container(other.kmer_container),
    kmer_count{other.kmer_count},
    consumer_count{other.consumer_count},
    kmers_read{other.kmers_read}
{}


template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::~Kmer_SPMC_Iterator()
{
    if(task_status != nullptr)
    {
        delete[] task_status;

        for(size_t id = 0; id < consumer_count; ++id)
            delete[] buffer[id];
    }
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::open_kmer_database(const std::string& db_path, CKMCFile& kmer_database)
{
    if(!kmer_database.open_for_listing_unbuffered(db_path))
    {
        std::cerr << "Error opening k-mer database with prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::close_kmer_database(CKMCFile& kmer_database)
{
    if(!kmer_database.Close())
    {
        std::cerr << "Error closing k-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::launch_production()
{
    // Initialize the buffers and the parsing data structures.

    task_status = new volatile Task_Status[consumer_count];

    buffer.resize(consumer_count, nullptr);
    for(size_t id = 0; id < consumer_count; ++id)
    {
        buffer[id] = new uint8_t[BUF_SZ_PER_CONSUMER];
        task_status[id] = Task_Status::pending;
    }

    pref_idx.resize(consumer_count);
    suff_idx.resize(consumer_count);
    kmers_available.resize(consumer_count, 0);
    kmers_parsed.resize(consumer_count, 0);


    // Open the underlying k-mer database.
    open_kmer_database(kmer_container->container_location(), kmer_database);


    // Launch the background disk-reader thread.
    reader.reset(
        new std::thread([this]()
            {
                read_raw_kmers();
            }
        )
    );
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::read_raw_kmers()
{
    std::cout << "Opened k-mer database with k-mer count: " << kmer_count << "\n";

    while(!kmer_database.Eof())
    {
        const size_t consumer_id = get_idle_consumer();
        std::cout << "Idle consumer ID: " << consumer_id << "\n";

        pref_idx[consumer_id] = kmer_database.curr_prefix_idx();
        suff_idx[consumer_id] = kmer_database.curr_suffix_idx();

        kmers_available[consumer_id] = kmer_database.read_raw_suffixes(buffer[consumer_id], BUF_SZ_PER_CONSUMER);
        if(!kmers_available[consumer_id])
        {
            std::cerr << "Error reading the suffix file. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        kmers_read += kmers_available[consumer_id];

        kmers_parsed[consumer_id] = 0;
        task_status[consumer_id] = Task_Status::available;
    }
}


template <uint16_t k>
inline size_t Kmer_SPMC_Iterator<k>::get_idle_consumer() const
{
    static size_t id{0};

    std::cout << "Looking for an idle consumer, starting from " << id << "\n";
    while(task_status[id] != Task_Status::pending)  // busy-wait
    {
        // id = (id + 1) % consumer_count;
        id++;
        if(id == consumer_count)
            id = 0;
    }

    return id;     
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::seize_production()
{
    // Wait for the disk-reads to be completed.
    if(!reader->joinable())
    {
        std::cerr << "Early termination encountered for the database reader thread. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    reader->join();

    std::cout << "Done production. Waiting for consumers to finish.\n";

    // Wait for the consumers to finish consumption, and signal them that the means of production have been seized.
    for(size_t id = 0; id < consumer_count; ++id)
    {
        while(task_status[id] != Task_Status::pending); // busy-wait
        
        task_status[id] = Task_Status::no_more;
    }


    // Close the underlying k-mer database.
    close_kmer_database(kmer_database);
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::value_at(const size_t consumer_id, Kmer<k>& kmer)
{
    if(kmers_parsed[consumer_id] == kmers_available[consumer_id])
    {
        std::cout << "Consumer " << consumer_id << " finished a chunk of size " << kmers_available[consumer_id] << "\n";
        task_status[consumer_id] = Task_Status::pending;
        return false;
    }

    kmer_database.parse_kmer<k>(pref_idx[consumer_id], suff_idx[consumer_id], buffer[consumer_id],
                                kmers_parsed[consumer_id] * kmer_database.suff_record_size(), kmer);
    kmers_parsed[consumer_id]++;

    return true;
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::operator==(const iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && kmers_read == rhs.kmers_read;
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::operator!=(const iterator& rhs) const
{
    return !operator==(rhs);
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::tasks_expected(const size_t consumer_id) const
{
    return task_status[consumer_id] != Task_Status::no_more;
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::task_available(const size_t consumer_id) const
{
    return task_status[consumer_id] == Task_Status::available;
}



#endif
