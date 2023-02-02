
#ifndef KMER_SPMC_ITERATOR_HPP
#define KMER_SPMC_ITERATOR_HPP



#include "Kmer.hpp"
#include "Kmer_Container.hpp"
#include "kmc_api/kmc_file.h"

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>
#include <string>
#include <thread>


// Data required by the consumers to correctly parse raw binary k-mers.
struct alignas(L1_CACHE_LINE_SIZE)
    Consumer_Data
{
    uint8_t* suff_buf{nullptr}; // Buffer for the raw binary suffixes of the k-mers.
    uint64_t kmers_available;   // Number of k-mers present in the current buffer.
    uint64_t kmers_parsed;      // Number of k-mers parsed from the current buffers.
    std::vector<std::pair<uint64_t, uint64_t>> pref_buf;    // Buffer for the raw binary prefixes of the k-mers, in the form: <prefix, #corresponding_suffix>
    std::vector<std::pair<uint64_t, uint64_t>>::iterator pref_it;   // Pointer to the prefix to start parsing k-mers from.
    // uint64_t pad_[1];           // Padding to avoid false-sharing.
};

// An "iterator" class to iterate over a k-mer database on disk, where a single producer thread
// (sequentially) reads the raw binary representations of the k-mers from disk, and a number of
// different consumer threads fetch (and parse) the raw binary k-mers.
// Note: in a technical sense, it's not an iterator.
template <uint16_t k>
class Kmer_SPMC_Iterator
{
    typedef Kmer_SPMC_Iterator iterator;


private:

    const Kmer_Container<k>* const kmer_container;  // The associated k-mer container over which to iterate.
    CKMC_DB kmer_database; // The k-mer database object.

    const uint64_t kmer_count;  // Number of k-mers present in the underlying database.
    const size_t consumer_count;  // Total number of consumer threads of the iterator.

    uint64_t kmers_read;    // Number of raw k-mers read (off disk) by the iterator.

    std::unique_ptr<std::thread> reader{nullptr};   // The thread doing the actual disk-read of the binary data, i.e. the producer thread.

    static constexpr size_t BUF_SZ_PER_CONSUMER = (1 << 24);   // Size of the consumer-specific buffers (in bytes): 16 MB.

    std::vector<Consumer_Data> consumer;   // Parsing data required for each consumer.

    // Status of the tasks for each consumer thread.
    enum class Task_Status: uint8_t
    {
        pending,    // k-mers yet to be provided;
        available,  // k-mers are available and waiting to be parsed and processed;
        no_more,    // no k-mers will be provided anymore.
    };

    std::atomic<Task_Status>* task_status{nullptr}; // Collection of the task statuses of the consumers.


    // Opens the k-mer database file with the path prefix `db_path`.
    void open_kmer_database(const std::string& db_path);

    // Closes the k-mer database file.
    void close_kmer_database();

    // Reads raw binary k-mer representations from the underlying k-mer database, and
    // makes those available for consumer threads. Reading continues until the database
    // has been depleted.
    void read_raw_kmers();

    // Returns the id (number) of an idle consumer thread.
    size_t get_idle_consumer() const;


public:

    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position, based on `at_begin` and `at_end`. The
    // iterator is to support `consumer_count` number of different consumers.
    Kmer_SPMC_Iterator(const Kmer_Container<k>* kmer_container, size_t consumer_count, bool at_begin = true, bool at_end = false);

    // Copy constructs an iterator from another one `other`.
    // Note: this should be prohibited, like the `operator=`. But the BBHash code
    // requires this to be implemented.
    Kmer_SPMC_Iterator(const iterator& other);

    // Destructs the iterator.
    ~Kmer_SPMC_Iterator();

    // Prohibits assignment-copying. This is a complex object with a background disk-reader
    // thread and a KMC database object in some *arbitrary* state. These should not be allowed
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

    // Whether production has been launched yet — and more specifically, whether the production data
    // structures have been instantiated yet. This must be found true before attempting any sort of
    // access into the data structure. The only exception is the `launch_production` invokation.
    bool launched() const;

    // Waits for the disk-reads of the raw k-mers to be completed, and then waits for the consumers
    // to finish their ongoing tasks; then signals them that no more data are to be provided, and
    // also closes the k-mer database. 
    void seize_production();

    // Returns `true` iff tasks might be provided to the consumer with id `consumer_id` in future.
    bool tasks_expected(size_t consumer_id) const;

    // Returns `true` iff a task is available for the consumer with id `consumer_id`.
    bool task_available(size_t consumer_id) const;

    // Returns the memory (in bytes) used by the iterator.
    std::size_t memory() const;

    // Returns the memory (in bytes) to be used by an iterator supporting `consumer_count` consumers.
    static std::size_t memory(std::size_t consumer_count);

    // Dummy methods.
    const iterator& operator++() { return *this; }
    Kmer<k> operator*() { return Kmer<k>(); }
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
           delete[] consumer[id].suff_buf;

        std::cerr << "\nCompleted a pass over the k-mer database.\n";
    }
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::open_kmer_database(const std::string& db_path)
{
    if(!kmer_database.open_for_cuttlefish_listing(db_path))
    {
        std::cerr << "Error opening k-mer database with prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::close_kmer_database()
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
    if(launched())
        return;

    // Initialize the buffers and the parsing data structures.

    task_status = new std::atomic<Task_Status>[consumer_count];

    consumer.resize(consumer_count);
    for(size_t id = 0; id < consumer_count; ++id)
    {
        auto& consumer_state = consumer[id];
        consumer_state.suff_buf = new uint8_t[BUF_SZ_PER_CONSUMER];
        consumer_state.kmers_available = 0;
        consumer_state.kmers_parsed = 0;
        consumer_state.pref_buf.clear();
        consumer_state.pref_it = consumer_state.pref_buf.begin();
        task_status[id] = Task_Status::pending;
    }

    // Open the underlying k-mer database.
    open_kmer_database(kmer_container->container_location());

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
inline bool Kmer_SPMC_Iterator<k>::launched() const
{
    return reader != nullptr;
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::read_raw_kmers()
{
    while(!kmer_database.Eof())
    {
        const size_t consumer_id = get_idle_consumer();
        Consumer_Data& consumer_state = consumer[consumer_id];

        consumer_state.kmers_available = kmer_database.read_raw_suffixes(consumer_state.suff_buf, consumer_state.pref_buf, BUF_SZ_PER_CONSUMER);
        consumer_state.pref_it = consumer_state.pref_buf.begin();

        if(!consumer_state.kmers_available)
        {
            std::cerr << "Error reading the suffix file. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        kmers_read += consumer_state.kmers_available;

        consumer_state.kmers_parsed = 0;
        task_status[consumer_id] = Task_Status::available;
    }
}


template <uint16_t k>
inline size_t Kmer_SPMC_Iterator<k>::get_idle_consumer() const
{
    size_t id{0};

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


    // Wait for the consumers to finish consumption, and signal them that the means of production have been seized.
    for(size_t id = 0; id < consumer_count; ++id)
    {
        while(task_status[id] != Task_Status::pending); // busy-wait
        
        task_status[id] = Task_Status::no_more;
    }

    // Close the underlying k-mer database.
    close_kmer_database();
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::value_at(const size_t consumer_id, Kmer<k>& kmer)
{
    // TODO: try to delay this `volatile` access as much as possible, as each access directly hits the actual memory location.
    if(!task_available(consumer_id))
        return false;

    auto& ts = consumer[consumer_id];
    if(ts.kmers_parsed == ts.kmers_available)
    {
        task_status[consumer_id] = Task_Status::pending;
        return false;
    }

    kmer_database.parse_kmer_buf<k>(ts.pref_it, ts.suff_buf, ts.kmers_parsed * kmer_database.suff_record_size(), kmer);
    ts.kmers_parsed++;

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


template <uint16_t k>
inline std::size_t Kmer_SPMC_Iterator<k>::memory() const
{
    return CKMC_DB::pref_buf_memory() + (consumer_count * BUF_SZ_PER_CONSUMER);
}


template <uint16_t k>
inline std::size_t Kmer_SPMC_Iterator<k>::memory(const std::size_t consumer_count)
{
    return CKMC_DB::pref_buf_memory() + (consumer_count * BUF_SZ_PER_CONSUMER);
}



#endif
