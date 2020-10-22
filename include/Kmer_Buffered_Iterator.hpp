
#ifndef KMER_BUFFERED_ITERATOR_HPP
#define KMER_BUFFERED_ITERATOR_HPP



#include "globals.hpp"
#include "Kmer_Container.hpp"
#include "kmc_api/kmc_file.h"
#include "readerwriterqueue/readerwriterqueue.h"

#include <thread>


// Branch prediction hints
#define OUR_LIKELY(x) __builtin_expect(x, 1)
#define OUR_UNLIKELY(x) __builtin_expect(x, 0)


// Wrapper class for a collection of k-mers (`Kmer<k>`).
template <uint16_t k>
class Kmer_Chunk
{
private:

    std::vector<Kmer<k>> group_;  // The collection of k-mers.
    const size_t capacity_; // Maximum number of k-mers this chunk can contain.
    size_t size_; // Number of k-mers this chunk currently contains.


public:

    // Constructs a `Kmer_Chunk` with the capacity to contain `capacity` k-mers.
    Kmer_Chunk(const size_t capacity) : group_(capacity), capacity_(capacity), size_(0) {}

    // Sets the number of k-mers in the chunk to `chunk_size`.
    void set_size(const size_t chunk_size) { size_ = chunk_size; }

    // Returns the number of k-mers in the chunk.
    size_t size() const { return size_; }

    // Returns the capacity of the chunk.
    size_t capacity() const { return capacity_; }

    // Returns a reference to the `idx`'th k-mer from the chunk.
    Kmer<k>& operator[](const size_t idx) { return group_[idx]; }

    // Returns a reference to the `idx`'th k-mer from the chunk.
    Kmer<k>& at(const size_t idx) { return this->operator[](idx); }

    // Returns an iterator pointing to the beginning of the chunk.
    typename std::vector<Kmer<k>>::iterator begin() { return group_.begin(); }

    // Returns an iterator pointing to the end of the chunk.
    typename std::vector<Kmer<k>>::iterator end() { return group_.begin() + size_; }
};



// Iterator class to iterate over KMC databases on disk, in a buffered manner.
template <uint16_t k>
class Kmer_Buffered_Iterator
{
    friend class Kmer_Container<k>;


public:

    typedef Kmer_Buffered_Iterator iterator;
    
    // iterator traits
    typedef std::input_iterator_tag iterator_category;
    typedef Kmer<k> value_type;
    typedef int difference_type;
    typedef Kmer<k>* pointer;
    typedef Kmer<k>& reference;

    // Can't support the `->` operator, as the contained k-mer isn't read off disk until after the first dereferencing.
    // typedef const Kmer<k>* const_ptr_t;


private:

    const Kmer_Container<k>* const kmer_container;    // The associated k-mer container on which to iterate on.
    const uint64_t kmer_count;  // Number of k-mers present at the underlying database.

    bool at_begin;  // Whether this iterator points to the beginning of the KMC database or not.
    bool at_end;  // Whether this iterator points to the ending of the KMC database or not.
    uint64_t offset{0};   // Offset position of the iterator relative to the beginning of the database.
    Kmer<k> kmer;   // The k-mer associated to the current iterator position.

    std::unique_ptr<std::thread> parser{nullptr};   // The thread doing the actual disk-reads.
    moodycamel::ReaderWriterQueue<Kmer_Chunk<k>*> rwq;   // A single-producer, single-consumer queue containing k-mer chunks read off the database.

    Kmer_Chunk<k>* curr_chunk{nullptr}; // Current k-mer chunk (parsed from the database) to read off k-mers from.
    typename std::vector<Kmer<k>>::iterator cur_chunk_it;   // Iterator pointing to the k-mer last read from the current chunk `curr_chunk`.

    static constexpr size_t QUEUE_SIZE = (k < 32 ? 100 : 500);   // Maximum number of k-mer chunks to contain at the chunk-queue `rwq`.
    static constexpr size_t CHUNK_CAPACITY = 40000; // Maximum number of k-mers to put in each k-mer chunk.


    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position based on the values of `at_begin` and `at_end`.
    Kmer_Buffered_Iterator(const Kmer_Container<k>* kmer_container, bool at_begin = true, bool at_end = false);

    // Opens the KMC database at path `db_path` into `kmer_database` for reading.
    static void open_kmer_database(const std::string& db_path, CKMCFile& kmer_database);

    // Closes the KMC database `kmer_database`.
    static void close_kmer_database(CKMCFile& kmer_database);

    // Launches background parsing of k-mers from the database, and returns the first parsed k-mer.
    value_type launch_parse();

    // Parses k-mers from the underlying KMC database `kmer_database` into k-mer chunks and puts the chunks
    // into the read-write queue `rwq`.
    void parse_kmers();

    // Pushes the k-mer chunk `kc` into the read-write queue `rwq` setting its size to `chunk_size`.
    void dump_chunk(Kmer_Chunk<k>* kc, size_t chunk_size);

    // Advances the iterator forward by offset one.
    // Note that, a precondition to invoke `advance` is that `curr_chunk` is pointing to a valid chunk,
    // i.e., parsing has already started and a chunk has been fetched from `rwq` to `curr_chunk`.
    void advance();


public:

    // Copy constructs an iterator from the another one `other`.
    Kmer_Buffered_Iterator(const iterator& other);

    // Destructs the iterator.
    ~Kmer_Buffered_Iterator();

    // Prohibits assignment-copying, as an iterator contains a KMC database object at some state and a
    // read-write queue of k-mer chunks. We shouldn't allow copying of these objects at arbitrary states.
    iterator& operator=(const iterator& rhs) = delete;

    // Returns the k-mer associated to the iterator position.
    value_type operator*() ;

    // Advances the iterator by offset one, and returns the new iterator. The returned iterator is not
    // actually usable as `operator=` has been explicitly prohibited for the class.
    const iterator& operator++();

    // Returns true iff this and `rhs` -- both the iterators refer to the same container and to the same
    // position of the underlying databse.
    bool operator==(const iterator& rhs) const;

    // Returns true iff the iterators this and `rhs` -- either they refer to different containers, or are
    // referring to different positions of the underlying database.
    bool operator!=(const iterator& rhs) const;
};


template <uint16_t k>
inline Kmer_Buffered_Iterator<k>::Kmer_Buffered_Iterator(const Kmer_Container<k>* const kmer_container, const bool at_begin, const bool at_end):
    kmer_container(kmer_container), kmer_count(kmer_container->size()),
    at_begin(at_begin), at_end(at_end), offset(at_begin ? 0 : kmer_count)
{
    if(!(at_begin ^ at_end))
    {
        std::cerr << "Invalid positions encountered during buffered iterator construction. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline Kmer_Buffered_Iterator<k>::Kmer_Buffered_Iterator(const iterator& other):
    kmer_container(other.kmer_container), kmer_count(other.kmer_count),
    at_begin(other.at_begin), at_end(other.at_end), offset(other.offset), kmer(other.kmer)
{}


template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::open_kmer_database(const std::string& db_path, CKMCFile& kmer_database)
{
    if(!kmer_database.OpenForListing(db_path))
    {
        std::cerr << "Error opening KMC database with prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::close_kmer_database(CKMCFile& kmer_database)
{
    if(!kmer_database.Close())
    {
        std::cerr << "Error closing KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
void Kmer_Buffered_Iterator<k>::parse_kmers()
{
    CKmerAPI kmer_api(k);
    Kmer_Chunk<k>* kc = new Kmer_Chunk<k>(CHUNK_CAPACITY);
    size_t chunk_size{0};

    CKMCFile kmer_database;
    open_kmer_database(kmer_container->container_location(), kmer_database);

    rwq = moodycamel::ReaderWriterQueue<Kmer_Chunk<k>*>(QUEUE_SIZE);


    while(kmer_database.ReadNextKmer(kmer_api))
    {
        if(chunk_size < CHUNK_CAPACITY)
            kc->at(chunk_size).from_CKmerAPI(kmer_api);
        else
        {
            dump_chunk(kc, chunk_size);

            // Get a new chunk of k-mers and put the latest-read k-mer there.
            kc = new Kmer_Chunk<k>(CHUNK_CAPACITY);
            chunk_size = 0;
            kc->at(chunk_size).from_CKmerAPI(kmer_api);
        }

        chunk_size++;
    }


    // If there is a remaining chunk, put it into the queue.
    if(chunk_size > 0)
        dump_chunk(kc, chunk_size);


    close_kmer_database(kmer_database);
}


template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::dump_chunk(Kmer_Chunk<k>* kc, size_t chunk_size)
{
    kc->set_size(chunk_size);
    while(!rwq.try_enqueue(kc));    // busy-wait
}


template <uint16_t k>
typename Kmer_Buffered_Iterator<k>::value_type Kmer_Buffered_Iterator<k>::launch_parse()
{
    // Launch the background parsing thread.
    parser.reset(
        new std::thread([this]()
        {
            parse_kmers();
        })
    );


    // Dequeue the first chunk to avoid a null check at dereferencing in advance.
    while(!rwq.try_dequeue(curr_chunk) || curr_chunk == nullptr);   // busy-wait

    cur_chunk_it = curr_chunk->begin();
    kmer = *cur_chunk_it;
    
    return kmer;
}


template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::advance()
{
    if(OUR_UNLIKELY(++cur_chunk_it == curr_chunk->end()))
    {
        // We've exhausted the current chunk; free its memory.
        delete curr_chunk;
        curr_chunk = nullptr;
        
        if(offset < kmer_count) // k-mer chunks are available.
        {
            while(!rwq.try_dequeue(curr_chunk) || curr_chunk == nullptr);    // busy-wait

            // Update the last-read k-mer iterator and read off the k-mer.
            cur_chunk_it = curr_chunk->begin();
            kmer = *cur_chunk_it;
        }
        else    // Parsing is completed.
            at_end = true;
    }
    else
        kmer = *cur_chunk_it;   // Read off the k-mer.
}


template <uint16_t k>
inline typename Kmer_Buffered_Iterator<k>::value_type Kmer_Buffered_Iterator<k>::operator*() 
{
    return parser != nullptr ? kmer : launch_parse();
}


template <uint16_t k>
inline const Kmer_Buffered_Iterator<k>& Kmer_Buffered_Iterator<k>::operator++()
{
    if(!at_end)
    {
        offset++;
        advance();

        at_begin = false;
    }

    return *this;
}


template <uint16_t k>
inline bool Kmer_Buffered_Iterator<k>::operator==(const iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && offset == rhs.offset;
}


template <uint16_t k>
inline bool Kmer_Buffered_Iterator<k>::operator!=(const iterator& rhs) const
{
    return !(this->operator==(rhs));
}


template <uint16_t k>
inline Kmer_Buffered_Iterator<k>::~Kmer_Buffered_Iterator()
{
    if(parser != nullptr)
    {
        if(!parser->joinable())
        {
            std::cerr << "Early termination encountered for the database parser thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        parser->join();

        std::cerr << "\nCompleted a pass over the k-mer database.\n";
    }
}



#endif
