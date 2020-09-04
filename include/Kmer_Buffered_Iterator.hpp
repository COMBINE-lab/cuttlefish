
#ifndef KMER_BUFFERED_ITERATOR_HPP
#define KMER_BUFFERED_ITERATOR_HPP

#include <thread>
#include <memory>
#include "globals.hpp"
#include "Kmer_Container.hpp"
#include "kmc_api/kmc_file.h"
#include "FastxParserThreadUtils.hpp"
#include "readerwriterqueue.h"

template <uint16_t T> class KmerChunk{
public:
  KmerChunk(size_t want) : group_(want), want_(want), have_(want) {}
  inline void have(size_t num) { have_ = num; }
  inline size_t size() { return have_; }
  inline size_t want() const { return want_; }
  Kmer<T>& operator[](size_t i) { return group_[i]; }
  typename std::vector<Kmer<T>>::iterator begin() { return group_.begin(); }
  typename std::vector<Kmer<T>>::iterator end() { return group_.begin() + have_; }

private:
  std::vector<Kmer<T>> group_;
  size_t want_;
  size_t have_;
};

/*
template <uint16_t T> class KmerGroup {
public:
  KmerGroup(moodycamel::ProducerToken&& pt, moodycamel::ConsumerToken&& ct)
      : pt_(std::move(pt)), ct_(std::move(ct)) {}
  moodycamel::ConsumerToken& consumerToken() { return ct_; }
  moodycamel::ProducerToken& producerToken() { return pt_; }
  // get a reference to the chunk this KmerGroup owns
  std::unique_ptr<KmerChunk<T>>& chunkPtr() { return chunk_; }
  // get a *moveable* reference to the chunk this KmerGroup owns
  std::unique_ptr<KmerChunk<T>>&& takeChunkPtr() { return std::move(chunk_); }
  inline void have(size_t num) { chunk_->have(num); }
  inline size_t size() { return chunk_->size(); }
  inline size_t want() const { return chunk_->want(); }
  T& operator[](size_t i) { return (*chunk_)[i]; }
  typename std::vector<T>::iterator begin() { return chunk_->begin(); }
  typename std::vector<T>::iterator end() {
    return chunk_->begin() + chunk_->size();
  }
  void setChunkEmpty() { chunk_.release(); }
  bool empty() const { return chunk_.get() == nullptr; }

private:
  std::unique_ptr<KmerChunk<T>> chunk_{nullptr};
  moodycamel::ProducerToken pt_;
  moodycamel::ConsumerToken ct_;
};
*/

// Iterator class to iterate over KMC databases on disk.
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

    typedef const Kmer<k>* const_ptr_t;



private:

    const Kmer_Container<k>* kmer_container; // The associated k-mer container on which to iterate on.
    CKMCFile kmer_database_input;   // The input reader object (from KMC databases).
    CKmerAPI kmer_object;   // Current KMC k-mer object that this iterator is holding.
    Kmer<k> kmer;   // K-mer present inside the `kmer_object` api.
    bool at_begin;  // Whether this iterator points to the beginning of the KMC database or not.
    bool at_end;  // Whether this iterator points to the beginning of the KMC database or not.
    bool started;
    std::atomic_uint_least64_t num_pushed{0};
    std::atomic_uint_least64_t num_popped{0};

    std::unique_ptr<std::thread> parsing_thread{nullptr};
    moodycamel::ReaderWriterQueue<KmerChunk<k>*> rwq;//(100000);       // Reserve space for at least 100,000 elements up front
    std::atomic_bool finished_parsing{false};
    int thread_result;
    bool was_advanced{false};
    uint64_t num_advances{0};
    KmerChunk<k>* cur_chunk{nullptr};
    typename std::vector<Kmer<k>>::iterator cur_chunk_it;

    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position based on the value of `at_begin`.
    Kmer_Buffered_Iterator(const Kmer_Container<k>* kmer_container, bool at_begin = true, bool at_end = false);

    // Opens the KMC database (internally buffered) to read k-mers.
    void open_kmer_database();

    // Advances the iterator forward by offset one.
    void advance();

    void start();

public:

    // Copy constructs an iterator from the another one `other`.
    Kmer_Buffered_Iterator(const iterator& other);

    ~Kmer_Buffered_Iterator();

    // Assigns the iterator `rhs` to this one, and returns the new iterator.
    const iterator& operator=(const iterator& rhs);

    // Returns a `Kmer` object corresponding to the iterator.
    value_type operator*() ;

    // Returns a const pointer to the `Kmer` object corresponding to the iterator.
    const_ptr_t operator->() ;

    // Advances the iterator by offset one, and returns the new iterator.
    const iterator& operator++();

    /*
    // Advances the iterator by offset one, and returns the old iterator.
    iterator operator++(int);
    */

    // Returns true iff this and `rhs` -- both the iterators refer to the same
    // container and have the same KMC k-mer object.
    bool operator==(const iterator& rhs) const;

    // Returns true iff the iterators this and `rhs` -- either they refer to
    // different containers, or have different KMC k-mer objects.
    bool operator!=(const iterator& rhs) const;
};


template <uint16_t k>
int parse_kmers(const Kmer_Container<k>* const kmer_container,
                CKMCFile& kmer_database_input,
                moodycamel::ReaderWriterQueue<KmerChunk<k>*>& rwq,
                std::atomic_bool& finished_parsing,
                std::atomic_uint_least64_t& num_pushed) {
    
    auto kmer_object = CKmerAPI(kmer_container->kmer_length());
    Kmer<k> kmer;
    bool more_to_read = true;
    static constexpr const size_t chunk_size = 10000;
    KmerChunk<k> *kc = new KmerChunk<k>(chunk_size);
    size_t cursize{0};
    while( (more_to_read = kmer_database_input.ReadNextKmer(kmer_object)) and !finished_parsing) {
        kmer.from_CKmerAPI(kmer_object);
        if (cursize < chunk_size) {
            (*kc)[cursize] = kmer;
            ++cursize;
        } else {
            kc->have(cursize);
            // dump what we currently have onto the queue
            while (!rwq.try_enqueue(kc) and !finished_parsing) {
                // busy wait
            }
            num_pushed += cursize;
            kc = new KmerChunk<k>(chunk_size);
            cursize = 0;
            (*kc)[cursize] = kmer;
            cursize++;
        }   
    }

    // if there is a remaining chunk, push it
    if (cursize > 0) {
        kc->have(cursize);
        while (!rwq.try_enqueue(kc) and !finished_parsing) {
                // busy wait
       }
       num_pushed += cursize;
    }

    std::cerr << "\n\n TOTAL PUSHED = " << num_pushed << "\n\n"; 
    //kmer_object = CKmerAPI();
    kmer_database_input.Close();
    finished_parsing = true;
    return 0;
}

template <uint16_t k>
inline Kmer_Buffered_Iterator<k>::~Kmer_Buffered_Iterator() {
    if (was_advanced) { std::cerr << "think about what to do here\n"; }
    finished_parsing = true;
    if (parsing_thread) { parsing_thread->join(); }
    if (num_popped + num_pushed > 0) { 
        std::cerr << "\n\n {DESTRUCTOR\n";
        std::cerr << "\t\tTOTAL PUSHED = " << num_pushed << "\n"; 
        std::cerr << "\t\tTOTAL POPPED = " << num_popped << "\n"; 
        std::cerr << " }\n\n";
    }
}

template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::start() {
    if (!started) {
    std::cerr << "\n\n actually starting to increment a unique iterator \n\n";
        rwq = moodycamel::ReaderWriterQueue<KmerChunk<k>*>(100);
        open_kmer_database();
        // start background thread
        parsing_thread.reset(new std::thread([this]() {
        this->thread_result = parse_kmers<k>(this->kmer_container, 
                   this->kmer_database_input,
                   this->rwq, this->finished_parsing, this->num_pushed);
        }));
        advance();
        started=true;
    }
}
 
template <uint16_t k>
inline Kmer_Buffered_Iterator<k>::Kmer_Buffered_Iterator(const Kmer_Container<k>* const kmer_container, const bool at_begin, const bool at_end):
    kmer_container(kmer_container), kmer_object(), at_begin(at_begin), at_end(at_end), started(false)
{
    if(at_begin)
    {

    } else if(at_end) {
        num_advances = std::numeric_limits<uint64_t>::max();
        finished_parsing = true;
    }
}


template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::open_kmer_database()
{
    if(!kmer_database_input.OpenForListing(kmer_container->container_location()))
    {
        std::cerr << "Error opening KMC database with prefix " << kmer_container->container_location() << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_Buffered_Iterator<k>::advance()
{
    if (cur_chunk == nullptr || (++cur_chunk_it == cur_chunk->end())) {

        // if we got here because we exhausted the current chunk,
        // then be sure to free the memory for it.
        if (cur_chunk != nullptr) { delete cur_chunk; cur_chunk = nullptr; }
        bool got_work{false};
        // current chunk is exhausted, get the next chunk
        while(!(got_work = rwq.try_dequeue(cur_chunk)) and !finished_parsing) {}
        if (got_work) { 
            // we "failed" the while because we got a chunk
            // set the current chunk iterator to the begin iterator of the chunk
            cur_chunk_it = cur_chunk->begin();
            // set the current k-mer to be the deref of that iterator
            kmer = *cur_chunk_it;
            num_popped += 1;
        } else {
            // if we failed the above while because parsing is done, then there is nothing to do
            num_advances = std::numeric_limits<uint64_t>::max(); at_end = true; kmer_object = CKmerAPI(); 
        }
    } else {
        // deref the next kmer
        // NOTE: we don't have to increment the iterator here b/c it 
        // happens in the `if` statement
        num_popped += 1;
        kmer = *cur_chunk_it;
    }
}


template <uint16_t k>
inline Kmer_Buffered_Iterator<k>::Kmer_Buffered_Iterator(const iterator& other):
    kmer_container(other.kmer_container), /*kmer_object(other.kmer_object),*/ kmer(other.kmer), num_advances(other.num_advances), 
    started(other.started), at_begin(other.at_begin), at_end(other.at_end)
{
    num_pushed.store(other.num_pushed.load());
    num_popped.store(other.num_popped.load());
    if(at_begin)
    {
        std::cerr << "copy constructed iterator after " << other.num_advances << " advances! with AT_BEGIN TRUE.\n";
    } else {
        std::cerr << "copy constructed iterator after " << other.num_advances << " advances!";
        if(other.at_end) { std::cerr << " at END!\n"; } else { std::cerr << " NOT at END!\n"; }
    }
}


template <uint16_t k>
inline const Kmer_Buffered_Iterator<k>& Kmer_Buffered_Iterator<k>::operator=(const iterator& rhs)
{
    kmer_container = rhs.kmer_container;
    //kmer_object = rhs.kmer_object;
    num_advances = rhs.num_advances;
    kmer = rhs.kmer;
    at_begin = rhs.at_begin;
    at_end = rhs.at_end;
    started = rhs.started;
    num_pushed = rhs.num_pushed;
    num_popped = rhs.num_popped;
    if(at_begin)
    {
        std::cerr << "copy assigned iterator after " << rhs.num_advances << " advances! with AT_BEGIN TRUE.\n";
    } else {
        //finished_parsing = true;
        std::cerr << "copy assigned iterator after " << rhs.num_advances << " advances!";
        if(rhs.at_end) { std::cerr << " at END!\n"; } else { std::cerr << " NOT at END!\n"; }
    }

    return *this;
}


template <uint16_t k>
inline typename Kmer_Buffered_Iterator<k>::value_type Kmer_Buffered_Iterator<k>::operator*() 
{
    if (!started) { start(); }
    return kmer;
}


template <uint16_t k>
inline typename Kmer_Buffered_Iterator<k>::const_ptr_t Kmer_Buffered_Iterator<k>::operator->() 
{
    if (!started) { start(); }
    return &kmer;
}


template <uint16_t k>
inline const Kmer_Buffered_Iterator<k>& Kmer_Buffered_Iterator<k>::operator++()
{
    if (!at_end) {
        ++num_advances;
        advance();
        //if (num_advances % 10000 == 1) {
        //    std::cerr << "approx queue size = " << rwq.size_approx() << " \n\n";
        //}
        at_begin = false;
    }
    return *this;
}

/*
template <uint16_t k>
inline Kmer_Buffered_Iterator<k> Kmer_Buffered_Iterator<k>::operator++(int)
{
    Kmer_Buffered_Iterator curr(*this);

    advance();
    if(at_begin)
        at_begin = false;

    return curr;
}
*/


template <uint16_t k>
inline bool Kmer_Buffered_Iterator<k>::operator==(const iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && num_advances == rhs.num_advances;
    // kmer == rhs.kmer;
}


template <uint16_t k>
inline bool Kmer_Buffered_Iterator<k>::operator!=(const iterator& rhs) const
{
    return !(this->operator==(rhs));
}



#endif
