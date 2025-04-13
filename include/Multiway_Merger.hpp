
#ifndef MULTIWAY_MERGER_HPP
#define MULTIWAY_MERGER_HPP



#include "Kmer_SPSC_Iterator.hpp"
#include "Kmer.hpp"
#include "Min_Heap.hpp"

#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <algorithm>


// =============================================================================


// A class to multiway-merge a number of k-mer databases and collect the list of
// the union k-mers in sorted order.
template <uint16_t k>
class Multiway_Merger
{
private:

    const uint32_t db_count;    // Number of input k-mer databases.
    const std::vector<std::string> db_list; // List of the paths to the input k-mer databases.
    std::vector<Kmer_SPSC_Iterator<k>> iterator;    // Collection of iterators over the input k-mer databases.

    typedef uint16_t source_id_t;   // Type of the source-identifier with each k-mer pulled from some database.
    std::unordered_map<std::string, source_id_t> db_id;    // Numerical ID of the input k-mer databases.

    class Kmer_Source_Pair; // Type of elements to be sorted (merged) through the min heap.
    Min_Heap<Kmer_Source_Pair> min_heap;    // Min heap of k-mers and their source database IDs.

    std::vector<uint32_t> kmer_count;   // Number of k-mers from each input database currently present in the heap.

    std::vector<Kmer_Source_Pair> kmer_source_pairs;    // Initial collection of elements to start the heap-merge from.
    std::vector<Kmer<k>> kmers; // Container to fetch k-mer chunks from each database separately.


    class Kmer_Source_Pair
    {
        friend class Multiway_Merger;

    private:

        Kmer<k> kmer;
        source_id_t source_id;

    public:
        Kmer_Source_Pair(const Kmer<k>& kmer, const source_id_t source_id):
            kmer(kmer),
            source_id(source_id)
        {}

        bool operator>(const Kmer_Source_Pair& rhs) const
        {
            return kmer != rhs.kmer ? kmer > rhs.kmer : source_id > rhs.source_id;
        }
    };


public:

    // Constructs a multiway merger for the k-mer databases at the paths-
    // collection `db_list`.
    Multiway_Merger(const std::vector<std::string>& db_list);

    // Destructs the multiway merger.
    ~Multiway_Merger();

    // Copy-construction is prohibited. This is a complex object containing
    // k-mer database iterator collections, each in some *arbitrary* state,
    // which should not be allowed to be copied.
    Multiway_Merger(const Multiway_Merger&) = delete;

    // Copy-assignment is prohibited. This is a complex object containing k-mer
    // database iterator collections, each in some *arbitrary* state, which
    // should not be allowed to be copied.
    const Multiway_Merger& operator=(const Multiway_Merger&) = delete;

    // Launches the multiway merge over the provided database list.
    void launch();

    // Fetches the next "minimum" k-mer from the multiway merger into `kmer`,
    // and collects its sources into the color representation `color`. Also
    // collects the source IDs whose k-mer count in the heap becomes empty due
    // to this k-mer extraction, into `empty_sources`. Returns `true` iff a
    // k-mer could be extracted, i.e. the min heap were not empty.
    template <typename T_color_>
    bool next(Kmer<k>& kmer, T_color_& color, std::vector<source_id_t>& empty_sources);

    // Given a list `sources` of databases, tops up the min heap by fetching
    // k-mer chunks from each database into the heap.
    void top_up(const std::vector<source_id_t>& sources);
};


template <uint16_t k>
template <typename T_color_>
inline bool Multiway_Merger<k>::next(Kmer<k>& kmer, T_color_& color, std::vector<source_id_t>& empty_sources)
{
    if(min_heap.empty())
        return false;


    color.clear();
    kmer = min_heap.top().kmer;

    while(!min_heap.empty() && min_heap.top().kmer == kmer)
    {
        const source_id_t source_id = min_heap.top().source_id;
        // color.insert(source_id);
        color.push_back(source_id); // TODO: placeholder until color-rep impls. are designed.

        min_heap.pop();
        if(--kmer_count[source_id] == 0)
            empty_sources.emplace_back(source_id);
    }

    return true;
}


template <uint16_t k>
inline void Multiway_Merger<k>::top_up(const std::vector<source_id_t>& sources)
{
    for(const source_id_t id : sources)
        if(iterator[id].parse_kmers_atomic(kmers))
        {
            kmer_source_pairs.clear();
            std::for_each(kmers.begin(), kmers.end(), [&](const Kmer<k>& kmer) { kmer_source_pairs.emplace_back(kmer, id); });

            kmer_count[id] = kmers.size();
            min_heap.push(kmer_source_pairs);
        }
}



#endif
