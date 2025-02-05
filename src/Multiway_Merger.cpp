
#include "Multiway_Merger.hpp"
#include "globals.hpp"


template <uint16_t k>
Multiway_Merger<k>::Multiway_Merger(const std::vector<std::string>& db_list):
    db_count(db_list.size()),
    db_list(db_list),
    iterator(db_count),
    db_id(db_count),
    kmer_count(db_count, 0)
{
    for(uint32_t i = 0; i < db_count; ++i)
    {
        iterator[i].init(db_list[i]);
        db_id.emplace(db_list[i], i);
    }
}


template <uint16_t k>
Multiway_Merger<k>::~Multiway_Merger()
{}


template <uint16_t k>
void Multiway_Merger<k>::launch()
{
    for(uint32_t i = 0; i < db_count; ++i)
    {
        iterator[i].launch();

        if(iterator[i].parse_kmers_atomic(kmers))
        {
            std::for_each(kmers.begin(), kmers.end(), [&](const Kmer<k>& kmer) { kmer_source_pairs.emplace_back(kmer, i); });

            kmer_count[i] = kmers.size();
        }
    }

    min_heap.init_heap(kmer_source_pairs);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Multiway_Merger)
