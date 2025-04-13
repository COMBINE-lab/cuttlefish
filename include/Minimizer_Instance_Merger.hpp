#ifndef MINIMIZER_INSTANCE_MERGER_HPP
#define MINIMIZER_INSTANCE_MERGER_HPP



#include "Minimizer_Instance.hpp"
#include "Minimizer_Instance_Iterator.hpp"
#include "Min_Heap.hpp"

#include <cstdint>
#include <vector>
#include <utility>


// A class to multiway-merge a number of sorted minimizer-instance containers
// of type `T_container_` and collect the union instances in sorted order.
template <typename T_container_>
class Minimizer_Instance_Merger
{
private:

    const uint32_t source_count;    // Number of source containers.
    std::vector<T_container_>& source;  // Collection of the input containers.
    std::vector<Minimizer_Instance_Iterator<T_container_>> iterator; // Collection of iterators over the sources.

    typedef std::pair<Minimizer_Instance, uint16_t> Min_Source_Pair;    // Type of elements to be sorted (merged) through the min heap.
    Min_Heap<Min_Source_Pair> min_heap; // Min heap of minimizer instances and their source container IDs.

public:

    // Constructs a multiway merger for the minimizer instance containers
    // collection `source`.
    Minimizer_Instance_Merger(std::vector<T_container_>& source);

    // Peeks into the next "minimum" minimizer instance from the multiway
    // merger into `min`. Returns `true` iff an instance could be extracted,
    // i.e. the min heap were not empty.
    bool peek(Minimizer_Instance& min) const;

    // Fetches the next "minimum" minimizer instance from the multiway merger
    // into `min`. Returns `true` iff an instance could be extracted,
    // i.e. the min heap were not empty.
    bool next(Minimizer_Instance& min);
};


template <typename T_container_>
inline Minimizer_Instance_Merger<T_container_>::Minimizer_Instance_Merger(std::vector<T_container_>& source):
    source_count(source.size()),
    source(source)
{
    for(uint32_t i = 0; i < source_count; ++i)
    {
        iterator.emplace_back(source[i]);
        if(iterator[i])
        {
            // min_heap.emplace(*iterator[i], i);
            min_heap.push(Min_Source_Pair(*iterator[i], i));
            ++iterator[i];
        }
    }
}


template <typename T_container_>
inline bool Minimizer_Instance_Merger<T_container_>::next(Minimizer_Instance& min)
{
    if(!peek(min))
        return false;

    const auto source_id = min_heap.top().second;
    min_heap.pop();

    if(iterator[source_id])
    {
        // min_heap.emplace(*iterator[source_id], source_id);
        min_heap.push(Min_Source_Pair(*iterator[source_id], source_id));
        ++iterator[source_id];
    }

    return true;
}


template <typename T_container_>
inline bool Minimizer_Instance_Merger<T_container_>::peek(Minimizer_Instance& min) const
{
    if(min_heap.empty())
        return false;

    min = min_heap.top().first;
    return true;
}



#endif