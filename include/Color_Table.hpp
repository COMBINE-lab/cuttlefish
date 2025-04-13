
#ifndef COLOR_TABLE_HPP
#define COLOR_TABLE_HPP



// #define USE_PARLAY_HASH


#include "Color_Encoding.hpp"

#ifndef USE_PARLAY_HASH
#include "boost/unordered/concurrent_flat_map.hpp"
#else
#include "parlayhash/parlay_hash/unordered_map.h"
#endif

#include <cstdint>
#include <cassert>


namespace cuttlefish
{

enum class Color_Status;    // Extraction-status of a color-set.


// Hashtable for color-sets. Keys are color-set hashes and values are color-set
// coordinates.
class Color_Table
{

    typedef uint64_t hash_t;
    typedef Color_Coordinate coord_t;

private:

    static constexpr uint64_t map_sz_init = 64 * 1024 * 1024;   // Map has preallocated memory for 64M color-hashes.

#ifndef USE_PARLAY_HASH
    boost::unordered::concurrent_flat_map<hash_t, coord_t> M;
#else
    parlay::parlay_unordered_map<hash_t, coord_t> M;
#endif

public:

    // Constructs an empty color-table.
    Color_Table();

    // Returns the size of the table. NB: work is proportional to the size of
    // the table if `parlayhash` is used.
    auto size() { return M.size(); }

    // Marks that the color with hash `h` is in the process of extraction by
    // the `w`'th worker, if a corresponding entry for `h` does not already
    // exist in the table. If it does, it is put in `c`. Returns the
    // extraction-status of the color prior to this invocation.
    Color_Status mark_in_process(hash_t h, uint64_t w, Color_Coordinate& c);

    // Updates the key `h` with value `c` if `h` is marked as in process of
    // extraction. Returns `true` iff the update is successful.
    bool update_if_in_process(hash_t h, Color_Coordinate c);

    // Returns the value associated to the key `h`.
    Color_Coordinate get(hash_t h);
};


// Extraction-status of a color-set.
enum class Color_Status
{
    undiscovered,   // has not been seen yet
    in_process,     // is in the process of extraction
    discovered,     // completely extracted
};


inline Color_Status Color_Table::mark_in_process(const hash_t h, const uint64_t w, Color_Coordinate& c)
{
    typedef Color_Status status_t;

#ifndef USE_PARLAY_HASH

    const auto r =  M.emplace_or_cvisit(h, Color_Coordinate(w), [&](const auto& p)
                    {
                        c = p.second;
                    });

    return r ? status_t::undiscovered : (c.is_in_process() ? status_t::in_process : status_t::discovered);

#else

    const auto r = M.Insert(h, Color_Coordinate(w));
    c = (r.value_or(c));

    return !r ? status_t::undiscovered : (c.is_in_process() ? status_t::in_process : status_t::discovered);

#endif
}


inline bool Color_Table::update_if_in_process(const hash_t h, const Color_Coordinate c)
{
    bool was_in_process = false;

#ifndef USE_PARLAY_HASH

    const auto r =  M.visit(h, [&](auto& p)
                    {
                        auto& color = p.second;
                        was_in_process = color.is_in_process();
                        color = (was_in_process ? c : color);
                    });
#else

    const auto r =  M.Upsert(h, [&](const auto& cur_c)
                    {
                        assert(cur_c);
                        was_in_process = cur_c->is_in_process();
                        return was_in_process ? c : *cur_c;
                    });

#endif

    assert(r); (void)r;

    return was_in_process;
}


inline Color_Coordinate Color_Table::get(const hash_t h)
{
#ifndef USE_PARLAY_HASH

    assert(M.contains(h));

    Color_Coordinate c;
    const auto r =  M.cvisit(h, [&](const auto& p)
                    {
                        c = p.second;
                    });
    assert(r); (void)r;
    return c;

#else

    const auto c = M.Find(h);
    assert(c);

    return *c;

#endif
}

}



#endif