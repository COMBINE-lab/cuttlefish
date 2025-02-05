
#ifndef PATH_INFO_HPP
#define PATH_INFO_HPP



#include "Kmer.hpp"
#include "globals.hpp"
#include "xxHash/xxhash.h"

#include <cstdint>


namespace cuttlefish
{

// =============================================================================
// Path-information of an object in a discontinuity graph: its path-ID, rank in
// a fixed traversal of the path, orientation in that traversal, and whether it
// actually forms a cycle (abusing notation). Path-IDs are k-mers.
template <uint16_t k>
class Path_Info
{
private:

    path_id_t<k> p_;    // The path-ID.
    weight_t r_;    // The rank.
    side_t o_;  // The orientation of the object in its specified rank—the path traversal exits it through the side `o`.
    bool is_cycle_; // Whether the path is a cycle (abusing notation).


public:

    Path_Info()
    {}


    // Constructs a path-info object for an object such that its path-ID is `p`
    // and rank in the path is `r` when the path is traversed in the
    // orientation such that the traversal exits the object through its side
    // `o`. `is_cycle` denotes whether the path is a cycle (abusing notation).
    Path_Info(const path_id_t<k> p, const weight_t r, const side_t o, const bool is_cycle):
          p_(p)
        , r_(r)
        , o_(o)
        , is_cycle_(is_cycle)
    {}


    // Returns the path-ID.
    auto p() const { return p_; }

    // Returns the rank.
    auto r() const { return r_; }

    // Returns whether the path is a cycle (abusing notation).
    auto is_cycle() const { return is_cycle_; }

    // Returns the orientation `o` of the object in its specified rank—the
    // path traversal exits the object through the side `o`.
    auto o() const { return o_; }

    // Returns `true` iff this information is the same as in `rhs`.
    bool operator==(const Path_Info& rhs) const { return p_ == rhs.p_ && r_ == rhs.r_ && o_ == rhs.o_; }

    // Returns `true` iff this information is lexicographically smaller than
    // `rhs`.
    bool operator<(const Path_Info& rhs) const { return p_ != rhs.p_ ? (p_ < rhs.p_) : (r_ < rhs.r_); }

    // Returns a 64-bit hash value of the path-information.
    uint64_t hash() const { return XXH3_64bits(&p_, sizeof(p_)) ^ XXH3_64bits(&r_, sizeof(r_)) ^ XXH3_64bits(&o_, sizeof(o_)); }
};


// An object and associated path-information. Path-IDs are k-mers.
template <typename T_, uint16_t k>
class Obj_Path_Info_Pair
{
private:

    T_ obj_;    // The object.
    Path_Info<k> path_info_;    // Path-information of the object.

public:

    Obj_Path_Info_Pair()    // TODO: consider removing.
    {}

    // For an object `obj`, constructs a pairing of it with its path-info
    // specified with its path-ID `p` and rank in the path `r` when the path is
    // traversed in the orientation such that the traversal exits the object
    // through its side `o`. `is_cycle` denotes whether the path is a cycle
    // (abusing notation).
    Obj_Path_Info_Pair(const T_ obj, const path_id_t<k> p, const weight_t r, const side_t o, const bool is_cycle):
          obj_(obj)
        , path_info_(p, r, o, is_cycle)
    {}

    // For an object `obj`, constructs a pairing of it with its path-info
    // specified with `path_info`.
    Obj_Path_Info_Pair(const T_ obj, const Path_Info<k> path_info):
          obj_(obj)
        , path_info_(path_info)
    {}

    // Returns the object.
    auto obj() const { return obj_; }

    // Returns the path-info of the object.
    auto path_info() const { return path_info_; }

    // Returns `true` iff this information is the same as in `rhs`.
    bool operator==(const Obj_Path_Info_Pair& rhs) const { return obj_ == rhs.obj_ && path_info_ == rhs.path_info_; }
};

}



#endif
