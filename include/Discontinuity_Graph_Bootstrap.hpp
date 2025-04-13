
#ifndef DISCONTINUITY_GRAPH_BOOTSTRAP_HPP
#define DISCONTINUITY_GRAPH_BOOTSTRAP_HPP



#include "Kmer.hpp"
#include "Edge_Matrix.hpp"

#include <cstdint>
#include <cstddef>
#include <string>



namespace cuttlefish
{

// =============================================================================
// Bootstrapper of the corresponding discontinuity-graph of a de Bruijn graph
// (of `k`-mers) from its compacted variant.
template <uint16_t k>
class Discontinuity_Graph_Bootstrap
{
private:

    const std::string cdbg_path;    // Path to the compacted de Bruijn graph FASTA file.
    const uint16_t l;   // Minimizer size.
    const uint64_t minimizer_seed;  // Seed used in computing minimizer-hashes.

    Edge_Matrix<k>& E;  // The edge-matrix of the discontinuity-graph.

    const std::string unitigs_path; // Path-prefix to the output files containing unitigs.
    const std::size_t unitig_buckets;   // Number of buckets for unitigs.


public:

    // Constructs a discontinuity graph bootstrapper from the compacted de
    // Bruijn graph at the path `cdbg_path`. `l`-minimizers are used for the
    // graph, and the seed-value `seed` is used in minimizer-hashing. The graph
    // is constructed at the edge-matrix `E`. The generated unitigs are stored
    // at path-prefix `unitigs_path` into `unitig_buckets` number of files.
    Discontinuity_Graph_Bootstrap(const std::string& cdbg_path, uint16_t l, Edge_Matrix<k>& E, const std::string& unitigs_path, std::size_t unitig_buckets, uint64_t seed = 0);

    // Generates the discontinuity-graph
    void generate();
};

}



#endif
