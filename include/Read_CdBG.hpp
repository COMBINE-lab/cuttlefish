
#ifndef READ_CDBG_HPP
#define READ_CDBG_HPP



#include "globals.hpp"
#include "Build_Params.hpp"
#include "Kmer_Hash_Table.hpp"
#include "dBG_Info.hpp"

#include <memory>


class kmer_Enumeration_Stats;


// Read de Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class Read_CdBG
{
private:

    const Build_Params params;  // Required parameters (wrapped inside).
    std::unique_ptr<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>> hash_table; // Hash table for the vertices (canonical k-mers) of the graph.

    dBG_Info<k> dbg_info;   // Wrapper object for structural information of the graph.


    // Enumerates the edges of the de Bruijn graph in a database at path `edge_db_path`,
    // and returns summary statistics of the enumearation.
    kmer_Enumeration_Stats enumerate_edges(const std::string& edge_db_path);

    // Enumerates the vertices of the de Bruijn graph in a database at path `vertex_db_path`,
    // from the edge database present at `edge_db_path`, using at most `max_memory` amount of
    // memory. Returns summary statistics of the enumeration.
    kmer_Enumeration_Stats enumerate_vertices(const std::string& edge_db_path, const std::string& vertex_db_path, std::size_t max_memory);

    // Constructs the Cuttlefish hash table for the `vertex_count` vertices in the database
    // at path `vertex_db_path`.
    void construct_hash_table(const std::string& vertex_db_path, uint64_t vertex_count);

    // Computes the states of the automata, i.e. the vertices of the graph having it edge
    // set present at the path prefix `edge_db_path`.
    void compute_DFA_states(const std::string& edge_db_path);

    // Extracts the maximal unitigs from the graph gaving its vertex set present at the
    // path prefix `vertex_db_path`.
    void extract_maximal_unitigs(const std::string& vertex_db_path);


public:

    // Constructs a `Read_CdBG` object with the parameters required for the construction of
    // the compacted representation of the underlying read de Bruijn graph wrapped in `params`.
    Read_CdBG(const Build_Params& params);

    // Constructs the compacted read de Bruijn graph, employing the parameters received
    // with the object-constructor.
    void construct();

    // Returns `true` iff the compacted de Bruijn graph to be built from the parameters
    // collection `params` had been constructed in an earlier execution.
    // NB: only the existence of the output meta-info file is checked for this purpose.
    static bool is_constructed(const Build_Params& params);
};



#endif
