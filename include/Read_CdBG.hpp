
#ifndef READ_CDBG_HPP
#define READ_CDBG_HPP



#include "globals.hpp"
#include "Build_Params.hpp"
#include "Data_Logistics.hpp"
#include "Kmer_Hash_Table.hpp"
#include "dBG_Info.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>


template <uint16_t k> class kmer_Enumeration_Stats;
template <uint16_t k> class Kmer_Index;


// Read de Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class Read_CdBG
{
private:

    const Build_Params params;  // Required parameters (wrapped inside).
    const Data_Logistics logistics; // Data logistics manager for the algorithm execution.
    std::unique_ptr<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>> hash_table; // Hash table for the vertices (canonical k-mers) of the graph.

    dBG_Info<k> dbg_info;   // Wrapper object for structural information of the graph.

    Kmer_Index<k>* const kmer_idx;  // Index over the k-mers of the graph.

    static constexpr double bits_per_vertex = 9.71; // Expected number of bits required per vertex by Cuttlefish 2.


    // Enumerates the edges of the de Bruijn graph and returns summary statistics of the
    // enumearation.
    kmer_Enumeration_Stats<k + 1> enumerate_edges() const;

    // Enumerates the vertices of the de Bruijn graph using at most `max_memory` amount of
    // memory, and returns summary statistics of the enumeration.
    kmer_Enumeration_Stats<k> enumerate_vertices(std::size_t max_memory) const;

    // Constructs the Cuttlefish hash table for the `vertex_count` vertices of the graph.
    // If `load` is specified, then it is loaded from disk.
    void construct_hash_table(uint64_t vertex_count, bool load = false);

    // Computes the states of the automata, i.e. the vertices of the graph.
    void compute_DFA_states();

    // Extracts the maximal unitigs from the graph.
    void extract_maximal_unitigs();

    // Returns `true` iff the compacted de Bruijn graph to be built from the parameters
    // collection `params` had been constructed in an earlier execution.
    // NB: only the existence of the output meta-info file is checked for this purpose.
    bool is_constructed() const;

    // Returns the maximum temporary disk-usage incurred by some execution of the algorithm,
    // that has its edges-enumeration stats in `edge_stats` and vertices-enumeration stats
    // in `vertex_stats`.
    static std::size_t max_disk_usage(const kmer_Enumeration_Stats<k + 1>& edge_stats, const kmer_Enumeration_Stats<k>& vertex_stats);


public:

    // Constructs a `Read_CdBG` object with the parameters required for the construction of
    // the compacted representation of the underlying read de Bruijn graph wrapped in `params`.
    Read_CdBG(const Build_Params& params);

    // Constructs a `Read_CdBG` object that is to be used to index the k-mers of the underlying
    // de Bruijn graph by `kmer_idx`, with the parameters required for the construction of the
    // compacted representation of the de Bruijn graph wrapped in `params`.
    Read_CdBG(const Build_Params& params, Kmer_Index<k>* kmer_idx);

    // Destructs the compacted graph builder object, freeing its hash table and dumping the
    // graph information to disk.
    ~Read_CdBG();

    // Constructs the compacted read de Bruijn graph, employing the parameters received
    // with the object-constructor.
    void construct();


// The following stuffs are not used anymore with the current algorithm.
/*
private:

    // Extracts the detached chordless cycles of the graph and appends the output to the
    // output file at path `output_file_path`. Specifying `rerun` implies that the graph
    // has been compacted earlier in a separate run of Cuttlefish; otherwise it's done
    // in this same run. Returns `true` iff either there is no DCC in the graph, or the
    // DCCs have already been extracted earlier.
    bool extract_DCCs(const std::string& output_file_path, bool rerun = false);

    // Returns `true` iff the graph contains detached chordless cycles and the current
    // execution is configured to extract those in this same run.
    bool extract_DCCs_this_run() const;

    // Returns `true` iff the data structures required for DCC-extraction is present
    // from an earlier execution of the algorithm.
    bool DCC_data_structs_exist() const;
*/
};



#endif
