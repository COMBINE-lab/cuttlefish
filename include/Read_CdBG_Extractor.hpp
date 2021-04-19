
#ifndef READ_CDBG_EXTRACTPR_HPP
#define READ_CDBG_EXTRACTOR_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Build_Params.hpp"
#include "Spin_Lock.hpp"
#include "Thread_Pool.hpp"


// A class to extract the vertices from a compacted de Bruin graph — which are the maximal unitigs of some ordinary de Bruijn graph.
template <uint16_t k>
class Read_CdBG_Extractor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (i.e. canonical k-mers) of the original (uncompacted) de Bruijn graph.

    // Members required to keep track of the total number of vertices processed across different worker (i.e. extractor) threads.
    mutable Spin_Lock lock;
    mutable uint64_t vertices_processed = 0;


    // Distributes the maximal unitigs extraction task — disperses the graph vertices (i.e. k-mers)
    // parsed by the parser `vertex_parser` to the worker threads in the thread pool `thread_pool`,
    // for the unitpath-flanking vertices to be identified and the corresponding unipaths to be extracted.
    void distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* vertex_parser, Thread_Pool<k>& thread_pool);

    // Processes the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`,
    // i.e. for each vertex `v` provided to that thread, identifies whether it is a unipath-flanking
    // vertex, and if it is, then piece-wise constructs the corresponding unipath.
    void process_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);


public:

    // Constructs a vertex-extractor object for some compacted read de Bruijn graph, with the required
    // parameters wrapped inside `params`, and uses the Cuttlefish hash table `hash_table`.
    Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Extracts the maximal unitigs of the de Bruijn graph.
    void extract_maximal_unitigs();
};



#endif
