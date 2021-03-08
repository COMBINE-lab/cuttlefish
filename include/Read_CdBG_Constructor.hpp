
#ifndef READ_CDBG_CONSTRUCTOR_HPP
#define READ_CDBG_CONSTRUCTOR_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Build_Params.hpp"
#include "Thread_Pool.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"


// A class to construct compacted read de Bruijn graphs.
template <uint16_t k>
class Read_CdBG_Constructor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (canonical k-mers) of the graph.
    const Kmer_Container<k + 1> edge_container; // Wrapper container for the edge-database.
    Kmer_SPMC_Iterator<k + 1> edge_parser;  // Parser for the edges from the edge-database.


    // Distributes the DFA-states computation task to the worker threads in the thread pool `thread_pool`.
    void distribute_states_computation(Thread_Pool<k>& thread_pool);

    // Processes the edges provided to the thread with id `thread_id`, i.e. makes state-transitions for
    // the DFA as per the edges provided to that thread.
    void process_edges(uint16_t thread_id);


public:

    // Consructs a read-CdBG builder object, with the required parameters wrapped in `params`, and uses
    // the Cuttlefish hash table `hash_table`.
    Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Computes the states of the DFA in the de Bruijn graph.
    void compute_DFA_states();
};



#endif
