
#ifndef CDBG_BUILDER_HPP
#define CDBG_BUILDER_HPP


#include "globals.hpp"
#include "Vertex.hpp"
#include "Vertex_Encoding.hpp"
#include "Kmer.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Annotated_Kmer.hpp"

#include <string>
#include <map>


class CdBG_Builder
{
private:
    std::string ref_file;   // Name of the file containing all the newline-separated references.
    uint16_t k; // The k parameter for the edge-centric de Bruijn graph to be compacted.
    Kmer_Hash_Table Vertices;   // The hash table for the vertices (canonical k-mers) of the de Bruijn graph.


    // Classifies the vertices into different types.
    void classify_vertices();

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the first k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, the canonical version of the next k-mer in the sequence is
    // `next_kmer_hat`, and the nucletiode succeeding the first k-mer is `next_nucl`.
    void process_first_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t next_nucl);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the last k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, and the nucletiode preceding the last k-mer is `prev_nucl`.
    void process_last_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::nucleotide_t prev_nucl);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // some internal k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, the canoninal version of the next k-mer in the sequence is
    // `next_kmer_hat`, the nucletiode preceding the k-mer is `prev_nucl`, and the
    // nucletiode succeeding the k-mer is `next_nucl`.
    void process_internal_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t prev_nucl, const cuttlefish::nucleotide_t next_nucl);

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` forms a
    // self loop with the canonical k-mer `next_kmer_hat` in the sequence. This
    // should only be used with the canonical versions of two adjacent k-mers.
    bool is_self_loop(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_t& next_kmer_hat) const;

    // Outputs all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) to a file named `output_file`.
    void output_maximal_unitigs(const std::string& output_file);

    // Returns a Boolean denoting whether a k-mer with state `state` traversed in
    // the direction `dir` starts a maximal unitig, where `prev_kmer_state` and
    // `prev_kmer_dir` are the state and the direction of the previous k-mer in
    // the sequence, respectively.
    bool is_unipath_start(const cuttlefish::state_t state, const cuttlefish::kmer_dir_t dir, const cuttlefish::state_t prev_kmer_state, const cuttlefish::kmer_dir_t prev_kmer_dir) const;

    // Returns a Boolean denoting whether a k-mer with state `state` traversed in
    // the direction `dir` ends a maximal unitig, where `next_kmer_state` and
    // `next_kmer_dir` are the state and the direction of the next k-mer in the
    // sequence, respectively.
    bool is_unipath_end(const cuttlefish::state_t state, const cuttlefish::kmer_dir_t dir, const cuttlefish::state_t next_kmer_state, const cuttlefish::kmer_dir_t next_kmer_dir) const;

    // Outputs the unitig at the k-mer range between the annotated k-mers
    // `start_kmer` and `end_kmer` of the sequence `seq` (if the unitig had not
    // been output already), to the stream `output`.
    void output_unitig(const char* ref, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, std::ofstream& output);
    
    // Writes the path in the sequence `seq` with its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively,
    // to the stream `output`. If `in_forward` is true, then the string spelled
    // by the path is written; otherwise its reverse complement is written.
    // Note that, the output operation appends a newline at the end.
    void write_path(const char* seq, const uint32_t start_kmer_idx, const uint32_t end_kmer_idx, const bool in_forward, std::ofstream& output) const;


public:
    CdBG_Builder()
    {}

    CdBG_Builder(const std::string& ref_file, const uint16_t k);

    // Constructs the compacted de Bruijn graph with its distinct k-mers collection
    // being present in the file named `kmers_file` (in 64-bit integers) and the
    // count of these k-mers being `kmers_count`, and outputs the maximal unitigs
    // into the file named `output_file`.
    void construct(const std::string& kmers_file, const uint64_t kmer_count, const std::string& output_file);
};


#endif
