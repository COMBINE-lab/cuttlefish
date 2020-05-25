
#ifndef CDBG_BUILDER_HPP
#define CDBG_BUILDER_HPP


#include "globals.hpp"
#include "Vertex.hpp"
#include "Vertex_Encoding.hpp"
#include "Kmer_Str.hpp"
#include "Kmer.hpp"
#include "Kmer_HashTable.hpp"

#include <string>
#include <map>


class CdBG_Builder
{
private:
    std::string ref_file;   // Name of the file containing all the newline-separated references.
    uint16_t k; // The k parameter for the edge-centric de Bruijn graph to be compacted.
    Kmer_HashTable Vertices;


    // Classifies the vertices into different types.
    void classify_vertices();

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the first k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, `next_kmer_hat` is the canoninal version of the next k-mer
    // in the sequence, and `next_nucl` is the nucletiode succeeding the first k-mer.
    void process_first_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t next_nucl);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the last k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, and `next_nucl` is the nucletiode preceding the last k-mer.
    void process_last_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::nucleotide_t prev_nucl);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // some k-mer of some sequence, where the k-mer is encountered in the direction
    // `dir`, `next_kmer_hat` is the canoninal version of the next k-mer in the
    // sequence, `prev_nucl` is the nucletiode preceding the k-mer, and `next_nucl`
    // is the nucletiode succeeding the k-mer.
    void process_internal_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t& dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t prev_nucl, const cuttlefish::nucleotide_t next_nucl);

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` forms a
    // self loop with the canonical k-mer `next_kmer_hat` in the sequence. This
    // should only be used with the canonical versions of two adjacent k-mers.
    bool is_self_loop(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_t& next_kmer_hat) const;

    // Outputs all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) to a file named `output_file`.
    void output_maximal_unitigs(const std::string& output_file);

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` traversed
    // in the direction `dir` starts a maximal unitig, where `prev_kmer_hat` and
    // `prev_kmer_dir` are the canonical version and the direction of the previous
    // k-mer in the sequence, respectively.
    bool is_unipath_start(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& prev_kmer_hat, const cuttlefish::kmer_dir_t prev_kmer_dir) const;

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` traversed
    // in the direction `dir` starts a maximal unitig, where `next_kmer_hat` and
    // `next_kmer_dir` are the canonical version and the direction of the next
    // k-mer in the sequence, respectively.
    bool is_unipath_end(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::kmer_dir_t next_kmer_dir) const;

    // Outputs the unitig at the index range [`start_idx`, `end_idx`] of the sequence
    // `seq` (if the unitig had not been output already), to the stream `output`.
    void output_unitig(const char* ref, const uint32_t start_idx, const uint32_t end_idx, std::ofstream& output);



public:
    CdBG_Builder()
    {}

    CdBG_Builder(const std::string& ref_file, const uint16_t k);

    // Construct the compacted de Bruijn graph, and output the maximal unitigs
    // into the file named `output_file`.
    void construct(const std::string& output_file);

    // For debugging.
    void print_vertices() const;
};


#endif
