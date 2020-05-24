
#ifndef CDBG_BUILDER_HPP
#define CDBG_BUILDER_HPP


#include "globals.hpp"
#include "Vertex.hpp"
#include "Vertex_Encoding.hpp"
#include "Kmer_Str.hpp"
#include "Kmer.hpp"

#include <string>
#include <map>


class CdBG_Builder
{
private:
    std::string ref_file;   // Name of the file containing all the newline-separated references.
    uint16_t k; // The k parameter for the edge-centric de Bruijn graph to be compacted.
    // std::map<cuttlefish::kmer_t, Vertex> Vertices;  // The set of vertices of the dBG.
    std::map<cuttlefish::kmer_t, Vertex_Encoding> Vertices; // The set of vertices of the dBG.


    // Classifies the vertices into different types.
    void classify_vertices();

    // Processes classification (partially) for the first k-mer of the sequence `seq`.
    void process_first_kmer(const char* seq);

    // Processes classification (partially) for the last k-mer of the sequence `seq`.
    void process_last_kmer(const char* seq, const uint32_t seq_len);

    // Processes classification (partially) for an internal k-mer of the sequence
    // `seq`, at index `kmer_idx`.
    void process_internal_kmer(const char* seq, const uint32_t kmer_idx);

    // Returns a Boolean denoting whether the k-mer `kmer` at index `kmer_idx` of
    // the sequence `seq` forms a self loop with its next k-mer in the sequence
    // (i.e. the k-mer at idx `kmer_idx` + 1). Expects that `kmer` is not the
    // last k-mer of `seq`.
    bool is_self_loop(const char* seq, const cuttlefish::kmer_t& kmer, const uint32_t kmer_idx) const;

    // Outputs all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) to a file named `output_file`.
    void output_maximal_unitigs(const std::string& output_file);

    // Returns a Boolean denoting whether the k-mer at index `kmer_idx` of
    // the sequence `seq` starts a maximal unitig.
    bool is_unipath_start(const char* seq, const uint32_t kmer_idx) const;

    // Returns a Boolean denoting whether the k-mer at index `kmer_idx` of
    // the sequence `seq` ends a maximal unitig.
    bool is_unipath_end(const char* seq, const uint32_t kmer_idx) const;

    // Outputs the unitig at the index range [`start_idx`, `end_idx`] of the sequence
    // `seq` (if the unitig had not been output already), to the stream `output`.
    void output_unitig(const char* ref, const uint32_t start_idx, const uint32_t end_idx, std::ofstream& output);



    // =========================================================================== //

    // Processes classification (partially) for the first k-mer of reference `ref`.
    // void process_first_kmer(const std::string& ref);

    // // Processes classification (partially) for the last k-mer of reference `ref`.
    // void process_last_kmer(const std::string& ref);

    // // Processes classification (partially) for an internal k-mer of reference
    // // `ref`, at index `kmer_idx`.
    // void process_internal_kmer(const std::string& ref, const uint32_t kmer_idx);

    // // Returns a Boolean denoting whether the k-mer `kmer` at index `kmer_idx` of
    // // reference `ref` forms a self loop with its next k-mer in the reference
    // // (i.e. k-mer at idx `kmer_idx` + 1). Expects that `kmer` is not the last
    // // k-mer of 'ref`.
    // bool is_self_loop(const std::string& ref, const cuttlefish::kmer_t& kmer, const uint32_t kmer_idx) const;

    // // Outputs the unitig at the index range [`start_idx`, `end_idx`] of reference
    // // `ref` (if the unitig had not been output already), to the stream `output`.
    // void output_unitig(const std::string& ref, const uint32_t start_idx, const uint32_t end_idx, std::ofstream& output);

    // // Returns a Boolean denoting whether the k-mer at index `kmer_idx` of
    // // reference `ref` starts a maximal unitig.
    // bool is_unipath_start(const std::string& ref, const uint32_t kmer_idx) const;

    // // Returns a Boolean denoting whether the k-mer at index `kmer_idx` of reference
    // // `ref` ends a maximal unitig.
    // bool is_unipath_end(const std::string& ref, const uint32_t kmer_idx) const;


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
