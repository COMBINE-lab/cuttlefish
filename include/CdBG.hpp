
#ifndef CdBG_HPP
#define CdBG_HPP


#include <cstdint>
#include <string>
#include <map>
#include <iostream>

#include "globals.hpp"
#include "Kmer.hpp"
#include "Vertex.hpp"


class CdBG
{
private:
    std::string ref_file;               // Name of the file containing all the input references.
    uint16_t k;                         // The `k` value.
    std::map<cuttlefish::kmer_t, Vertex> Vertices;    // The set of vertices of the dBG.

    // Return a Boolean denoting whether the k-mer `kmer` at index `kmer_idx` of
    // reference `ref` forms a self loop with its next k-mer in the reference
    // (i.e. k-mer at idx `kmer_idx` + 1). Expects that `kmer` is not the last
    // k-mer of 'ref`.
    bool is_self_loop(const std::string &ref, const cuttlefish::kmer_t& kmer, const uint32_t kmer_idx);

    // Process classification (partially) for the first k-mer of reference `ref`.
    void process_first_kmer(const std::string& ref);

    // Process classification (partially) for the last k-mer of reference `ref`.
    void process_last_kmer(const std::string& ref);

    // Process classification (partially) for an internal k-mer `kmer` of reference
    // `ref`, at index `kmer_idx`.
    void process_internal_kmer(const std::string& ref, const uint32_t kmer_idx);

    // Return a Boolean denoting whether the k-mer at index `kmer_idx` of reference
    // `ref` starts a maximal unitig.
    bool is_unipath_start(const std::string& ref, const uint32_t kmer_idx) const;

    // Return a Boolean denoting whether the k-mer at index `kmer_idx` of reference
    // `ref` ends a maximal unitig.
    bool is_unipath_end(const std::string &ref, const uint32_t kmer_idx) const;

    // Output the unitig at the index range [`start_idx`, `end_idx`] of reference
    // `ref`(if the unitig had not been output already), to the stream `output`.
    void output_unipath(const std::string& ref, std::ofstream &output, const uint32_t start_idx, const uint32_t end_idx);

public:
    CdBG()
    {}

    CdBG(std::string& ref_file, uint8_t k):
        ref_file(ref_file), k(k)
    {}

    // Construct the compacted de Bruijn graph (performs the vertex classification).
    void construct();

    // Output all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) to a file named `output_file`.
    void output_maximal_unitigs(const std::string& output_file);

    // For debugging purposes: print the vertices information.
    void print_vertices() const;
};


#endif
