
#ifndef CdBG_HPP
#define CdBG_HPP


#include <cstdint>
#include <string>
#include <map>
#include <iostream>

#include "Kmer.hpp"
#include "Vertex.hpp"


class CdBG
{
private:
    std::string ref_file;
    uint16_t k;
    std::map<Kmer, Vertex> Vertices;

    bool is_self_loop(const std::string &ref, const Kmer& kmer, const uint32_t kmer_idx);

    void process_first_kmer(const std::string& ref);

    void process_last_kmer(const std::string& ref);

    void process_internal_kmer(const std::string& ref, const uint32_t kmer_idx);

    bool is_unipath_start(const std::string& ref, const uint32_t kmer_idx) const;

    bool is_unipath_end(const std::string &ref, const uint32_t kmer_idx) const;

    void output_unipath(const std::string& ref, std::ofstream &output, const uint32_t start_idx, const uint32_t end_idx);

public:
    CdBG()
    {}

    CdBG(std::string& ref_file, uint8_t k):
        ref_file(ref_file), k(k)
    {}

    void construct();

    void output_maximal_unitigs(const std::string& output_file);

    void print_vertices() const;
};


#endif
