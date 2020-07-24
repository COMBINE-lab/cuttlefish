
#include "CdBG.hpp"


// Initialize the static fields required for the GFA output.
const std::string CdBG::GFA1_HEADER = "H\tVN:Z:1.0";
const std::string CdBG::GFA2_HEADER = "H\tVN:Z:2.0";
const std::string CdBG::PATH_OUTPUT_PREFIX = "tmp-path-output-";
const std::string CdBG::OVERLAP_OUTPUT_PREFIX = "tmp-overlap-output-";


CdBG::CdBG(const std::string& ref_file, const uint16_t k):
    ref_file(ref_file), k(k)
{
    Kmer::set_k(k);
}


void CdBG::construct(const std::string& kmc_file_name, const std::string& bbhash_file_name, const uint16_t thread_count, const std::string& output_file_name)
{
    std::cout << "Constructing the minimal perfect hash function.\n";
    Vertices.construct(kmc_file_name, bbhash_file_name, thread_count);

    std::cout << "Classifying the vertices.\n";
    classify_vertices(thread_count);

    std::cout << "Outputting the maximal unitigs.\n";
    // output_maximal_unitigs(output_file_name, thread_count);
    output_maximal_unitigs_gfa(output_file_name, 2, thread_count);

    Vertices.clear();
}


size_t CdBG::search_valid_kmer(const char* const seq, const size_t left_end, const size_t right_end)
{
    size_t valid_start_idx;
    uint16_t nucl_count;
    

    size_t idx = left_end;
    while(idx <= right_end)
    {
        // Go over the contiguous subsequence of 'N's.
        for(; idx <= right_end && seq[idx] == cuttlefish::PLACEHOLDER_NUCLEOTIDE; idx++);

        // Go over the contiguous subsequence of non-'N's.
        if(idx <= right_end)
        {
            valid_start_idx = idx;
            nucl_count = 0;

            for(; idx <= right_end + k - 1 && seq[idx] != cuttlefish::PLACEHOLDER_NUCLEOTIDE; ++idx)
                if(++nucl_count == k)
                    return valid_start_idx;
        }
    }


    return right_end + 1;
}
