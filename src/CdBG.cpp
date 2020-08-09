
#include "CdBG.hpp"


// Initialize the static fields required for the GFA output.
const std::string CdBG::GFA1_HEADER = "H\tVN:Z:1.0";
const std::string CdBG::GFA2_HEADER = "H\tVN:Z:2.0";
std::string CdBG::PATH_OUTPUT_PREFIX = "cuttlefish-path-output-";
std::string CdBG::OVERLAP_OUTPUT_PREFIX = "cuttlefish-overlap-output-";


CdBG::CdBG(const Build_Params& params):
    params(params), k(params.k())
{
    cuttlefish::kmer_t::set_k(k);
}


void CdBG::construct()
{
    std::cout << "Constructing the minimal perfect hash function.\n";
    Vertices.construct(params.kmc_db_path(), params.thread_count(), params.mph_file_path());

    std::cout << "Classifying the vertices.\n";
    classify_vertices();

    std::cout << "Outputting the maximal unitigs.\n";
    params.output_format() == 0 ? output_maximal_unitigs() : output_maximal_unitigs_gfa();

    Vertices.clear();
}


size_t CdBG::search_valid_kmer(const char* const seq, const size_t left_end, const size_t right_end) const
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
