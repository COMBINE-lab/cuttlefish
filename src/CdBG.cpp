
#include "CdBG.hpp"
#include "utility.hpp"


template <uint16_t k> 
CdBG<k>::CdBG(const Build_Params& params):
    params(params)
{
    Kmer<k>::set_k(params.k());
}


template <uint16_t k> 
void CdBG<k>::construct()
{
    std::cout << "\nConstructing the minimal perfect hash function (MPHF).\n";
    Vertices.construct(params.vertex_db_path(), params.thread_count(), params.working_dir_path(), params.mph_file_path());

    if(params.remove_kmc_db())
    {
        remove_kmer_set(params.vertex_db_path());
        std::cout << "\nRemoved the KMC database from disk.\n";
    }

    if(!params.output_file_path().empty())
    {
        std::cout << "\nComputing the vertex-states.\n";
        classify_vertices();

        std::cout << "\nOutputting the compacted de Bruijn graph.\n";
        output_maximal_unitigs();
    }

    Vertices.clear();
}


template <uint16_t k> 
size_t CdBG<k>::search_valid_kmer(const char* const seq, const size_t left_end, const size_t right_end) const
{
    size_t valid_start_idx;
    uint16_t base_count;
    

    size_t idx = left_end;
    while(idx <= right_end)
    {
        // Go over the contiguous subsequence of 'N's.
        for(; idx <= right_end && Kmer<k>::is_placeholder(seq[idx]); idx++);

        // Go over the contiguous subsequence of non-'N's.
        if(idx <= right_end)
        {
            valid_start_idx = idx;
            base_count = 0;

            for(; idx <= right_end + k - 1 && !Kmer<k>::is_placeholder(seq[idx]); ++idx)
                if(++base_count == k)
                    return valid_start_idx;
        }
    }


    return right_end + 1;
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
