
#include "CdBG.hpp"
#include "kmer_Enumerator.hpp"
#include "Kmer_Container.hpp"
#include "kmer_Enumeration_Stats.hpp"


template <uint16_t k> 
CdBG<k>::CdBG(const Build_Params& params):
    params(params),
    logistics(this->params),
    hash_table(nullptr),
    dbg_info(params.json_file_path())
{
    Kmer<k>::set_k(params.k());
}


template <uint16_t k>
CdBG<k>::~CdBG()
{
    if(hash_table != nullptr)
        hash_table->clear();

    dbg_info.dump_info();
}


template <uint16_t k> 
void CdBG<k>::construct()
{
    if(is_constructed())
    {
        std::cout << "\nThe compacted de Bruijn graph has been constructed earlier. Check " << dbg_info.file_path() << " for results.\n";
        return;
    }

    dbg_info.add_build_params(params);

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    std::cout << "\nEnumerating the vertices of the de Bruijn graph.\n";
    kmer_Enumeration_Stats<k> vertex_stats = enumerate_vertices();
    vertex_stats.log_stats();

    std::chrono::high_resolution_clock::time_point t_vertex = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the vertex set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_vertex - t_start).count() << " seconds.\n";


    const uint64_t vertex_count = vertex_stats.counted_kmer_count();
    std::cout << "Number of vertices: " << vertex_count << ".\n";


    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    construct_hash_table(vertex_count);

#ifdef CF_DEVELOP_MODE
    if(params.vertex_db_path().empty())
#endif
    if(!params.save_vertices())
        Kmer_Container<k>::remove(logistics.vertex_db_path());

    std::chrono::high_resolution_clock::time_point t_mphf = std::chrono::high_resolution_clock::now();
    std::cout << "Constructed the minimal perfect hash function for the vertices. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_mphf - t_vertex).count() << " seconds.\n";
    

    std::cout << "\nComputing the DFA states.\n";
    classify_vertices(short_refs);
    dbg_info.add_short_refs_info(short_refs);

    std::chrono::high_resolution_clock::time_point t_dfa = std::chrono::high_resolution_clock::now();
    std::cout << "Computed the states of the automata. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_dfa - t_mphf).count() << " seconds.\n";


    std::cout << "\nExtracting the maximal unitigs.\n";
    output_maximal_unitigs();

    std::chrono::high_resolution_clock::time_point t_extract = std::chrono::high_resolution_clock::now();
    std::cout << "Extracted the maximal unitigs. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_extract - t_dfa).count() << " seconds.\n";


    const double max_disk = static_cast<double>(max_disk_usage(vertex_stats)) / (1024.0 * 1024.0 * 1024.0);
    std::cout << "\nMaximum temporary disk-usage: " << max_disk << "GB.\n";
}


template <uint16_t k>
kmer_Enumeration_Stats<k> CdBG<k>::enumerate_vertices() const
{
    const KMC::InputFileType ip_type = KMC::InputFileType::MULTILINE_FASTA;
    return kmer_Enumerator<k>().enumerate(
        ip_type, logistics.input_paths_collection(), 1, params.thread_count(),
        params.max_memory(), params.strict_memory(), params.strict_memory(), bits_per_vertex,
        logistics.working_dir_path(), logistics.vertex_db_path()
    );
}


template <uint16_t k>
void CdBG<k>::construct_hash_table(const uint64_t vertex_count)
{
    std::size_t max_memory = std::max(process_peak_memory(), params.max_memory() * 1024U * 1024U * 1024U);
    max_memory = (max_memory > parser_memory ? max_memory - parser_memory : 0);

    hash_table = (params.strict_memory() ?
                    std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_REF_KMER>>(logistics.vertex_db_path(), vertex_count, max_memory) :
                    std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_REF_KMER>>(logistics.vertex_db_path(), vertex_count, max_memory, std::numeric_limits<double>::max()));

    hash_table->construct(params.thread_count(), logistics.working_dir_path(), params.mph_file_path(), params.save_mph());
}


template <uint16_t k>
bool CdBG<k>::is_constructed() const
{
    return file_exists(params.json_file_path());
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


template <uint16_t k>
std::size_t CdBG<k>::max_disk_usage(const kmer_Enumeration_Stats<k>& vertex_stats)
{
    return std::max(vertex_stats.temp_disk_usage(), vertex_stats.db_size());
}


template <uint16_t k>
const Unipaths_Meta_info<k>& CdBG<k>::unipaths_meta_info() const
{
    return unipaths_meta_info_;
}


template <uint16_t k>
uint64_t CdBG<k>::vertex_count() const
{
    return hash_table->size();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
