
#include "Build_Params.hpp"
#include "utility.hpp"

#include <iostream>


Build_Params::Build_Params( const bool is_read_graph,
                            const bool is_ref_graph,
                            const std::optional<std::vector<std::string>>& seq_paths,
                            const std::optional<std::vector<std::string>>& list_paths,
                            const std::optional<std::vector<std::string>>& dir_paths,
                            const uint16_t k,
                            const std::optional<uint32_t> cutoff,
                            const bool color,
                            const std::size_t subgraph_count,
                            const std::size_t vertex_part_count,
                            const std::size_t lmtig_bucket_count,
                            const std::size_t gmtig_bucket_count,
                            const std::string& vertex_db_path,
                            const std::string& edge_db_path,
                            const uint16_t thread_count,
                            const std::optional<std::size_t> max_memory,
                            const bool strict_memory,
                            const bool idx,
                            const uint16_t min_len,
                            const std::string& output_file_path,
                            const std::optional<cuttlefish::Output_Format> output_format,
                            const bool track_short_seqs,
                            const bool poly_n_stretch,
                            const std::string& working_dir_path,
                            const bool path_cover,
                            const bool save_mph,
                            const bool save_buckets,
                            const bool save_vertices
#ifdef CF_DEVELOP_MODE
                            , const double gamma
#endif
                    ):
    is_read_graph_(is_read_graph),
    is_ref_graph_(is_ref_graph),
    seq_input_(seq_paths, list_paths, dir_paths),
    k_(k),
    cutoff_(cutoff),
    color_(color),
    subgraph_count_(subgraph_count),
    vertex_part_count_(vertex_part_count),
    lmtig_bucket_count_(lmtig_bucket_count),
    gmtig_bucket_count_(gmtig_bucket_count),
    vertex_db_path_(vertex_db_path),
    edge_db_path_(edge_db_path),
    thread_count_(thread_count),
    max_memory_(max_memory),
    strict_memory_(strict_memory),
    idx_(idx),
    min_len_(min_len),
    output_file_path_(output_file_path),
    output_format_(output_format),
    track_short_seqs_(track_short_seqs),
        poly_n_stretch_(poly_n_stretch),
    working_dir_path_(working_dir_path.back() == '/' ? working_dir_path : working_dir_path + "/"),
    path_cover_(path_cover),
    save_mph_(save_mph),
    save_buckets_(save_buckets),
    save_vertices_(save_vertices)
#ifdef CF_DEVELOP_MODE
    , gamma_(gamma)
#endif
{}


const std::string Build_Params::output_file_ext() const
{
    if(is_read_graph() || is_ref_graph())
        return cuttlefish::file_ext::unipaths_ext;

    switch(output_format())
    {
    case cuttlefish::Output_Format::fa:
        return cuttlefish::file_ext::unipaths_ext;

    case cuttlefish::Output_Format::gfa1:
        return cuttlefish::file_ext::gfa1_ext;

    case cuttlefish::Output_Format::gfa2:
        return cuttlefish::file_ext::gfa2_ext;

    default:
        break;
    }


    return "";
}


bool Build_Params::is_valid() const
{
    bool valid = true;


    // Input data need to be non-empty.
    if(seq_input_.empty())
    {
        std::cout << "No sequence input provided for compacted de Bruijn graph construction.\n";
        valid = false;
    }


    // Even `k` values are not consistent with the theory.
    // Also, `k` needs to be in the range `(1, MAX_K]`.
    if((k_ % 2) == 0 || k_ == 1 || k_ > cuttlefish::MAX_K)
    {
        std::cout << "The k-mer length (k) needs to be odd and within " << cuttlefish::MAX_K << ".\n";
        valid = false;
    }


    // Unsupported thread counts are to be discarded.
    const auto num_threads = std::thread::hardware_concurrency();
    if(num_threads > 0 && thread_count_ > num_threads)
    {
        std::cout << "At most " << num_threads << " concurrent threads are supported by the machine.\n";
        valid = false;
    }


    // l-minimizer length must be at most k-mer length, if indexing is requested.
    if(idx_ && min_len_ >= k_)
    {
        std::cout << "l-minimizer length can be at most k-mer length.\n";
        valid = false;
    }


    // Output directory must exist.
    const std::string op_dir = dirname(output_file_path_);
    if(!dir_exists(op_dir))
    {
        std::cout << "Output directory " << op_dir << " does not exist.\n";
        valid = false;
    }


    // Working directory must exist.
    const std::string work_dir = dirname(working_dir_path_);
    if(!dir_exists(work_dir))
    {
        std::cout << "Working directory " << work_dir << " does not exist.\n";
        valid = false;
    }


    // Memory budget options should not be mixed with.
    if(max_memory_  && !strict_memory_)
        std::cout << "Both a memory bound and the option for unrestricted memory usage specified. Unrestricted memory mode will be used.\n";


    if(is_read_graph_ || is_ref_graph_) // Validate Cuttlefish 2 specific arguments.
    {
        // Read and reference de Bruijn graph parameters can not be mixed with.
        if(is_read_graph_ && is_ref_graph_)
        {
            std::cout << "Both read and reference de Bruijn graph specified. Please select only one for Cuttlefish 2, or none to use Cuttlefish 1.\n";
            valid = false;
        }


        // A cutoff frequency of 0 is theoretically inconsistent.
        if(cutoff() == 0)
        {
            std::cout << "Cutoff frequency specified to be 0, which is theoretically inconsistent. Please use 1 if you wish to retain all the k-mers without filtering.\n";
            valid = false;
        }

        // Cutoff frequency _should be_ 1 for reference de Bruijn graphs.
        if(is_ref_graph_ && cutoff() != 1)
            std::cout << "WARNING: cutoff frequency specified not to be 1 on reference sequences.\n";


        // Cuttlefish 1 specific arguments can not be specified.
        if(output_format_)
        {
            std::cout << "Cuttlefish 1 specific arguments specified while using Cuttlefish 2.\n";
            valid = false;
        }
    }
    else    // Validate Cuttlefish 1 specific arguments.
    {
        // Invalid output formats are to be discarded.
        if(output_format() >= cuttlefish::num_op_formats)
        {
            std::cout << "Invalid output file format.\n";
            valid = false;
        }


        // Cuttlefish 2 specific arguments can not be specified.
        if(cutoff_ || path_cover_)
        {
            std::cout << "Cuttlefish 2 specific arguments specified while using Cuttlefish 1.\n";
            valid = false;
        }
    }


    // Develop-mode options can not to be provided in regular use.
#ifndef CF_DEVELOP_MODE
    if(!vertex_db_path_.empty() || !edge_db_path_.empty())
    {
        std::cout << "Paths to vertex- and edge-databases are supported only in develop mode.\n";
        valid = false;
    }
#endif


    return valid;
}
