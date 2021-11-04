
#include "Build_Params.hpp"
#include "utility.hpp"


Build_Params::Build_Params( const bool is_read_graph,
                            const std::vector<std::string>& seq_paths,
                            const std::vector<std::string>& list_paths,
                            const std::vector<std::string>& dir_paths,
                            const uint16_t k,
                            const uint32_t cutoff,
                            const std::string& vertex_db_path,
                            const std::string& edge_db_path,
                            const uint16_t thread_count,
                            const std::size_t max_memory,
                            const bool strict_memory,
                            const std::string& output_file_path,
                            const uint8_t output_format,
                            const std::string& working_dir_path,
                            const bool remove_kmc_db,
                            const std::string& mph_file_path,
                            const std::string& buckets_file_path,
                            const bool save_vertices,
                            const std::string& json_file_path
#ifdef CF_DEVELOP_MODE
                    , const double gamma
#endif
                    ):
        is_read_graph_(is_read_graph),
        seq_input_(seq_paths, list_paths, dir_paths),
        k_(k),
        cutoff_(cutoff),
        vertex_db_path_(vertex_db_path),
        edge_db_path_(edge_db_path),
        thread_count_(thread_count),
        max_memory_(max_memory),
        strict_memory_(strict_memory),
        output_file_path_(output_file_path),
        output_format_(cuttlefish::Output_Format(output_format)),
        working_dir_path_(working_dir_path.back() == '/' ? working_dir_path : working_dir_path + "/"),
        remove_kmc_db_(remove_kmc_db),
        mph_file_path_(mph_file_path),
        buckets_file_path_(buckets_file_path),
        save_vertices_(save_vertices),
        json_file_path_(json_file_path)
#ifdef CF_DEVELOP_MODE
        , gamma_(gamma)
#endif
    {}


bool Build_Params::is_valid() const
{
    // TODO: do better — is a mess.
    
    bool valid = true;


    if(seq_input_.empty())
    {
        std::cout << "No sequence input provided for compacted de Bruijn graph construction.\n";
        valid = false;
    }


    // Check if read and reference de Bruijn graph parameters are being mixed with.
    if(is_read_graph_)  // Is a read de Bruijn graph.
    {
        if(output_format_ != cuttlefish::Output_Format::txt)
        {
            std::cout << "(Currently) Unsupported output file format requested for the compacted read de Bruijn graph.\n";
            valid = false;
        }
    }
    else    // Is a reference de Bruijn graph.
    {
        if(!edge_db_path_.empty())
        {
            std::cout << "No edge (i.e. (k + 1)-mer) database is required for a compacted reference de Bruijn graph construction.\n";
            valid = false;
        }
    }


    // Even `k` values are not consistent with the theory.
    // Also, `k` needs to be in the range `[1, MAX_K]`.
    if((k_ & 1U) == 0 || (k_ > cuttlefish::MAX_K))
    {
        std::cout << "The k-mer length (k) needs to be odd and within " << cuttlefish::MAX_K << ".\n";
        valid = false;
    }


    // Discard unsupported thread counts.
    const auto num_threads = std::thread::hardware_concurrency();
    if(num_threads > 0 && thread_count_ > num_threads)
    {
        std::cout << "At most " << num_threads << " concurrent threads are supported at the machine.\n";
        valid = false;
    }


    // Discard invalid output formats.
    if(output_format_ >= cuttlefish::num_op_formats)
    {
        std::cout << "Invalid output file format.\n";
        valid = false;
    }


    return valid;
}
