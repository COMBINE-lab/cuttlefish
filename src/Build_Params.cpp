
#include "Build_Params.hpp"
#include "utility.hpp"


bool Build_Params::is_valid() const
{
    // TODO: do better — is a mess.
    
    bool valid = true;


    // Check if read and reference de Bruijn graph parameters are being mixed with.
    if(is_read_graph_)  // Is a read de Bruijn graph.
    {
        if(!seq_input_.empty())
        {
            std::cout << "No reference is to be provided for a compacted read de Bruijn graph construction.\n";
            valid = false;
        }

        if(edge_db_path_.empty())
        {
            std::cout << "The path prefix to the KMC-database for edges (i.e. (k + 1)-mers) is required.\n";
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

        // if(dcc_opt_ || extract_cycles_)
        // {
        //     std::cout << "Existence of detached chordless cycles are impossible for reference de Bruijn graphs by definition.\n";
        //     valid = false;
        // }
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
