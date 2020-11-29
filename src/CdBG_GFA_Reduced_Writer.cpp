
#include "CdBG.hpp"



template <uint16_t k>
void CdBG<k>::write_segment(const uint16_t thread_id, const char* const seq, const uint64_t segment_name, const size_t start_kmer_idx, const size_t end_kmer_idx, const cuttlefish::dir_t dir)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    std::string& buffer = output_buffer[thread_id];
    const size_t segment_len = end_kmer_idx - start_kmer_idx + k;

    
    ensure_buffer_space(buffer, segment_len + 22, output_[thread_id]);
    
    
    // The 'Name' field.
    buffer += fmt::format_int(segment_name).c_str();
    
    // The segment field.
    buffer += "\t";
    if(dir == cuttlefish::FWD)
        for(size_t offset = 0; offset < segment_len; ++offset)
            buffer += Kmer<k>::upper(seq[start_kmer_idx + offset]);
    else
        for(size_t offset = 0; offset < segment_len; ++offset)
            buffer += Kmer<k>::complement(seq[end_kmer_idx + k - 1 - offset]);


    // End the segment line.
    buffer += "\n";

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    seg_write_time[thread_id] += elapsed_seconds;


    // Mark buffer size increment.
    check_output_buffer(thread_id);
}


template <uint16_t k>
void CdBG<k>::write_sequence_tiling(const std::string& path_id)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const uint16_t thread_count = params.thread_count();
    const std::string& seq_file_path = params.output_file_path() + SEQ_FILE_EXT;
    

    // Search the very first GFA edge in the sequence, as that is not inferrable from the path outputs.
    Oriented_Unitig left_unitig;
    Oriented_Unitig right_unitig;
    search_first_connection(left_unitig, right_unitig);

    // The sequence does not contain any unitig (possible if there's no valid k-mer in the sequence).
    if(!left_unitig.is_valid())
        return;

    
    // Open the output file in append mode.
    std::ofstream output(seq_file_path.c_str(), std::ios_base::app);

    // Write the path ID.
    output << path_id;

    // Write the path members.
    output << "\t";
    
    // The first vertex of the path (not inferrable from the path output files).
    output << left_unitig.unitig_id << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // Copy the thread-specific path output file contents to the sequence-tiling file.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const std::string path_file_name = (path_file_prefix + std::to_string(t_id));
        std::ifstream input(path_file_name.c_str(), std::ios_base::in);

        if(input.fail())
        {
            std::cerr << "Error opening temporary path output file " << path_file_name << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
        
        // Copy the path output for thread number `t_id` to the end of the sequence-tiling file.
        if(input.peek() != EOF)
            output << input.rdbuf();

        input.close();
    }


    // End the path.
    output << "\n";

    output.close();

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    (void)elapsed_seconds;
    // std::cout << "Time taken to write paths = " << elapsed_seconds << " seconds.\n";
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
