
#include "CdBG.hpp"
#include "DNA_Utility.hpp"
#include "Job_Queue.hpp"
#include "fmt/format.h"


template <uint16_t k>
void CdBG<k>::write_segment(const uint16_t thread_id, const char* const seq, const uint64_t segment_name, const size_t start_kmer_idx, const size_t end_kmer_idx, const cuttlefish::dir_t dir)
{
    std::string& buffer = output_buffer[thread_id];
    const size_t segment_len = end_kmer_idx - start_kmer_idx + k;

    
    ensure_buffer_space(buffer, segment_len + 22, output_[thread_id]);
    
    
    // The 'Name' field.
    buffer += fmt::format_int(segment_name).c_str();
    
    // The segment field.
    buffer += "\t";
    if(dir == cuttlefish::FWD)
        for(size_t offset = 0; offset < segment_len; ++offset)
            buffer += DNA_Utility::upper(seq[start_kmer_idx + offset]);
    else
        for(size_t offset = 0; offset < segment_len; ++offset)
            buffer += DNA_Utility::complement(seq[end_kmer_idx + k - 1 - offset]);


    // End the segment line.
    buffer += "\n";


    if(kmer_idx != nullptr)
    {
        const char* const unitig_seq = buffer.data() + (buffer.size() - (segment_len + 1)); // The +1 length is to account for the ending line-break.
        kmer_idx->deposit(token[thread_id], unitig_seq, segment_len);
    }


    // Mark buffer size increment.
    check_output_buffer(thread_id);
}


template <uint16_t k>
void CdBG<k>::write_sequence_tiling(Job_Queue<std::string, Oriented_Unitig>& job_queue)
{
    const uint16_t thread_count = params.thread_count();
    const std::string& seq_file_path = params.sequence_file_path();
    const bool poly_n_stretch = params.poly_n_stretch();

    // Open the output file in append mode.
    std::ofstream output(seq_file_path.c_str(), std::ios_base::app);


    while(true)
    {
        // Busy-wait for jobs.
        while(!job_queue.job_available())
            if(!job_queue.jobs_remain())
            {
                output.close();

                return;
            }


        // Fetch the next tiling ID, and its very first GFA edge in the sequence (it is not inferrable from the tilings).
        std::string path_id;
        Oriented_Unitig left_unitig;
        
        job_queue.fetch_job(path_id, left_unitig);

        // The sequence does not contain any unitig (possible if there's no valid k-mer in the sequence).
        if(!left_unitig.is_valid())
        {
            remove_temp_files(job_queue.next_job_to_finish());
            job_queue.finish_job();
            continue;
        }

        // Write the path ID.
        output << path_id;

        // Write the path members.
        output << "\t";

        if(poly_n_stretch && left_unitig.start_kmer_idx > 0)
            output << "N" << left_unitig.start_kmer_idx << " ";
        
        // The first vertex of the path (not inferrable from the path output files).
        output << left_unitig.unitig_id << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");

        // Copy the thread-specific path output file contents to the sequence-tiling file.
        for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        {
            const std::string path_file_name_(path_file_name(t_id, job_queue.next_job_to_finish()));
            std::ifstream input(path_file_name_.c_str(), std::ios_base::in);

            if(input.fail())
            {
                std::cerr << "Error opening temporary path output file " << path_file_name_ << ". Aborting.\n";
                std::exit(EXIT_FAILURE);
            }
            
            // Copy the path output for thread number `t_id` to the end of the sequence-tiling file.
            if(input.peek() != EOF)
                output << input.rdbuf();

            input.close();
        }

        // End the path.
        output << "\n";

        // Remove the thread-specific path output files (for this tiling job).
        remove_temp_files(job_queue.next_job_to_finish());
        

        // Mark the tiling as completed.
        job_queue.finish_job();
    }
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
