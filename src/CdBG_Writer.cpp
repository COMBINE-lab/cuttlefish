
#include "CdBG.hpp"
#include "Parser.hpp"
#include "Output_Format.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <iomanip>
#include <numeric>


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs()
{
    const uint8_t output_format = params.output_format();

    if(output_format == cuttlefish::txt)
        output_maximal_unitigs_plain();
    else
        output_maximal_unitigs_gfa();
}


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs_plain()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    
    const Reference_Input& reference_input = params.reference_input();
    const uint16_t thread_count = params.thread_count();
    // const std::string& output_file_path = params.output_file_path();

    // Open a parser for the FASTA / FASTQ file containing the reference.
    Parser parser(reference_input);


    // Clear the output file and initialize the output loggers.
    clear_output_file();
    init_output_loggers(); 
    
    // Allocate output buffers for each thread.
    allocate_output_buffers();


    // Construct a thread pool.
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_plain);

    // Debug
    seg_write_time.resize(thread_count);
    buff_flush_time.resize(thread_count);


    // Track the maximum sequence buffer size used.
    size_t max_buf_sz = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    uint32_t seq_count = 0;
    uint64_t ref_len = 0;
    while(parser.read_next_seq())
    {
        const char* const seq = parser.seq();
        const size_t seq_len = parser.seq_len();
        const size_t seq_buf_sz = parser.buff_sz();

        seq_count++;
        ref_len += seq_len;
        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);
        std::cerr << "\rProcessing sequence " << seq_count << ", with length " << std::setw(10) << seq_len << ".";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Single-threaded writing.
        // output_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        distribute_output_plain(seq, seq_len, thread_pool);
        thread_pool.wait_completion();
    }

    std::cerr << "\rProcessed " << seq_count << " sequences. Total reference length is " << ref_len << " bases.\n";
    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";


    // Close the thread-pool.
    thread_pool.close();


    // Flush the buffers.
    flush_output_buffers();

    // Clear and close the loggers.
    spdlog::drop_all();


    // Close the parser.
    parser.close();


    // Debug
    std::cout << "Unitigs writing time to in-memory buffers (total): " <<
        std::accumulate(seg_write_time.begin(), seg_write_time.end(), 0) << "\n";
    std::cout << "Unitigs flush time to disk (total): " <<
        std::accumulate(buff_flush_time.begin(), buff_flush_time.end(), 0) << "\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs (in plain text). Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::distribute_output_plain(const char* const seq, const size_t seq_len, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();
    const size_t task_size = (seq_len - k + 1) / thread_count;
    const uint16_t partition_count = (task_size < PARTITION_SIZE_THRESHOLD ? 1 : thread_count);

    size_t left_end = 0;
    size_t right_end;

    for(uint16_t task_id = 0; task_id < partition_count; ++task_id)
    {
        right_end = (task_id == partition_count - 1 ? seq_len - k : left_end + task_size - 1);
        
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_output_task(idle_thread_id, seq, seq_len, left_end, right_end);

        left_end += task_size;
    }
}


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs_gfa()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const Reference_Input& reference_input = params.reference_input();
    const uint16_t thread_count = params.thread_count();
    // const std::string& output_file_path = params.output_file_path();
    const std::string& working_dir_path = params.working_dir_path();


    // Open a parser for the FASTA / FASTQ file containing the reference.
    Parser parser(reference_input);


    // Clear the output file and write the GFA header.
    clear_output_file();
    write_gfa_header();

    // Set the prefixes of the temporary path output files. This is to avoid possible name
    // conflicts in the file system.
    set_temp_file_prefixes(working_dir_path);


    // Allocate the output buffers for each thread.
    allocate_output_buffers();

    // Allocate the path output buffers for each thread.
    allocate_path_buffers();

    // Allocate entries for the first, the second, and the last unitigs seen by each thread.
    first_unitig.resize(thread_count);
    second_unitig.resize(thread_count);
    last_unitig.resize(thread_count);


    // Construct a thread pool.
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_gfa);

    // Debug
    seg_write_time.resize(thread_count);
    link_write_time.resize(thread_count);
    buff_flush_time.resize(thread_count);
    path_write_time.resize(thread_count);
    path_flush_time.resize(thread_count);


    // Track the maximum sequence buffer size used.
    size_t max_buf_sz = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    // uint32_t seq_count = 0;
    // TODO: replace usage of `seq_count` with a parser method.
    seq_count = 0;
    uint64_t ref_len = 0;
    while(parser.read_next_seq())
    {
        const char* const seq = parser.seq();
        const size_t seq_len = parser.seq_len();
        const size_t seq_buf_sz = parser.buff_sz();


        seq_count++;
        ref_len += seq_len;
        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);
        std::cerr << "\rProcessing sequence " << seq_count << ", with length " << std::setw(10) << seq_len << ".";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Initialize the output loggers.
        // Note: `spdlog` appends to the output file by default, so the results are accumulated into the
        // same output file; thus repeated initializations do not result in clearing out the output file.
        init_output_loggers();

        // Reset the path output streams for each thread.
        reset_path_streams();

        // Reset the first, the second, and the last unitigs seen for each thread.
        std::fill(first_unitig.begin(), first_unitig.end(), Oriented_Unitig());
        std::fill(second_unitig.begin(), second_unitig.end(), Oriented_Unitig());
        std::fill(last_unitig.begin(), last_unitig.end(), Oriented_Unitig());
        

        // Single-threaded writing.
        // output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        distribute_output_gfa(seq, seq_len, thread_pool);
        thread_pool.wait_completion();
        write_inter_thread_connections();

        
        // Flush all the buffered content (segments and links), as the GFA path to be appended to
        // the same output sink file is written in a different way than using the `spdlog` logger.
        std::chrono::high_resolution_clock::time_point t_s = std::chrono::high_resolution_clock::now();
        close_output_loggers();
        std::chrono::high_resolution_clock::time_point t_e = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_e - t_s).count();
        logger_flush_time += elapsed_seconds;
        

        // Write the GFA path for this sequence.
        // flush_path_buffers();
        // params.output_format() == 1 ? write_gfa_path() : write_gfa_ordered_group();
    }

    std::cerr << "\rProcessed " << seq_count << " sequences. Total reference length is " << ref_len << " bases.\n";
    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";

    
    // Close the thread pool.
    thread_pool.close();


    // Flush the buffers.
    init_output_loggers();
    flush_output_buffers();

    // Remove the temporary files.
    remove_temp_files();

    // Close `spdlog`.
    spdlog::drop_all();


    // Close the parser.
    parser.close();


    // Debug
    std::cout << "Segments writing time to in-memory buffers (total): " <<
        std::accumulate(seg_write_time.begin(), seg_write_time.end(), 0) << "\n";
    std::cout << "Links writing time to in-memory buffers (total): " <<
        std::accumulate(link_write_time.begin(), link_write_time.end(), 0) << "\n";
    std::cout << "Segments and links flush time to disk (total): " <<
        std::accumulate(buff_flush_time.begin(), buff_flush_time.end(), 0) << "\n";
    std::cout << "Paths writing time to in-memory buffers (total): " <<
        std::accumulate(path_write_time.begin(), path_write_time.end(), 0) << "\n";
    std::cout << "Paths flush time to disk (total): " <<
        std::accumulate(path_flush_time.begin(), path_flush_time.end(), 0) << "\n";
    std::cout << "Total paths concatenation time: " << path_concat_time << "\n";
    std::cout << "Total logger shutdown time: " << logger_flush_time << "\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs (in GFA format). Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::distribute_output_gfa(const char* const seq, const size_t seq_len, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();
    const size_t task_size = (seq_len - k + 1) / thread_count;
    const uint16_t partition_count = (task_size < PARTITION_SIZE_THRESHOLD ? 1 : thread_count);

    size_t left_end = 0;
    size_t right_end;

    for(uint16_t task_id = 0; task_id < partition_count; ++task_id)
    {
        right_end = (task_id == partition_count - 1 ? seq_len - k : left_end + task_size - 1);
        
        thread_pool.get_thread(task_id);
        thread_pool.assign_output_task(task_id, seq, seq_len, left_end, right_end);

        left_end += task_size;
    }
}


template <uint16_t k>
void CdBG<k>::clear_output_file() const
{
    const std::string& output_file_path = params.output_file_path();

    std::ofstream output(output_file_path.c_str(), std::ofstream::out);
    if(!output)
    {
        std::cerr << "Error opening output file " << output_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    output.close();
}


template <uint16_t k>
void CdBG<k>::init_output_loggers()
{
    const std::string& output_file_path = params.output_file_path();
    const uint16_t thread_count = params.thread_count();

    // Open an asynchronous logger to write into the output file, and set its log message pattern.
    output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
    output->set_pattern("%v");

    output_.resize(thread_count);
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        output_[t_id] = output;
}


template <uint16_t k>
void CdBG<k>::close_output_loggers()
{
    // Note: If using an async logger, `logger->flush()` posts a message to the queue requesting the flush
    // operation, so the function returns immediately. Hence a forceful eviction is necessary by shutdown.
    output->flush();
    spdlog::shutdown();
}


template <uint16_t k>
void CdBG<k>::allocate_output_buffers()
{
    const uint16_t thread_count = params.thread_count();

    output_buffer.resize(thread_count);
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        output_buffer[t_id].reserve(BUFFER_CAPACITY);
}


template <uint16_t k>
void CdBG<k>::allocate_path_buffers()
{
    const uint16_t thread_count = params.thread_count();
    const uint8_t gfa_v = params.output_format();


    path_buffer.resize(thread_count);
    if(gfa_v == cuttlefish::gfa1)
        overlap_buffer.resize(thread_count);

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        path_buffer[t_id].reserve(BUFFER_CAPACITY);
        if(gfa_v == cuttlefish::gfa1)
            overlap_buffer[t_id].reserve(BUFFER_CAPACITY);
    }
}


template <uint16_t k>
bool CdBG<k>::is_unipath_start(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::dir_t dir, const cuttlefish::Vertex_Class prev_kmer_class, const cuttlefish::dir_t prev_kmer_dir) const
{
    if(vertex_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(vertex_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;


    // assert(kmer_idx > 0);


    if(prev_kmer_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD)
    {
        if(prev_kmer_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;
    }
    else    // prev_kmer_dir == cuttlefish::BWD
        if(prev_kmer_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;

    
    return false;
}


template <uint16_t k>
bool CdBG<k>::is_unipath_end(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::dir_t dir, const cuttlefish::Vertex_Class next_kmer_class, const cuttlefish::dir_t next_kmer_dir) const
{
    if(vertex_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(vertex_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;


    // assert(kmer_idx < ref.length() - k);


    if(next_kmer_class == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;

    if(next_kmer_dir == cuttlefish::FWD)
    {
        if(next_kmer_class == cuttlefish::Vertex_Class::multi_in_single_out)
            return true;
    }
    else    // next_kmer_dir == cuttlefish::BWD
        if(next_kmer_class == cuttlefish::Vertex_Class::single_in_multi_out)
            return true;


    return false;
}


template <uint16_t k>
void CdBG<k>::check_output_buffer(const uint16_t thread_id)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    if(output_buffer[thread_id].size() >= BUFFER_THRESHOLD)
    {
        write(thread_id, output_buffer[thread_id]);
        output_buffer[thread_id].clear();
    }

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    buff_flush_time[thread_id] += elapsed_seconds;
}


template <uint16_t k>
void CdBG<k>::check_path_buffer(const uint16_t thread_id)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    if(path_buffer[thread_id].size() >= BUFFER_THRESHOLD)
    {
        write(path_output[thread_id], path_buffer[thread_id]);
        path_buffer[thread_id].clear();
    }


    const uint8_t gfa_v = params.output_format();
    if(gfa_v == cuttlefish::gfa1 && overlap_buffer[thread_id].size() >= BUFFER_THRESHOLD)
    {
        write(overlap_output[thread_id], overlap_buffer[thread_id]);
        overlap_buffer[thread_id].clear();
    }

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    path_flush_time[thread_id] += elapsed_seconds;
}


template <uint16_t k>
void CdBG<k>::write(const uint16_t thread_id, const std::string& str)
{
    output_[thread_id]->info("{}", str);
}


template <uint16_t k>
void CdBG<k>::write(std::ofstream& op, const std::string& str)
{
    op.write(str.c_str(), str.size());
}


template <uint16_t k>
void CdBG<k>::flush_output_buffers()
{
    const uint16_t thread_count = params.thread_count();

    for (uint16_t t_id = 0; t_id < thread_count; ++t_id)
        if(!output_buffer[t_id].empty())
        {
            std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

            write(t_id, output_buffer[t_id]);
            output_buffer[t_id].clear();

            std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
            double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

            buff_flush_time[t_id] += elapsed_seconds;
        }
}


template <uint16_t k>
void CdBG<k>::flush_path_buffers()
{
    const uint16_t thread_count = params.thread_count();
    const uint8_t gfa_v = params.output_format();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

        if(!path_buffer[t_id].empty())
        {
            write(path_output[t_id], path_buffer[t_id]);
            path_buffer[t_id].clear();
        }

        if(gfa_v == cuttlefish::gfa1 && !overlap_buffer[t_id].empty())
        {
            write(overlap_output[t_id], overlap_buffer[t_id]);
            overlap_buffer[t_id].clear();
        }

        std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

        path_flush_time[t_id] += elapsed_seconds;
    }
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
