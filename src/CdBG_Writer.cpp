
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
    const std::string& output_file_path = params.output_file_path();

    // Open a parser for the FASTA / FASTQ file containing the reference.
    Parser parser(reference_input);


    // Clear the output file.
    clear_output_file();

    // Open an asynchronous logger to write into the output file, and set its log message pattern.
    cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
    output->set_pattern("%v");

    
    // Allocate output buffers for each thread.
    allocate_output_buffers();

    // Construct a thread pool.
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_plain);


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
        distribute_output_plain(seq, seq_len, output,thread_pool);
        thread_pool.wait_completion();
    }

    std::cerr << "\rProcessed " << seq_count << " sequences. Total reference length is " << ref_len << " bases.\n";
    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";


    // Close the thread-pool.
    thread_pool.close();


    // Flush the buffers.
    flush_output_buffers(output);

    // Clear and close the loggers.
    spdlog::drop_all();


    // Close the parser.
    parser.close();


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs (in plain text). Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::distribute_output_plain(const char* const seq, const size_t seq_len, cuttlefish::logger_t output, Thread_Pool<k>& thread_pool)
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
        thread_pool.assign_output_task(idle_thread_id, seq, seq_len, left_end, right_end, output);

        left_end += task_size;
    }
}


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs_gfa()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const Reference_Input& reference_input = params.reference_input();
    const uint16_t thread_count = params.thread_count();
    const std::string& output_file_path = params.output_file_path();
    const std::string& working_dir_path = params.working_dir_path();

    // Open a parser for the FASTA / FASTQ file containing the reference.
    Parser parser(reference_input);


    // Clear the output file and write the GFA header.
    std::ofstream op_stream(output_file_path.c_str(), std::ofstream::out);
    if(!op_stream)
    {
        std::cerr << "Error opening output file " << output_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    write_gfa_header(op_stream);
    op_stream.close();

    // Open an asynchronous logger to write into the output file, and set its log message pattern.
    // Note: `spdlog` appends to the output file by default, so the results for the sequences are accumulated into the same output file.
    // cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
    // output->set_pattern("%v");

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


        // Open an asynchronous logger to write into the output file, and set its log message pattern.
        // Note: `spdlog` appends to the output file by default, so the results for the sequences are accumulated into the same output file.
        cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
        output->set_pattern("%v");

        // Reset the path output streams for each thread.
        reset_path_streams();

        
        // Reset the first, the second, and the last unitigs seen for each thread.
        std::fill(first_unitig.begin(), first_unitig.end(), Oriented_Unitig());
        std::fill(second_unitig.begin(), second_unitig.end(), Oriented_Unitig());
        std::fill(last_unitig.begin(), last_unitig.end(), Oriented_Unitig());
        

        // Single-threaded writing.
        // output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        distribute_output_gfa(seq, seq_len, output, thread_pool);
        thread_pool.wait_completion();
        write_inter_thread_connections(output);

        
        // Flush all the buffered content (segments and links), as the GFA path to be appended to
        // the same output sink file is written in a different way than using the `spdlog` logger.
        // Note: If using an async logger, `logger->flush()` posts a message to the queue requesting the flush
        // operation, so the function returns immediately. Hence a forceful eviction is necessary by shutdown.
        std::chrono::high_resolution_clock::time_point t_s = std::chrono::high_resolution_clock::now();
        output->flush();
        spdlog::shutdown();
        std::chrono::high_resolution_clock::time_point t_e = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_e - t_s).count();
        logger_flush_time += elapsed_seconds;
        

        // Write the GFA path for this sequence.
        flush_path_buffers();
        params.output_format() == 1 ? write_gfa_path() : write_gfa_ordered_group();
    }

    std::cerr << "\rProcessed " << seq_count << " sequences. Total reference length is " << ref_len << " bases.\n";
    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";

    
    // Close the thread pool.
    thread_pool.close();


    // Flush the buffers.
    cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
    output->set_pattern("%v");
    flush_output_buffers(output);

    // Remove the temporary files.
    remove_temp_files();

    // Close the loggers?
    spdlog::drop_all();


    // Close the parser.
    parser.close();


    // Debug
    std::cout << "Segments writing time to in-memory buffers (avg): " <<
        std::accumulate(seg_write_time.begin(), seg_write_time.end(), 0) / thread_count << "\n";
    std::cout << "Links writing time to in-memory buffers (avg): " <<
        std::accumulate(link_write_time.begin(), link_write_time.end(), 0) / thread_count << "\n";
    std::cout << "Segments and links flush time to disk (avg): " <<
        std::accumulate(buff_flush_time.begin(), buff_flush_time.end(), 0) / thread_count << "\n";
    std::cout << "Paths writing time to in-memory buffers (avg): " <<
        std::accumulate(path_write_time.begin(), path_write_time.end(), 0) / thread_count << "\n";
    std::cout << "Paths flush time to disk (avg): " <<
        std::accumulate(path_flush_time.begin(), path_flush_time.end(), 0) / thread_count << "\n";
    std::cout << "Total paths concatenation time: " << path_concat_time << "\n";
    std::cout << "Total logger shutdown time: " << logger_flush_time << "\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs (in GFA format). Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::distribute_output_gfa(const char* const seq, const size_t seq_len, cuttlefish::logger_t output, Thread_Pool<k>& thread_pool)
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
        thread_pool.assign_output_task(task_id, seq, seq_len, left_end, right_end, output);

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
void CdBG<k>::check_output_buffer(const uint16_t thread_id, cuttlefish::logger_t output)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    if(output_buffer[thread_id].size() >= BUFFER_THRESHOLD)
    {
        // TODO: Avoid the presence of the same content in the memory simultaneously.
        // Done, by replacing `stringstream` with `string`.
        // E.g. skipping this line results in consuming memory:
        // "fixed" HG38: 176 instead of 346 MB;
        // HG38: 17 instead of 93 MB;
        // Gorgor3: 0.5 instead of 40MB.

        write(output, output_buffer[thread_id]);
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
void CdBG<k>::write(cuttlefish::logger_t output, const std::string& str)
{
    output->info("{}", str);
}


template <uint16_t k>
void CdBG<k>::write(std::ofstream& output, const std::string& str)
{
    output.write(str.c_str(), str.size());
}


template <uint16_t k>
void CdBG<k>::flush_output_buffers(cuttlefish::logger_t output)
{
    const uint16_t thread_count = params.thread_count();

    for (uint16_t t_id = 0; t_id < thread_count; ++t_id)
        if(!output_buffer[t_id].empty())
        {
            std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

            write(output, output_buffer[t_id]);
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
