
#include "CdBG.hpp"
#include "Ref_Parser.hpp"
#include "Output_Format.hpp"
#include "Thread_Pool.hpp"
#include "Job_Queue.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <iomanip>


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs()
{
    const uint8_t output_format = params.output_format();
    unipaths_info_local.resize(params.thread_count());

    if(output_format == cuttlefish::fa)
        output_maximal_unitigs_plain();
    else if(output_format == cuttlefish::gfa1 || output_format == cuttlefish::gfa2)
        output_maximal_unitigs_gfa();
    else if(output_format == cuttlefish::gfa_reduced)
        output_maximal_unitigs_gfa_reduced();

    
    for(uint16_t t_id = 0; t_id < params.thread_count(); ++t_id)
        unipaths_meta_info_.aggregate(unipaths_info_local[t_id]);

    dbg_info.add_unipaths_info(*this);
}


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs_plain()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    
    const Seq_Input& reference_input = params.sequence_input();
    const uint16_t thread_count = params.thread_count();

    // Open a parser for the FASTA / FASTQ file containing the reference.
    Ref_Parser parser(reference_input);


    // Clear the output file and initialize the output loggers.
    clear_output_file();
    init_output_loggers(); 
    
    // Allocate output buffers for each thread.
    allocate_output_buffers();


    // Construct a thread pool.
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_plain);


    // Track the maximum sequence buffer size used and the total length of the references.
    size_t max_buf_sz = 0;
    uint64_t ref_len = 0;
    uint64_t seq_count = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    while(parser.read_next_seq())
    {
        const char* const seq = parser.seq();
        const size_t seq_len = parser.seq_len();
        const size_t seq_buf_sz = parser.buff_sz();

        seq_count++;
        ref_len += seq_len;
        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);
        std::cerr << "\rProcessing sequence " << parser.seq_id() << ", with length:\t" << std::setw(10) << seq_len << ".";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Single-threaded writing.
        // output_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        distribute_output_plain(seq, seq_len, thread_pool);
        thread_pool.wait_completion();
    }

    std::cout << "\nProcessed " << seq_count << " sequences. Total reference length: " << ref_len << " bases.\n";
    std::cout << "Maximum input sequence buffer size used: " << max_buf_sz / (1024 * 1024) << " MB.\n";


    // Close the thread-pool.
    thread_pool.close();


    // Flush the buffers.
    flush_output_buffers();

    // Clear `spdlog`.
    spdlog::drop_all();


    // Close the parser.
    parser.close();


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done writing the maximal unitigs (in plain text). Time taken = " << elapsed_seconds << " seconds.\n";
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


    const Seq_Input& reference_input = params.sequence_input();
    const uint16_t thread_count = params.thread_count();
    const std::string& working_dir_path = params.working_dir_path();



    // Clear the output file and write the GFA header.
    clear_output_file();
    write_gfa_header();

    // Set the prefixes of the temporary path output files. This is to avoid possible name
    // conflicts in the file system.
    set_temp_file_prefixes(working_dir_path);


    // Allocate the output buffers and the path buffers for each thread.
    allocate_output_buffers();
    allocate_path_buffers();


    // Construct a thread pool.
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_gfa);


    // Open a parser for the FASTA / FASTQ file containing the reference.
    Ref_Parser parser(reference_input);

    // Track the maximum sequence buffer size used and the total length of the references.
    size_t max_buf_sz = 0;
    uint64_t ref_len = 0;
    uint64_t seq_count = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    while(parser.read_next_seq())
    {
        const char* const seq = parser.seq();
        const size_t seq_len = parser.seq_len();
        const size_t seq_buf_sz = parser.buff_sz();


        seq_count++;
        ref_len += seq_len;
        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);
        std::cerr << "\rProcessing sequence " << parser.seq_id() << ", with length:\t" << std::setw(10) << seq_len << ".";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Initialize the output loggers.
        // Note: `spdlog` appends to the output file by default, so the results are accumulated into the
        // same output file; thus repeated initializations do not result in clearing out the output file.
        init_output_loggers();

        // Reset the path output streams for each thread.
        reset_path_loggers();

        // Reset the first, the second, and the last unitigs seen for each thread.
        reset_extreme_unitigs();
        

        // Single-threaded writing.
        // output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        distribute_output_gfa(seq, seq_len, thread_pool);
        thread_pool.wait_completion();

        write_inter_thread_connections();

        flush_path_buffers();

        // Force flush all the in-memory logs (segments, connections, and broken paths), as the GFA path to be
        // appended to the same output sink file is written using a different mechanism (copy with `rdbuf()`)
        // than using the `spdlog` logger.
        close_loggers();

        // Write the GFA path for this sequence.
        const std::string path_name =   std::string("Reference:") + std::to_string(parser.ref_id()) +
                                        std::string("_Sequence:") + remove_whitespaces(parser.seq_name());
        params.output_format() == 1 ? write_gfa_path(path_name) : write_gfa_ordered_group(path_name);
    }

    std::cout << "\nProcessed " << seq_count << " sequences. Total reference length: " << ref_len << " bases.\n";
    std::cout << "Maximum input sequence buffer size used: " << max_buf_sz / (1024 * 1024) << " MB.\n";

    
    // Close the thread pool.
    thread_pool.close();


    // Flush the buffers.
    init_output_loggers();
    flush_output_buffers();

    // Close `spdlog`.
    spdlog::drop_all();

    // Remove the temporary files.
    remove_temp_files();

    // Close the parser.
    parser.close();


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout <<    "Done writing the compacted graph"
                    " (in GFA " << (params.output_format() == cuttlefish::Output_Format::gfa1 ? 1 : 2) << " format)."
                    " Time taken = " << elapsed_seconds << " seconds.\n";
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
void CdBG<k>::output_maximal_unitigs_gfa_reduced()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const Seq_Input& reference_input = params.sequence_input();
    const uint16_t thread_count = params.thread_count();
    const std::string& working_dir_path = params.working_dir_path();


    // Clear the output file and initilize the output loggers.
    clear_output_file();
    init_output_loggers();

    // Set the prefixes of the temporary path output files. This is to avoid possible name
    // conflicts in the file system.
    set_temp_file_prefixes(working_dir_path);


    // Allocate the output buffers and the path buffers for each thread.
    allocate_output_buffers();
    allocate_path_buffers();


    // Construct a thread pool.
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_gfa_reduced);

    // Dedicated thread and job-queue to concatenate thread-specific tilings.
    std::unique_ptr<std::thread> concatenator{nullptr};
    Job_Queue<std::string, Oriented_Unitig> job_queue;

    // Launch the background tilings-concatenator thread.
    concatenator.reset(
        new std::thread([this, &job_queue]()
        {
            write_sequence_tiling(job_queue);
        })
    );


    // Open a parser for the FASTA / FASTQ file containing the reference.
    Ref_Parser parser(reference_input);

    // Track the maximum sequence buffer size used and the total length of the references.
    size_t max_buf_sz = 0;
    uint64_t ref_len = 0;
    uint64_t seq_count = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    while(parser.read_next_seq())
    {
        const char* const seq = parser.seq();
        const size_t seq_len = parser.seq_len();
        const size_t seq_buf_sz = parser.buff_sz();


        seq_count++;
        ref_len += seq_len;
        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);
        std::cerr << "\rProcessing sequence " << parser.seq_id() << ", with length:\t" << std::setw(10) << seq_len << ".";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Reset the path output streams for each thread.
        reset_path_loggers(job_queue.next_job_to_post());

        // Reset the first, the second, and the last unitigs seen for each thread.
        reset_extreme_unitigs();
        

        // Single-threaded writing.
        // output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        distribute_output_gfa(seq, seq_len, thread_pool);
        thread_pool.wait_completion();

        write_inter_thread_connections();

        flush_path_buffers();

        // Force-flush all the path loggers, as the thread-specific files for the path outputs are to be reused
        // for the next sequence. A problem with having a new set of files for the next sequence (and thus avoiding
        // this force-flush) is that an input reference collection may have millions of sequences (e.g. the conifers),
        // thus exploding the limits of the underlying file system.
        close_path_loggers();

        
        // Write the GFA path for this sequence.
        const std::string path_name =   std::string("Reference:") + std::to_string(parser.ref_id()) +
                                        std::string("_Sequence:") + remove_whitespaces(parser.seq_name());
        

        // Post a tiling-concatenation job.
        
        // Search the very first GFA edge in the sequence, as that is not inferrable from the path outputs.
        Oriented_Unitig left_unitig, right_unitig;
        search_first_connection(left_unitig, right_unitig);

        job_queue.post_job(path_name, left_unitig);
    }

    std::cout << "\nProcessed " << seq_count << " sequences. Total reference length: " << ref_len << " bases.\n";
    std::cout << "Maximum input sequence buffer size used: " << max_buf_sz / (1024 * 1024) << " MB.\n";

    
    // Close the thread pool.
    thread_pool.close();

    // Signal the end of the tiling concatenations.
    job_queue.signal_end();
    if(!concatenator->joinable())
    {
        std::cerr << "Early termination encountered for the sequence-tilings concatenator thread. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    concatenator->join();


    // Flush the buffers.
    flush_output_buffers();

    // Close `spdlog`.
    spdlog::drop_all();

    // Close the parser.
    parser.close();


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done writing the compacted graph (in GFA-reduced format). Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::clear_output_file() const
{
    const cuttlefish::Output_Format op_format = params.output_format();

    if(op_format == cuttlefish::fa || op_format == cuttlefish::gfa1 || op_format == cuttlefish::gfa2)
        clear_file(params.output_file_path());
    else if(op_format == cuttlefish::gfa_reduced)
    {
        const std::string seg_file_path(params.segment_file_path());
        const std::string seq_file_path(params.sequence_file_path());

        clear_file(seg_file_path);
        clear_file(seq_file_path);
    }
}


template <uint16_t k>
void CdBG<k>::init_output_loggers()
{
    const cuttlefish::Output_Format gfa_v = params.output_format();
    const std::string& output_file_path = (gfa_v == cuttlefish::Output_Format::gfa_reduced ?
                                            params.segment_file_path() : params.output_file_path());
    const uint16_t thread_count = params.thread_count();


    // Initialize the global thread-pool of `spdlog` with restricting its default queue size from 8192 to
    // a suitable one for us. This is required to restrict the memory-usage during the output step, as
    // `spdlog` can continue on accumulating logs (each of max length `BUFFER_THRESHOLD`), i.e. `spdlog`
    // itself can take up memory up-to (`BUFFER_THRESHOLD x #pending_log_message`), besides our own buffer
    // memory of (`BUFFER_CAPACITY x #output_threads`).
    
    // NB: Possibly, `spdlog::shutdown()` resets the global `spdlog` thread pool — so if using the global
    // thread pool, the following initializer needs to be invoked repeatedly after each shutdown.
    // In the current GFA (and -reduced) scheme, the segments and edges logger and the thread-specific
    // path loggers use separate thread pool instances, so we do not drop (and initialize) the global one.
    // The global pool was used in an earlier implementation.
    // spdlog::init_thread_pool(ASYNC_LOG_QUEUE_SZ, ASYNC_LOG_N_THREADS);

    // Instantiate a `spdlog` thread pool for outputting unitigs (segments) and links / edges (for GFA).
    tp_output = std::make_shared<spdlog::details::thread_pool>(ASYNC_LOG_QUEUE_SZ, ASYNC_LOG_N_THREADS);

    // Open an asynchronous logger to write into the output file, and set its log message pattern.
    // output = spdlog::basic_logger_mt<spdlog::async_factory>("async_output_logger", output_file_path);
    auto sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(output_file_path);
    output = std::make_shared<spdlog::async_logger>("async_output_logger", sink, tp_output, spdlog::async_overflow_policy::block);
    output->set_pattern("%v");

    output_.resize(thread_count);
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        output_[t_id] = output;
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
bool CdBG<k>::is_unipath_start(const cuttlefish::State_Class state_class, const cuttlefish::dir_t dir, const cuttlefish::State_Class prev_kmer_class, const cuttlefish::dir_t prev_kmer_dir) const
{
    if(state_class == cuttlefish::State_Class::multi_in_multi_out)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(state_class == cuttlefish::State_Class::multi_in_single_out)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(state_class == cuttlefish::State_Class::single_in_multi_out)
            return true;


    // assert(kmer_idx > 0);


    if(prev_kmer_class == cuttlefish::State_Class::multi_in_multi_out)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD)
    {
        if(prev_kmer_class == cuttlefish::State_Class::single_in_multi_out)
            return true;
    }
    else    // prev_kmer_dir == cuttlefish::BWD
        if(prev_kmer_class == cuttlefish::State_Class::multi_in_single_out)
            return true;

    
    return false;
}


template <uint16_t k>
bool CdBG<k>::is_unipath_end(const cuttlefish::State_Class state_class, const cuttlefish::dir_t dir, const cuttlefish::State_Class next_kmer_class, const cuttlefish::dir_t next_kmer_dir) const
{
    if(state_class == cuttlefish::State_Class::multi_in_multi_out)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(state_class == cuttlefish::State_Class::single_in_multi_out)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(state_class == cuttlefish::State_Class::multi_in_single_out)
            return true;


    // assert(kmer_idx < ref.length() - k);


    if(next_kmer_class == cuttlefish::State_Class::multi_in_multi_out)
        return true;

    if(next_kmer_dir == cuttlefish::FWD)
    {
        if(next_kmer_class == cuttlefish::State_Class::multi_in_single_out)
            return true;
    }
    else    // next_kmer_dir == cuttlefish::BWD
        if(next_kmer_class == cuttlefish::State_Class::single_in_multi_out)
            return true;


    return false;
}


template <uint16_t k>
void CdBG<k>::ensure_buffer_space(std::string& buf, const size_t log_len, const cuttlefish::logger_t& log)
{
    if(buf.size() + log_len >= BUFFER_CAPACITY - 1)
        flush_buffer(buf, log);
}


template <uint16_t k>
void CdBG<k>::flush_buffer(std::string& str, const cuttlefish::logger_t& log)
{
    write(str, log);

    str.clear();
}


template <uint16_t k>
void CdBG<k>::write(const std::string& str, const cuttlefish::logger_t& log)
{
    log->info("{}", str);
}


template <uint16_t k>
void CdBG<k>::check_output_buffer(const uint16_t thread_id)
{
    if(output_buffer[thread_id].size() >= BUFFER_THRESHOLD)
        flush_buffer(output_buffer[thread_id], output_[thread_id]);
}


template <uint16_t k>
void CdBG<k>::flush_output_buffers()
{
    const uint16_t thread_count = params.thread_count();

    for (uint16_t t_id = 0; t_id < thread_count; ++t_id)
        if(!output_buffer[t_id].empty())
            flush_buffer(output_buffer[t_id], output_[t_id]);
}


template <uint16_t k>
void CdBG<k>::close_loggers()
{
    flush_loggers();

    // Note: If using an async logger, `logger->flush()` posts a message to the queue requesting the flush
    // operation, so the function returns immediately. Hence a forceful eviction is necessary by shutdown.
    spdlog::shutdown();

    // Drop the `spdlog` thread pools. Required as `shutdown()` force-flushes messages from the global pool only.
    tp_output.reset();
    tp_path.reset();
}


template <uint16_t k>
void CdBG<k>::close_path_loggers()
{
    flush_path_loggers();

    // Note: If using an async logger, `logger->flush()` posts a message to the queue requesting the flush
    // operation, so the function returns immediately. Hence a forceful eviction is necessary by shutdown.
    // The following shuts down the global thread pool — which is no longer used in the current implementation.
    // spdlog::shutdown();

    // Drop the `spdlog` thread pools for path-outputs to force-flush the pending logs. Required as `shutdown()`
    // force-flushes messages from the global pool only.
    tp_path.reset();
}


template <uint16_t k>
void CdBG<k>::flush_loggers()
{
    flush_output_logger();

    flush_path_loggers();
}


template <uint16_t k>
void CdBG<k>::flush_output_logger()
{
    output->flush();
}


template<uint16_t k>
void CdBG<k>::flush_path_loggers()
{
    const uint16_t thread_count = params.thread_count();
    const uint8_t gfa_v = params.output_format();
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        path_output_[t_id]->flush();
        if(gfa_v == cuttlefish::gfa1)
            overlap_output_[t_id]->flush();
    }
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
