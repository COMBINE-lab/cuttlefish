
#include "CdBG.hpp"
#include "kseq/kseq.h"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs()
{
    const uint8_t output_format = params.output_format();

    // TODO: enum.
    if(output_format == 0)
        output_maximal_unitigs_plain();
    else
        output_maximal_unitigs_gfa();
}


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs_plain()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    
    const std::string& ref_file_path = params.ref_file_path();
    const uint16_t thread_count = params.thread_count();
    const std::string& output_file_path = params.output_file_path();

    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* const input = fopen(ref_file_path.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Initialize the parser.
    kseq_t* const parser = kseq_init(fileno(input));


    // Clear the output file.
    std::ofstream op_stream(output_file_path.c_str(), std::ofstream::out);
    if(!op_stream)
    {
        std::cerr << "Error opening output file " << output_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    op_stream.close();


    // Open an asynchronous logger to write into the output file.
    cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);

    // Set the log message pattern for the writer.
    output->set_pattern("%v");

    
    // Allocate output buffers for each thread.
    output_buffer.resize(thread_count);
    buffer_size.resize(thread_count);

    // Track the maximum sequence buffer size used.
    size_t max_buf_sz = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    uint32_t seqCount = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* const seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;
        const size_t seq_buf_sz = parser->seq.m;

        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Single-threaded writing.
        // output_off_substring(0, seq, seq_len, 0, seq_len - k, output);


        // Multi-threaded writing.
        size_t task_size = (seq_len - k + 1) / thread_count;
        if(!task_size)
            output_plain_off_substring(0, seq, seq_len, 0, seq_len - k, output);
        else
        {
            std::vector<std::thread> task;
            size_t left_end = 0;
            size_t right_end;

            for(uint16_t task_id = 0; task_id < thread_count; ++task_id)
            {
                right_end = (task_id == thread_count - 1 ? seq_len - k : left_end + task_size - 1);
                task.emplace_back(&CdBG::output_plain_off_substring, this, task_id, seq, seq_len, left_end, right_end, output);
                left_end += task_size;
            }

            for(std::thread& t: task){
                if(t.joinable())
                    t.join();
                else
                {
                    std::cerr << "Early termination of a worker thread encountered during writing of maximal unitigs. Aborting.\n";
                    std::exit(EXIT_FAILURE);
                }
            }
        }
    }

    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";


    // Flush the buffers.
    flush_buffers(output);

    
    // Close the loggers?
    spdlog::drop_all();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::output_maximal_unitigs_gfa()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const std::string& ref_file_path = params.ref_file_path();
    const uint16_t thread_count = params.thread_count();
    const std::string& output_file_path = params.output_file_path();
    const std::string& working_dir_path = params.working_dir_path();

    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* const input = fopen(ref_file_path.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << params.ref_file_path() << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Initialize the parser.
    kseq_t* const parser = kseq_init(fileno(input));


    // Clear the output file and write the GFA header.
    std::ofstream op_stream(output_file_path.c_str(), std::ofstream::out);
    if(!op_stream)
    {
        std::cerr << "Error opening output file " << output_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    write_gfa_header(op_stream);
    op_stream.close();


    // Allocate the output buffers for each thread.
    output_buffer.resize(thread_count);
    buffer_size.resize(thread_count);

    // Allocate entries for the first, the second, and the last unitigs seen by each thread.
    first_unitig.resize(thread_count);
    second_unitig.resize(thread_count);
    last_unitig.resize(thread_count);

    // Set the prefixes of the temporary path output files. This is to avoid possible name
    // conflicts in the file system.
    set_temp_file_prefixes(working_dir_path);


    // Track the maximum sequence buffer size used.
    size_t max_buf_sz = 0;

    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    // uint32_t seq_count = 0;
    seq_count = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* const seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;
        const size_t seq_buf_sz = parser->seq.m;

        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);

        std::chrono::high_resolution_clock::time_point t_s = std::chrono::high_resolution_clock::now();
        std::cout << "Processing sequence " << ++seq_count << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Open an asynchronous logger to write into the output file, and set its log message pattern.
        // Note: `spdlog` appends to the output file by default, so the results for the sequences are accumulated into the same output file.
        cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
        output->set_pattern("%v");

        
        // Reset the first, the second, and the last unitigs seen for each thread.
        std::fill(first_unitig.begin(), first_unitig.end(), Oriented_Unitig());
        std::fill(second_unitig.begin(), second_unitig.end(), Oriented_Unitig());
        std::fill(last_unitig.begin(), last_unitig.end(), Oriented_Unitig());

        // Reset the path output streams for each thread.
        reset_path_streams();
        

        // Single-threaded writing.
        // output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        size_t task_size = (seq_len - k + 1) / thread_count;
        if(!task_size)
            output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, output);
        else
        {
            std::vector<std::thread> task;
            size_t left_end = 0;
            size_t right_end;

            for(uint16_t task_id = 0; task_id < thread_count; ++task_id)
            {
                right_end = (task_id == thread_count - 1 ? seq_len - k : left_end + task_size - 1);
                task.emplace_back(&CdBG::output_gfa_off_substring, this, task_id, seq, seq_len, left_end, right_end, output);
                left_end += task_size;
            }

            for(std::thread& t: task)
                if(t.joinable())
                    t.join();
                else
                {
                    std::cerr << "Early termination of a worker thread encountered during writing of maximal unitigs. Aborting.\n";
                    std::exit(EXIT_FAILURE);
                }


            write_inter_thread_connections(output);

            std::chrono::high_resolution_clock::time_point t_e = std::chrono::high_resolution_clock::now();
            double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_e - t_s).count();
            (void)elapsed_seconds;
            // std::cout << "Time taken to write segments and links = " << elapsed_seconds << " seconds.\n";
        }
        
        // Flush all the buffered content (segments and links), as the GFA path to be appended to
        // the same output sink file is written in a different way than using the `spdlog` logger.
        // Note: If using an async logger, `logger->flush()` posts a message to the queue requesting the flush
        // operation, so the function returns immediately. Hence a forceful eviction is necessary by shutdown.
        output->flush();
        spdlog::shutdown();
        

        // Write the GFA path for this sequence.
        params.output_format() == 1 ? write_gfa_path() : write_gfa_ordered_group();
    }

    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";


    // Flush the buffers.
    cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file_path);
    output->set_pattern("%v");
    flush_buffers(output);

    // Remove the temporary files.
    remove_temp_files();

    
    // Close the loggers?
    spdlog::drop_all();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
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
void CdBG<k>::fill_buffer(const uint16_t thread_id, const uint64_t fill_amount, cuttlefish::logger_t output)
{
    buffer_size[thread_id] += fill_amount;

    if(buffer_size[thread_id] > MAX_BUFF_SIZE)
    {
        // TODO: Avoid the presence of the same content in the memory simultaneously.
        // E.g. skipping this line results in consuming memory:
        // "fixed" HG38: 176 instead of 346 MB;
        // HG38: 17 instead of 93 MB;
        // Gorgor3: 0.5 instead of 40MB.

        write(output, output_buffer[thread_id].str());
        
        output_buffer[thread_id].str("");
        buffer_size[thread_id] = 0;
    }
}


template <uint16_t k>
void CdBG<k>::write(cuttlefish::logger_t output, const std::string& str)
{
    output->info("{}", str);
}


template <uint16_t k>
void CdBG<k>::flush_buffers(cuttlefish::logger_t output)
{
    const uint16_t thread_count = params.thread_count();

    for (uint16_t t_id = 0; t_id < thread_count; ++t_id)
        if(buffer_size[t_id] > 0)
        {
            write(output, output_buffer[t_id].str());
            output_buffer[t_id].str("");
            buffer_size[t_id] = 0;
        }
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
