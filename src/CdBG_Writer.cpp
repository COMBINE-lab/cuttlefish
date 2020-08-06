
#include "CdBG.hpp"
#include "kseq/kseq.h"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <chrono>
#include <thread>
#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


void CdBG::output_maximal_unitigs(const std::string& output_file, const uint16_t thread_count)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* input = fopen(ref_file.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Initialize the parser.
    kseq_t* parser = kseq_init(fileno(input));


    // Clear the output file.
    std::ofstream op_stream(output_file.c_str(), std::ofstream::out);
    if(!op_stream)
    {
        std::cerr << "Error opening output file " << output_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    op_stream.close();


    // Open an asynchronous logger to write into the output file.
    cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", output_file.c_str());

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
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;
        const size_t seq_buf_sz = parser->seq.m;

        max_buf_sz = std::max(max_buf_sz, seq_buf_sz);

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Single-threaded writing.
        // output_off_substring(seq, seq_len, 0, seq_len - k, output);


        // Multi-threaded writing.
        size_t task_size = (seq_len - k + 1) / thread_count;
        if(!task_size)
            output_off_substring(0, seq, seq_len, 0, seq_len - k, output);
        else
        {
            std::vector<std::thread> task;
            size_t left_end = 0;
            size_t right_end;

            for(uint16_t task_id = 0; task_id < thread_count; ++task_id)
            {
                right_end = (task_id == thread_count - 1 ? seq_len - k : left_end + task_size - 1);
                task.emplace_back(&CdBG::output_off_substring, this, static_cast<uint64_t>(task_id), seq, seq_len, left_end, right_end, output);
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
    flush_buffers(thread_count, output);

    
    // Close the loggers?
    spdlog::drop_all();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


void CdBG::output_off_substring(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end, cuttlefish::logger_t output)
{
    size_t kmer_idx = left_end;
    while(kmer_idx <= right_end)
    {
        kmer_idx = search_valid_kmer(seq, kmer_idx, right_end);

        // No valid k-mer remains in the sequence.
        if(kmer_idx > right_end)
            break;

        // Process a maximal valid contiguous subsequence, and advance to the index following it.
        kmer_idx = output_maximal_unitigs(thread_id, seq, seq_len, right_end, kmer_idx, output);
    }
}


size_t CdBG::output_maximal_unitigs(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx, cuttlefish::logger_t output)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Annotated_Kmer curr_kmer(cuttlefish::kmer_t(seq, kmer_idx), kmer_idx, Vertices);

    // The subsequence contains only an isolated k-mer, i.e. there's no valid left or right
    // neighboring k-mer to this k-mer. So it's a maximal unitig by itself.
    if((kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE) &&
        (kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE))
        output_unitig(thread_id, seq, curr_kmer, curr_kmer, output);
    else    // At least one valid neighbor exists, either to the left or to the right, or on both sides.
    {
        // No valid right neighbor exists for the k-mer.
        if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
        {
            // A valid left neighbor exists as it's not an isolated k-mer.
            Annotated_Kmer prev_kmer(cuttlefish::kmer_t(seq, kmer_idx - 1), kmer_idx, Vertices);
            
            if(is_unipath_start(curr_kmer.vertex_class(), curr_kmer.dir(), prev_kmer.vertex_class(), prev_kmer.dir()))
                // A maximal unitig ends at the ending of a maximal valid subsequence.
                output_unitig(thread_id, seq, curr_kmer, curr_kmer, output);

            // The contiguous sequence ends at this k-mer.
            return kmer_idx + k;
        }


        // A valid right neighbor exists for the k-mer.
        Annotated_Kmer next_kmer = curr_kmer;
        next_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);

        bool on_unipath = false;
        Annotated_Kmer unipath_start_kmer;
        Annotated_Kmer prev_kmer;

        // No valid left neighbor exists for the k-mer.
        if(kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
        {
            // A maximal unitig starts at the beginning of a maximal valid subsequence.
            on_unipath = true;
            unipath_start_kmer = curr_kmer;
        }
        // Both left and right valid neighbors exist for this k-mer.
        else
        {
            prev_kmer = Annotated_Kmer(cuttlefish::kmer_t(seq, kmer_idx - 1), kmer_idx, Vertices);
            if(is_unipath_start(curr_kmer.vertex_class(), curr_kmer.dir(), prev_kmer.vertex_class(), prev_kmer.dir()))
            {
                on_unipath = true;
                unipath_start_kmer = curr_kmer;
            }
        }

        if(on_unipath && is_unipath_end(curr_kmer.vertex_class(), curr_kmer.dir(), next_kmer.vertex_class(), next_kmer.dir()))
        {
            output_unitig(thread_id, seq, unipath_start_kmer, curr_kmer, output);
            on_unipath = false;
        }


        // Process the rest of the k-mers of this contiguous subsequence.
        for(kmer_idx++; on_unipath || kmer_idx <= right_end; ++kmer_idx)
        {
            prev_kmer = curr_kmer;
            curr_kmer = next_kmer;

            if(is_unipath_start(curr_kmer.vertex_class(), curr_kmer.dir(), prev_kmer.vertex_class(), prev_kmer.dir()))
            {
                on_unipath = true;
                unipath_start_kmer = curr_kmer;
            }


            // No valid right neighbor exists for the k-mer.
            if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
            {
                // A maximal unitig ends at the ending of a maximal valid subsequence.
                if(on_unipath)
                {
                    output_unitig(thread_id, seq, unipath_start_kmer, curr_kmer, output);
                    on_unipath = false;
                }

                // The contiguous sequence ends at this k-mer.
                return kmer_idx + k;
            }
            else    // A valid right neighbor exists.
            {
                next_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);
                
                if(on_unipath && is_unipath_end(curr_kmer.vertex_class(), curr_kmer.dir(), next_kmer.vertex_class(), next_kmer.dir()))
                {
                    output_unitig(thread_id, seq, unipath_start_kmer, curr_kmer, output);
                    on_unipath = false;
                }
            }
        }
    }
    
    
    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


bool CdBG::is_unipath_start(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::dir_t dir, const cuttlefish::Vertex_Class prev_kmer_class, const cuttlefish::dir_t prev_kmer_dir) const
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


bool CdBG::is_unipath_end(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::dir_t dir, const cuttlefish::Vertex_Class next_kmer_class, const cuttlefish::dir_t next_kmer_dir) const
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


void CdBG::output_unitig(const uint64_t thread_id, const char* const seq, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, cuttlefish::logger_t output)
{
    // This is to avoid race conditions that may arise while multi-threading.
    // If two threads try to output the same unitig at the same time but
    // encounter it in the opposite orientations, then data races may arise.
    // For a particular unitig, always query the same well-defined canonical flanking
    // k-mer, irrespective of which direction the unitig may be traversed at.
    const cuttlefish::kmer_t min_flanking_kmer = std::min(start_kmer.canonical(), end_kmer.canonical());
    Kmer_Hash_Entry_API hash_table_entry = Vertices[min_flanking_kmer];
    State& state = hash_table_entry.get_state();

    if(state.is_outputted())
        return;
    

    state = state.outputted();

    // If the hash table update is successful, only then this thread may output this unitig.
    if(Vertices.update(hash_table_entry))
        write_path(thread_id, seq, start_kmer.idx(), end_kmer.idx(), start_kmer.kmer() < end_kmer.rev_compl(), output);
}


void CdBG::write_path(const uint64_t thread_id, const char* const seq, const size_t start_kmer_idx, const size_t end_kmer_idx, const cuttlefish::dir_t dir, cuttlefish::logger_t output) 
{
    std::stringstream& buffer = output_buffer[thread_id];
    const size_t path_len = end_kmer_idx - start_kmer_idx + k;

    if(dir == cuttlefish::FWD)
        for(size_t offset = 0; offset < path_len; ++offset)
            buffer << seq[start_kmer_idx + offset];
    else    // dir == cuttlefish::BWD
        for(size_t offset = 0; offset < path_len; ++offset)
            buffer << cuttlefish::kmer_t::complement(seq[end_kmer_idx + k - 1 - offset]);

    // End the path.
    buffer << "\n";

    
    // TODO: Fix a max memory scheme for buffers instead of a fixed line count.
    // Mark buffer size increment.
    fill_buffer(thread_id, 1, output);
}


void CdBG::fill_buffer(const uint64_t thread_id, const uint64_t fill_amount, cuttlefish::logger_t output)
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


void CdBG::write(cuttlefish::logger_t output, const std::string& str)
{
    output->info("{}", str);
}


void CdBG::flush_buffers(const uint16_t thread_count, cuttlefish::logger_t output)
{
    for (uint16_t t_id = 0; t_id < thread_count; ++t_id)
        if(buffer_size[t_id] > 0)
        {
            write(output, output_buffer[t_id].str());
            output_buffer[t_id].str("");
            buffer_size[t_id] = 0;
        }
}
