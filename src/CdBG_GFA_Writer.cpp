
#include "CdBG.hpp"
#include "kseq/kseq.h"
#include "utility.hpp"

#include <chrono>
#include <thread>
#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


void CdBG::output_maximal_unitigs_gfa(const std::string& gfa_file_name, const uint8_t gfa_v, const uint16_t thread_count, const std::string& working_dir)
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


    // Clear the output file and write the GFA header.
    std::ofstream op_stream(gfa_file_name.c_str(), std::ofstream::out);
    if(!op_stream)
    {
        std::cerr << "Error opening output file " << gfa_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    write_gfa_header(gfa_v, op_stream);
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
    set_temp_file_prefixes(working_dir);


    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    // uint32_t seq_count = 0;
    seq_count = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::chrono::high_resolution_clock::time_point t_s = std::chrono::high_resolution_clock::now();
        std::cout << "Processing sequence " << ++seq_count << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            continue;


        // Open an asynchronous logger to write into the output file, and set its log message pattern.
        // Note: `spdlog` appends to the output file by default, so the results for the sequences are accumulated into the same output file.
        cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", gfa_file_name.c_str());
        output->set_pattern("%v");

        
        // Reset the first, the second, and the last unitigs seen for each thread.
        std::fill(first_unitig.begin(), first_unitig.end(), Oriented_Unitig());
        std::fill(second_unitig.begin(), second_unitig.end(), Oriented_Unitig());
        std::fill(last_unitig.begin(), last_unitig.end(), Oriented_Unitig());

        // Reset the path output streams for each thread.
        reset_path_streams(gfa_v, thread_count);
        

        // Single-threaded writing.
        // output_off_substring(seq, seq_len, 0, seq_len - k, output);

        // Multi-threaded writing.
        size_t task_size = (seq_len - k + 1) / thread_count;
        if(!task_size)
            output_gfa_off_substring(0, seq, seq_len, 0, seq_len - k, gfa_v, output);
        else
        {
            std::vector<std::thread> task;
            size_t left_end = 0;
            size_t right_end;

            for(uint16_t task_id = 0; task_id < thread_count; ++task_id)
            {
                right_end = (task_id == thread_count - 1 ? seq_len - k : left_end + task_size - 1);
                task.emplace_back(&CdBG::output_gfa_off_substring, this, static_cast<uint64_t>(task_id), seq, seq_len, left_end, right_end, gfa_v, output);
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


            write_inter_thread_connections(thread_count, gfa_v, output);

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
        gfa_v == 1 ? write_gfa_path(thread_count, gfa_file_name) : write_gfa_ordered_group(thread_count, gfa_file_name);
    }


    // Flush the buffers.
    cuttlefish::logger_t output = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", gfa_file_name.c_str());
    output->set_pattern("%v");
    flush_buffers(thread_count, output);

    // Remove the temporary files.
    remove_temp_files(thread_count, gfa_v);

    
    // Close the loggers?
    spdlog::drop_all();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


void CdBG::set_temp_file_prefixes(const std::string& working_dir)
{
    // Check if temporary files can be created with random names.
    const uint64_t RETRY_COUNT = 10;
    std::string temp_file_prefix;
    
    for(uint64_t attempt = 0; attempt < RETRY_COUNT; ++attempt)
    {
        temp_file_prefix = get_random_string(TEMP_FILE_PREFIX_LEN);
        if(!file_prefix_exists(working_dir, temp_file_prefix))
        {
            PATH_OUTPUT_PREFIX = working_dir + "/" + PATH_OUTPUT_PREFIX + temp_file_prefix;
            OVERLAP_OUTPUT_PREFIX = working_dir + "/" + OVERLAP_OUTPUT_PREFIX + temp_file_prefix;

            std::cout << "Temporary file name prefixes: " << PATH_OUTPUT_PREFIX << "\n";

            return;
        }
    }

    std::cerr << "Failed to find any random prefix for temporary file names. Aborting.\n";
    std::exit(EXIT_FAILURE);
}


void CdBG::reset_path_streams(const uint8_t gfa_v, const uint16_t thread_count)
{
    path_output.clear();
    if(gfa_v == 1)
        overlap_output.clear();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        path_output.emplace_back((PATH_OUTPUT_PREFIX + std::to_string(t_id)).c_str(), std::ofstream::out);

        if(gfa_v == 1)
            overlap_output.emplace_back((OVERLAP_OUTPUT_PREFIX + std::to_string(t_id)).c_str(), std::ofstream::out);
    }
}


void CdBG::output_gfa_off_substring(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end, const uint8_t gfa_v, cuttlefish::logger_t output)
{
    size_t kmer_idx = left_end;
    while(kmer_idx <= right_end)
    {
        kmer_idx = search_valid_kmer(seq, kmer_idx, right_end);

        // No valid k-mer remains in the sequence.
        if(kmer_idx > right_end)
            break;

        // Process a maximal valid contiguous subsequence, and advance to the index following it.
        kmer_idx = output_maximal_unitigs_gfa(thread_id, seq, seq_len, right_end, kmer_idx, gfa_v, output);
    }
}


size_t CdBG::output_maximal_unitigs_gfa(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx, const uint8_t gfa_v, cuttlefish::logger_t output)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Annotated_Kmer curr_kmer(cuttlefish::kmer_t(seq, kmer_idx), kmer_idx, Vertices);

    // The subsequence contains only an isolated k-mer, i.e. there's no valid left or right
    // neighboring k-mer to this k-mer. So it's a maximal unitig by itself.
    if((kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE) &&
        (kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE))
        output_unitig_gfa(thread_id, seq, curr_kmer, curr_kmer, gfa_v, output);
    else    // At least one valid neighbor exists, either to the left or to the right, or on both sides.
    {
        // No valid right neighbor exists for the k-mer.
        if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
        {
            // A valid left neighbor exists as it's not an isolated k-mer.
            Annotated_Kmer prev_kmer(cuttlefish::kmer_t(seq, kmer_idx - 1), kmer_idx, Vertices);
            
            if(is_unipath_start(curr_kmer.vertex_class(), curr_kmer.dir(), prev_kmer.vertex_class(), prev_kmer.dir()))
                // A maximal unitig ends at the ending of a maximal valid subsequence.
                output_unitig_gfa(thread_id, seq, curr_kmer, curr_kmer, gfa_v, output);

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
            output_unitig_gfa(thread_id, seq, unipath_start_kmer, curr_kmer, gfa_v, output);
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
                    output_unitig_gfa(thread_id, seq, unipath_start_kmer, curr_kmer, gfa_v, output);
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
                    output_unitig_gfa(thread_id, seq, unipath_start_kmer, curr_kmer, gfa_v, output);
                    on_unipath = false;
                }
            }
        }
    }
    
    
    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


void CdBG::output_unitig_gfa(const uint64_t thread_id, const char* const seq, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, const uint8_t gfa_v, cuttlefish::logger_t output)
{
    // This is to avoid race conditions that may arise while multi-threading.
    // If two threads try to output the same unitig at the same time but
    // encounter it in the opposite orientations, then data races may arise.
    // For a particular unitig, always query the same well-defined canonical flanking
    // k-mer, irrespective of which direction the unitig may be traversed at.
    const cuttlefish::kmer_t min_flanking_kmer = std::min(start_kmer.canonical(), end_kmer.canonical());
    const uint64_t bucket_id = Vertices.bucket_id(min_flanking_kmer);
    Kmer_Hash_Entry_API hash_table_entry = Vertices[bucket_id];
    State& state = hash_table_entry.get_state();

    // Name the GFA segment with the hash value of the first k-mer of the canonical form unitig.
    const uint64_t unitig_id = bucket_id;
    const cuttlefish::dir_t unitig_dir = (start_kmer.kmer() < end_kmer.rev_compl() ? cuttlefish::FWD : cuttlefish::BWD);
    const Oriented_Unitig current_unitig(unitig_id, unitig_dir, start_kmer.idx(), end_kmer.idx());


    // Output a possible GFA segment.

    if(!state.is_outputted())
    {
        state = state.outputted();

        // If the hash table update is successful, only then this thread may output this unitig.
        if(Vertices.update(hash_table_entry))
            write_gfa_segment(thread_id, seq, unitig_id, start_kmer.idx(), end_kmer.idx(), unitig_dir, gfa_v, output);
    }


    // Output a possible GFA connection (link, edge, or gap).

    if(!first_unitig[thread_id].is_valid())
        first_unitig[thread_id] = current_unitig;
    else if(!second_unitig[thread_id].is_valid())
        second_unitig[thread_id] = current_unitig;

    Oriented_Unitig& prev_unitig = last_unitig[thread_id];
    if(prev_unitig.is_valid())
        write_gfa_connection(thread_id, prev_unitig, current_unitig, gfa_v, output);
    
    prev_unitig = current_unitig;
}


void CdBG::write_gfa_header(const uint8_t gfa_v, std::ofstream& output) const
{
    // The GFA header record.
    if(gfa_v == 1)
        output << GFA1_HEADER;
    else if(gfa_v == 2)
        output << GFA2_HEADER;
    else
    {
        std::cerr << "Unsupported / nonexistent GFA version number provided. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // End the header line.
    output << "\n";
}


void CdBG::write_gfa_segment(const uint64_t thread_id, const char* const seq, const uint64_t segment_name, const size_t start_kmer_idx, const size_t end_kmer_idx, const cuttlefish::dir_t dir, const uint8_t gfa_v, cuttlefish::logger_t output)
{
    std::stringstream& buffer = output_buffer[thread_id];
    const size_t segment_len = end_kmer_idx - start_kmer_idx + k;

    // The 'RecordType' field for segment lines.
    buffer << "S";
    
    // The 'Name' field.
    buffer << "\t" << segment_name;

    // The 'SegmentLength' field (required for GFA2).
    if(gfa_v == 2)
        buffer << "\t" << segment_len;
    
    // The segment field.
    buffer << "\t";
    if(dir == cuttlefish::FWD)
        for(size_t offset = 0; offset < segment_len; ++offset)
            buffer << seq[start_kmer_idx + offset];
    else
        for(size_t offset = 0; offset < segment_len; ++offset)
            buffer << complement(seq[end_kmer_idx + k - 1 - offset]);


    // Write some optional fields that are trivially inferrable here.
    if(gfa_v == 1)  // No need of the length tag for GFA2 here.
        buffer << "\tLN:i:" << segment_len; // The segment length.

    buffer << "\tKC:i:" << (end_kmer_idx - start_kmer_idx + 1); // The k-mer count.


    // End the segment line.
    buffer << "\n";


    // Mark buffer size increment.
    fill_buffer(thread_id, 1, output);
}


void CdBG::write_gfa_connection(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, const uint8_t gfa_v, cuttlefish::logger_t output)
{
    if(gfa_v == 1)
        write_gfa_link(thread_id, left_unitig, right_unitig, output);
    else    // GFA2
        if(right_unitig.start_kmer_idx == left_unitig.end_kmer_idx + 1)
            write_gfa_edge(thread_id, left_unitig, right_unitig, output);
        else
            write_gfa_gap(thread_id, left_unitig, right_unitig, output);
}


void CdBG::write_gfa_link(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output)
{
    std::stringstream& buffer = output_buffer[thread_id];

    // The 'RecordType' field for link lines.
    buffer << "L";

    // The 'From' fields.
    buffer << "\t" << left_unitig.unitig_id << "\t" << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // The 'To' fields.
    buffer << "\t" << right_unitig.unitig_id << "\t" << (right_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // The 'Overlap' field.
    const uint16_t overlap = (right_unitig.start_kmer_idx == left_unitig.end_kmer_idx + 1 ? k - 1 : 0);
    buffer << "\t" << overlap << "M";

    // End the link line.
    buffer << "\n";


    // Mark buffer size increment.
    fill_buffer(thread_id, 1, output);

    
    // Append a link to the growing path for this thread.
    append_link_to_path(thread_id, left_unitig, right_unitig);
}


void CdBG::write_gfa_edge(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output)
{
    std::stringstream& buffer = output_buffer[thread_id];

    // The 'RecordType' field for edge lines.
    buffer << "E";

    // The 'Edge-ID' field.
    buffer << "\t*";

    // The 'Segment-ID' fields.
    buffer << "\t" << left_unitig.unitig_id << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");
    buffer << "\t" << right_unitig.unitig_id << (right_unitig.dir == cuttlefish::FWD ? "+" : "-");


    // The 'Begin' and 'End' fields for the first segment.
    size_t unitig_len = left_unitig.length();
    if(left_unitig.dir == cuttlefish::FWD)
        buffer << "\t" << (unitig_len - (k - 1)) << "\t" << unitig_len << "$";
    else
        buffer << "\t" << 0 << "\t" << (k - 1);

    // The 'Begin' and 'End' fields for the second segment.
    unitig_len = right_unitig.length();
    if(right_unitig.dir == cuttlefish::FWD)
        buffer << "\t" << 0 << "\t" << (k - 1);
    else
        buffer << "\t" << (unitig_len - (k - 1)) << "\t" << unitig_len << "$";

    
    // The 'Alignment' field.
    buffer << "\t*";

    // End the edge line.
    buffer << "\n";


    // Mark buffer size increment.
    fill_buffer(thread_id, 1, output);

    
    // Append an edge to the growing path for this thread.
    append_edge_to_path(thread_id, left_unitig, right_unitig);
}


void CdBG::write_gfa_gap(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output)
{
    std::stringstream& buffer = output_buffer[thread_id];

    // The 'RecordType' field for gap lines.
    buffer << "G";

    // The 'Gap-ID' field.
    buffer << "\t*";

    // The 'Segment-ID' fields.
    buffer << "\t" << left_unitig.unitig_id << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");
    buffer << "\t" << right_unitig.unitig_id << (right_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // Write the 'Distance' field.
    buffer << "\t" << (right_unitig.start_kmer_idx - (left_unitig.end_kmer_idx + k));

    // Write the `Variance` field.
    buffer << "\t*";

    // End the gap line.
    buffer << "\n";


    // Mark buffer size increment.
    fill_buffer(thread_id, 1, output);


    // Append an edge to the growing path for this thread.
    append_edge_to_path(thread_id, left_unitig, right_unitig);
}


void CdBG::append_link_to_path(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig)
{
    // The destination vertex (unitig) is written for each link.
    // Note that, the very first vertex of the path tiling for the sequence is thus missing in the path outputs.

    path_output[thread_id] << "," << right_unitig.unitig_id << (right_unitig.dir == cuttlefish::FWD ? "+" : "-");
    overlap_output[thread_id] << "," << (right_unitig.start_kmer_idx == left_unitig.end_kmer_idx + 1 ? k - 1 : 0) << "M";

    // Error checking for writing failures is done while closing the streams.
}


void CdBG::append_edge_to_path(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig)
{
    // The destination vertex (unitig) is written for each edge.
    // Note that, the very first vertex of the path tiling for the sequence is thus missing in the path outputs.

    (void)left_unitig;
    path_output[thread_id] << " " << right_unitig.unitig_id << (right_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // Error checking for writing failures is done while closing the streams.
}


void CdBG::write_inter_thread_connections(const uint16_t thread_count, const uint8_t gfa_v, cuttlefish::logger_t output)
{
    Oriented_Unitig left_unitig;
    uint16_t left_t_id;

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        if(!left_unitig.is_valid())
        {
            left_unitig = last_unitig[t_id];
            left_t_id = t_id;
        }
        else
            if(first_unitig[t_id].is_valid())
            {
                // A link exists between the last unitig of the thread number `left_t_id`
                // and the first unitig of the thread number `t_id`.
                const Oriented_Unitig& right_unitig = first_unitig[t_id];

                write_gfa_connection(left_t_id, left_unitig, right_unitig, gfa_v, output);

                // There definitely exists a last unitig for the thread number `t_id`, as it has a first unitig.
                left_unitig = last_unitig[t_id];
                left_t_id = t_id;
            }
}


void CdBG::search_first_connection(const uint16_t thread_count, Oriented_Unitig& left_unitig, Oriented_Unitig& right_unitig) const
{
    left_unitig = Oriented_Unitig();
    right_unitig = Oriented_Unitig();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        if(first_unitig[t_id].is_valid())
        {
            if(!left_unitig.is_valid())
                left_unitig = first_unitig[t_id];
            else
            {
                right_unitig = first_unitig[t_id];
                return;
            }
        }

        if(second_unitig[t_id].is_valid())
        {
            // Obviously, `first_unitig[t_id]` must also be valid for the thread number `t_id`;
            // so `left_unitig` is already set to a valid value at this point.
            right_unitig = second_unitig[t_id];
            return;
        }
    }
}


void CdBG::write_gfa_path(const uint16_t thread_count, const std::string& gfa_file_name)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    // Close the path output streams.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        if(path_output[t_id].fail() || overlap_output[t_id].fail())
        {
            std::cerr << "Errors encountered while writing to the temporary path output files. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        path_output[t_id].close();
        overlap_output[t_id].close();
    }
    

    // Search the very first GFA link in the sequence, as that is not inferrable from the path outputs.
    Oriented_Unitig left_unitig;
    Oriented_Unitig right_unitig;
    search_first_connection(thread_count, left_unitig, right_unitig);

    // The sequence does not contain any unitig (possible if there's no valid k-mer in the sequence).
    if(!left_unitig.is_valid())
        return;

    
    // Open the output file in append mode.
    std::ofstream output(gfa_file_name.c_str(), std::ios_base::app);

    // The 'RecordType' field for the path lines.
    output << "P";

    // The 'PathName' field.
    output << "\tP" << seq_count;

    // The 'SegmentNames' field.
    output << "\t";
    
    // The first vertex of the path (not inferrable from the path output files).
    output << left_unitig.unitig_id << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // Copy the thread-specific path output file contents to the GFA output file.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const std::string path_file_name = (PATH_OUTPUT_PREFIX + std::to_string(t_id));
        std::ifstream input(path_file_name.c_str(), std::ios_base::in);

        if(input.fail())
        {
            std::cerr << "Error opening temporary path output file " << path_file_name << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
        
        // Copy the path output for thread number `t_id` to the end of the output GFA file.
        if(input.peek() != EOF)
            output << input.rdbuf();

        input.close();
    }


    // The 'Overlaps' field.
    output << "\t";

    // The sequence contains only one unitig.
    if(!right_unitig.is_valid())
        output << "*";  // Write an empty CIGAR string at the 'Overlaps' field.
    else
    {
        // The first overlap of the path (not inferrable from the path output files).
        const uint16_t overlap = (right_unitig.start_kmer_idx == left_unitig.end_kmer_idx +  1 ? k - 1 : 0);
        output << overlap << "M";

        // Copy the thread-specific overlap output file contents to the GFA output file.
        for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        {
            const std::string overlap_file_name = (OVERLAP_OUTPUT_PREFIX + std::to_string(t_id));
            std::ifstream input(overlap_file_name.c_str(), std::ifstream::in);

            if(input.fail())
            {
                std::cerr << "Error opening temporary path output file " << overlap_file_name << ". Aborting.\n";
                std::exit(EXIT_FAILURE);
            }

            // Copy the overlaps output for thread number `t_id` to the end of the output GFA file.
            if(input.peek() != EOF)
                output << input.rdbuf();

            input.close();
        }
    }


    // End the path line.
    output << "\n";

    output.close();

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    (void)elapsed_seconds;
    // std::cout << "Time taken to write paths = " << elapsed_seconds << " seconds.\n";
}


void CdBG::write_gfa_ordered_group(const uint16_t thread_count, const std::string& gfa_file_name)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    // Close the path output streams.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        if(path_output[t_id].fail())
        {
            std::cerr << "Errors encountered while writing to the temporary path output files. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        path_output[t_id].close();
    }
    

    // Search the very first GFA edge in the sequence, as that is not inferrable from the path outputs.
    Oriented_Unitig left_unitig;
    Oriented_Unitig right_unitig;
    search_first_connection(thread_count, left_unitig, right_unitig);

    // The sequence does not contain any unitig (possible if there's no valid k-mer in the sequence).
    if(!left_unitig.is_valid())
        return;

    
    // Open the output file in append mode.
    std::ofstream output(gfa_file_name.c_str(), std::ios_base::app);

    // The 'RecordType' field for the ordered group line.
    output << "O";

    // The 'Group-ID' field.
    output << "\tP" << seq_count;

    // The 'Members' field.
    output << "\t";
    
    // The first vertex of the path (not inferrable from the path output files).
    output << left_unitig.unitig_id << (left_unitig.dir == cuttlefish::FWD ? "+" : "-");

    // Copy the thread-specific path output file contents to the GFA output file.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const std::string path_file_name = (PATH_OUTPUT_PREFIX + std::to_string(t_id));
        std::ifstream input(path_file_name.c_str(), std::ios_base::in);

        if(input.fail())
        {
            std::cerr << "Error opening temporary path output file " << path_file_name << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
        
        // Copy the path output for thread number `t_id` to the end of the output GFA file.
        if(input.peek() != EOF)
            output << input.rdbuf();

        input.close();
    }


    // End the path line.
    output << "\n";

    output.close();

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    (void)elapsed_seconds;
    // std::cout << "Time taken to write paths = " << elapsed_seconds << " seconds.\n";
}


void CdBG::remove_temp_files(const uint16_t thread_count, const uint8_t gfa_v) const
{
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const std::string path_file_name = PATH_OUTPUT_PREFIX + std::to_string(t_id);
        const std::string overlap_file_name = OVERLAP_OUTPUT_PREFIX + std::to_string(t_id);

        if(remove(path_file_name.c_str()) != 0 || (gfa_v == 1 && remove(overlap_file_name.c_str()) != 0))
        {
            std::cerr << "Error deleting temporary files. Aborting\n";
            std::exit(EXIT_FAILURE);
        }
    }
}
