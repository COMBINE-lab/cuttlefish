
#include "CdBG_Builder.hpp"
#include "kseq/kseq.h"
#include "Kmer_Iterator.hpp"

#include <fstream>
#include <thread>
#include <cstdlib>
#include <cassert>
#include <chrono>
#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


CdBG_Builder::CdBG_Builder(const std::string& ref_file, const uint16_t k):
    ref_file(ref_file), k(k)
{
    Kmer::set_k(k);
}


void CdBG_Builder::construct(const std::string& kmc_file_name, const std::string& bbhash_file_name, const uint16_t thread_count, const std::string& output_file_name)
{
    std::cout << "Constructing the minimal perfect hash function.\n";
    Vertices.construct(kmc_file_name, bbhash_file_name, thread_count);

    std::cout << "Classifying the vertices.\n";
    classify_vertices(thread_count);

    std::cout << "Outputting the maximal unitigs.\n";
    output_maximal_unitigs(output_file_name);

    Vertices.clear();
}


void CdBG_Builder::classify_vertices(const uint16_t thread_count)
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

    // Parse sequences one-by-one, and continue partial classification of the k-mers through them.
    uint32_t seqCount = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            break;


        // Single-threaded classification.
        // process_substring(seq, seq_len, 0, seq_len - k);


        // Multi-threaded classification.
        size_t task_size = (seq_len - k + 1) / thread_count;
        if(!task_size)
            process_substring(seq, seq_len, 0, seq_len - k);
        else
        {
            std::vector<std::thread> task;
            size_t left_end = 0;
            size_t right_end;

            for(uint16_t task_id = 0; task_id < thread_count; ++task_id)
            {
                right_end = (task_id == thread_count - 1 ? seq_len - k : left_end + task_size - 1);
                task.emplace_back(&CdBG_Builder::process_substring, this, seq, seq_len, left_end, right_end);
                left_end += task_size;
            }

            for(std::thread& t: task)
                if(t.joinable())
                    t.join();
                else
                {
                    std::cerr << "Early termination of a worker thread encountered. Aborting.\n";
                    std::exit(EXIT_FAILURE);
                }
        } 
    }


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done classifying the vertices. Time taken = " << elapsed_seconds << " seconds.\n";
}


void CdBG_Builder::process_substring(const char* seq, const size_t seq_len, const size_t left_end, const size_t right_end)
{
    size_t kmer_idx = left_end;
    while(kmer_idx <= right_end)
    {
        kmer_idx = search_valid_kmer(seq, kmer_idx, right_end);

        // No valid k-mer remains in the substring anymore.
        if(kmer_idx > right_end)
            break;

        // Process a maximal valid contiguous subsequence, and advance to the index following it.
        kmer_idx = process_contiguous_subseq(seq, seq_len, right_end, kmer_idx);
    }
}


size_t CdBG_Builder::search_valid_kmer(const char* seq, const size_t left_end, const size_t right_end)
{
    size_t valid_start_idx;
    uint16_t nucl_count;
    

    size_t idx = left_end;
    while(idx <= right_end)
    {
        // Go over the contiguous subsequence of 'N's.
        for(; idx <= right_end && seq[idx] == 'N'; idx++);

        // Go over the contiguous subsequence of non-'N's.
        if(idx <= right_end)
        {
            valid_start_idx = idx;
            nucl_count = 0;

            for(; idx <= right_end + k - 1 && seq[idx] != 'N'; ++idx)
                if(++nucl_count == k)
                    return valid_start_idx;
        }
    }


    return right_end + 1;
}


size_t CdBG_Builder::process_contiguous_subseq(const char* seq, const size_t seq_len, const size_t right_end, const size_t start_idx)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Directed_Kmer curr_kmer(cuttlefish::kmer_t(seq, kmer_idx));

    // The subsequence contains only an isolated k-mer,
    // i.e. there's no valid left or right neighboring k-mer to this k-mer.
    if((kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE) &&
        (kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE))
        while(!process_isolated_kmer(curr_kmer.canonical));
    else    // At least one valid neighbor exists, either to the left or to the right, or on both sides.
    {
        // Process the leftmost k-mer of this contiguous subsequence.

        // No valid right neighbor exists for the k-mer.
        if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
        {
            // A valid left neighbor exists at it's not an isolated k-mer.
            while(!process_rightmost_kmer(curr_kmer.canonical, curr_kmer.dir, seq[kmer_idx - 1]));

            // The contiguous sequence ends at this k-mer.
            return kmer_idx + k;
        }

        // A valid right neighbor exists for the k-mer.
        Directed_Kmer next_kmer = curr_kmer;
        next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);
        
        // No valid left neighbor exists for the k-mer.
        if(kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
            while(!process_leftmost_kmer(curr_kmer.canonical, curr_kmer.dir, next_kmer.canonical, seq[kmer_idx + k]));
        // Both left and right valid neighbors exist for this k-mer.
        else
            while(!process_internal_kmer(curr_kmer.canonical, curr_kmer.dir, next_kmer.canonical, seq[kmer_idx - 1], seq[kmer_idx + k]));
        

        // Process the internal k-mers of this contiguous subsequence.
        // Each of these k-mers have valid neighbors to their left and right.
        for(kmer_idx++; kmer_idx < right_end && seq[kmer_idx + k] != cuttlefish::PLACEHOLDER_NUCLEOTIDE; ++kmer_idx)
        {
            curr_kmer = next_kmer;
            next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

            while(!process_internal_kmer(curr_kmer.canonical, curr_kmer.dir, next_kmer.canonical, seq[kmer_idx - 1], seq[kmer_idx + k]));
        }


        // Process the rightmost k-mer of this contiguous subsequence. This does not coincide with the leftmost k-mer.

        if(kmer_idx <= right_end)   // Required for the cases where the provided range has just length 1, so `start_idx = right_end`.
        {
            curr_kmer = next_kmer;
        
            // No valid right neighbor exists for the k-mer.
            if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
                while(!process_rightmost_kmer(curr_kmer.canonical, curr_kmer.dir, seq[kmer_idx - 1]));
            // A valid right neighbor exists for the k-mer.
            else
            {
                next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

                while(!process_internal_kmer(curr_kmer.canonical, curr_kmer.dir, next_kmer.canonical, seq[kmer_idx - 1], seq[kmer_idx + k]));
            }
        }
        else
            kmer_idx--; // `kmer_idx` has to be the index of the last valid k-mer encountered.
    }


    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


bool CdBG_Builder::is_self_loop(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_t& next_kmer_hat) const
{
    return kmer_hat == next_kmer_hat;
}


bool CdBG_Builder::process_leftmost_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t next_nucl)
{
    // Fetch the entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    Vertex_Encoding& vertex_encoding = hash_table_entry.get_vertex_encoding();
    // const Vertex_Encoding old_encoding = vertex_encoding;


    // The k-mer is already classified as a complex node.
    if(vertex_encoding.is_visited() && vertex_encoding.vertex_class() == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
        vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::multi_in_multi_out));
    else if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::multi_in_single_out, next_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.vertex_class == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.exit == next_nucl)
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_single_out;
                else
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit != next_nucl)
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::single_in_multi_out, complement(next_nucl)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.vertex_class == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter == complement(next_nucl))
                    vertex.vertex_class = cuttlefish::Vertex_Class::single_in_multi_out;
                else
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter != complement(next_nucl))
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }


    // Ideally, we can skip calling the `update` method when `vertex_encoding` has not changed here,
    // i.e. `vertex_encoding` remains equal to the value with which it had been initialized. But in a
    // multi-threaded environment, the actual location of `vertex_encoding` in the hash table could
    // have been changed in between the initialization of the variable at the beginning and this exit
    // point of the method. For such cases, the `update` would fail, signalling the caller of the event.
    // TODO: We may also get away without updating the same value again, as the edge-ordering should not
    // matter in the algorithm.
    return Vertices.update(hash_table_entry);
}


bool CdBG_Builder::process_rightmost_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::nucleotide_t prev_nucl)
{
    // Fetch the entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    Vertex_Encoding& vertex_encoding = hash_table_entry.get_vertex_encoding();
    // const Vertex_Encoding old_encoding = vertex_encoding;


    // The k-mer is already classified as a complex node.
    if(vertex_encoding.is_visited() && vertex_encoding.vertex_class() == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::single_in_multi_out, prev_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.vertex_class == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter == prev_nucl)
                    vertex.vertex_class = cuttlefish::Vertex_Class::single_in_multi_out;
                else
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;
                
                vertex_encoding = Vertex_Encoding(vertex);
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter != prev_nucl)
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::multi_in_single_out, complement(prev_nucl)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.vertex_class == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.exit == complement(prev_nucl))
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_single_out;
                else
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding= Vertex_Encoding(vertex);
            }
            else if(vertex.vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
        }
    }


    // Ideally, we can skip calling the `update` method when `vertex_encoding` has not changed here,
    // i.e. `vertex_encoding` remains equal to the value with which it had been initialized. But in a
    // multi-threaded environment, the actual location of `vertex_encoding` in the hash table could
    // have been changed in between the initialization of the variable at the beginning and this exit
    // point of the method. For such cases, the `update` would fail, signalling the caller of the event.
    // TODO: We may also get away without updating the same value again, as the edge-ordering should not
    // matter in the algorithm.
    return Vertices.update(hash_table_entry);
}


bool CdBG_Builder::process_internal_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t prev_nucl, const cuttlefish::nucleotide_t next_nucl)
{
    // Fetch the hash table entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    Vertex_Encoding& vertex_encoding = hash_table_entry.get_vertex_encoding();
    // const Vertex_Encoding old_encoding = vertex_encoding;


    // The k-mer is already classified as a complex node.
    if(vertex_encoding.is_visited() && vertex_encoding.vertex_class() == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.


    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
        vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::multi_in_multi_out));
    else if(dir == cuttlefish::FWD)
    {
        // The k-mer is encountered for the first time, and in the forward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::single_in_single_out, prev_nucl, next_nucl));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.vertex_class == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter == prev_nucl && vertex.exit == next_nucl)
                    ; // return true; // Probable race if returned early from here
                    // TODO: Maybe not. May get away without updating the same value again,
                    // as the edge-ordering should not matter in the algorithm.
                else if(vertex.enter != prev_nucl && vertex.exit != next_nucl)
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;
                else if(vertex.enter != prev_nucl)
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_single_out;
                else    // vertex.exit != next_nucl
                    vertex.vertex_class = cuttlefish::Vertex_Class::single_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit != next_nucl)
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter != prev_nucl)
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }
    else
    {
        // The k-mer is encountered for the first time, and in the backward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::single_in_single_out, complement(next_nucl), complement(prev_nucl)));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.vertex_class == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter == complement(next_nucl) && vertex.exit == complement(prev_nucl))
                    ; // return true; // Probable race if returned early from here
                    // TODO: Maybe not. May get away without updating the same value again,
                    // as the edge-ordering should not matter in the algorithm.
                else if(vertex.enter != complement(next_nucl) && vertex.exit != complement(prev_nucl))
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;
                else if(vertex.enter != complement(next_nucl))
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_single_out;
                else    // vertex.exit != complement(prev_nucl)
                    vertex.vertex_class = cuttlefish::Vertex_Class::single_in_multi_out;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter != complement(next_nucl))
                {
                    vertex.vertex_class = cuttlefish::Vertex_Class::multi_in_multi_out;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }


    // Ideally, we can skip calling the `update` method when `vertex_encoding` has not changed here,
    // i.e. `vertex_encoding` remains equal to the value with which it had been initialized. But in a
    // multi-threaded environment, the actual location of `vertex_encoding` in the hash table could
    // have been changed in between the initialization of the variable at the beginning and this exit
    // point of the method. For such cases, the `update` would fail, signalling the caller of the event.
    // TODO: We may also get away without updating the same value again, as the edge-ordering should not
    // matter in the algorithm.
    return Vertices.update(hash_table_entry);
}


bool CdBG_Builder::process_isolated_kmer(const cuttlefish::kmer_t& kmer_hat)
{
    // Fetch the hash table entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    Vertex_Encoding& vertex_encoding = hash_table_entry.get_vertex_encoding();


    // The k-mer is already classified as a complex node.
    if(vertex_encoding.is_visited() && vertex_encoding.vertex_class() == cuttlefish::Vertex_Class::multi_in_multi_out)
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.
    
    // Classify the isolated k-mer as a complex node.
    vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::Vertex_Class::multi_in_multi_out));
    return Vertices.update(hash_table_entry);
}


// TODO: Separate the classification and the outputting task implementations into different source files.
// Reason: This source file is getting too large.


void CdBG_Builder::output_maximal_unitigs(const std::string& output_file)
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


    // Open the output file.
    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if(!output)
    {
        std::cerr << "Error opening output file " << output_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    uint32_t seqCount = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            break;


        // Look for contiguous subsequences of k-mers without the placeholder nucleotide 'N'.
        size_t kmer_idx = 0;
        while(kmer_idx <= seq_len - k)
        {
            kmer_idx = search_valid_kmer(seq, kmer_idx, seq_len - k);

            // No valid k-mer remains in the sequence.
            if(kmer_idx > seq_len - k)
                break;

            // Process a maximal valid contiguous subsequence, and advance to the index following it.
            kmer_idx = output_maximal_unitigs(seq, seq_len, kmer_idx, output);
        }
    }

    
    // Close the output file.
    if(output.fail())
        std::cerr << "Errors had been encountered for the output stream.\n";

    output.close();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


size_t CdBG_Builder::output_maximal_unitigs(const char* seq, const size_t seq_len, const size_t start_idx, std::ofstream& output)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Annotated_Kmer curr_annot_kmer(cuttlefish::kmer_t(seq, kmer_idx), kmer_idx, Vertices);

    // The k-mer is an isolated one, so is a maximal unitig by itself.
    if(kmer_idx + k == seq_len || seq[kmer_idx + k] == 'N')
        output_unitig(seq, curr_annot_kmer, curr_annot_kmer, output);
    else    // At least two adjacent k-mers are present from the index `kmer_idx`.
    {
        Annotated_Kmer prev_annot_kmer;
        Annotated_Kmer next_annot_kmer = curr_annot_kmer;
        next_annot_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);

        Annotated_Kmer unipath_start_kmer;


        // Process the first k-mer of this subsequence.

        // A maximal unitig starts at the beginning of a maximal valid subsequence.
        unipath_start_kmer = curr_annot_kmer;
        if(is_unipath_end(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, next_annot_kmer.vertex_class, next_annot_kmer.dir))
            output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);


        // Process the internal k-mers of this subsequence.
        for(kmer_idx++; kmer_idx < seq_len - k && seq[kmer_idx + k] != 'N'; ++kmer_idx)
        {
            prev_annot_kmer = curr_annot_kmer;
            curr_annot_kmer = next_annot_kmer;
            next_annot_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);


            if(is_unipath_start(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, prev_annot_kmer.vertex_class, prev_annot_kmer.dir))
                unipath_start_kmer = curr_annot_kmer;

            if(is_unipath_end(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, next_annot_kmer.vertex_class, next_annot_kmer.dir))
                output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);
        }

        prev_annot_kmer = curr_annot_kmer;
        curr_annot_kmer = next_annot_kmer;


        // Process the last k-mer of this subsequence.

        if(is_unipath_start(curr_annot_kmer.vertex_class, curr_annot_kmer.dir, prev_annot_kmer.vertex_class, prev_annot_kmer.dir))
            unipath_start_kmer = curr_annot_kmer;

        // A maximal unitig ends at the ending of a maximal valid subsequence.
        output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);
    }


    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


bool CdBG_Builder::is_unipath_start(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::kmer_dir_t dir, const cuttlefish::Vertex_Class prev_kmer_class, const cuttlefish::kmer_dir_t prev_kmer_dir) const
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


bool CdBG_Builder::is_unipath_end(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::kmer_dir_t dir, const cuttlefish::Vertex_Class next_kmer_class, const cuttlefish::kmer_dir_t next_kmer_dir) const
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


void CdBG_Builder::output_unitig(const char* seq, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, std::ofstream& output)
{
    // This is to avoid race conditions that may arise while multi-threading.
    // If two threads try to output the same unitig at the same time but
    // encounter it in the opposite orientations, then data races may arise.
    // For a particular unitig, always query the same well-defined canonical flanking
    // k-mer, irrespective of which direction the unitig may be traversed at.
    const cuttlefish::kmer_t min_flanking_kmer = std::min(start_kmer.canonical, end_kmer.canonical);
    Kmer_Hash_Entry_API hash_table_entry = Vertices[min_flanking_kmer];
    Vertex_Encoding& vertex_encoding = hash_table_entry.get_vertex_encoding();

    if(vertex_encoding.is_outputted())
        return;
    

    vertex_encoding = vertex_encoding.outputted();

    // If the hash table update is successful, only then this thread may output this unitig.
    if(Vertices.update(hash_table_entry))
        write_path(seq, start_kmer.idx, end_kmer.idx, start_kmer.kmer < end_kmer.rev_compl, output);
}


void CdBG_Builder::write_path(const char* seq, const uint32_t start_kmer_idx, const uint32_t end_kmer_idx, const bool in_forward, std::ofstream& output) const
{
    if(in_forward)
        for(uint32_t idx = start_kmer_idx; idx <= end_kmer_idx + k - 1; ++idx)
            output << seq[idx];
    else
    {
        // To avoid underflow of unsigned integers, the flanking indices are incremented by 1.
        uint32_t idx = end_kmer_idx + k;
        while(idx > start_kmer_idx)
        {
            idx--;
            output << complement(seq[idx]);
        }
    }

    output << "\n";
}


void CdBG_Builder::print_vertex_class_dist(const std::string& kmc_file_name) const
{
    Kmer_Container kmers(kmc_file_name);
    auto it_beg = kmers.begin();
    auto it_end = kmers.end();
    size_t C[4] = {0, 0, 0, 0};

    for(auto it = it_beg; it != it_end; ++it)
        C[uint8_t(Vertices[*it].decode().vertex_class)]++;


    std::cout << "Single-In-Single-Out:\t" << C[0] << "\n";
    std::cout << "Multi-In-Single-Out:\t" << C[1] << "\n";
    std::cout << "Single-In-Multi-Out:\t" << C[2] << "\n";
    std::cout << "Multi-In-Multi-Out:\t" << C[3] << "\n";
}
