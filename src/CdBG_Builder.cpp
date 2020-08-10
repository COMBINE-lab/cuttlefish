
#include "CdBG.hpp"
#include "kseq/kseq.h"
#include "Kmer_Iterator.hpp"

#include <fstream>
#include <thread>
#include <cassert>
#include <chrono>
#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


void CdBG::classify_vertices()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
    const std::string& ref_file_path =  params.ref_file_path();
    FILE* const input = fopen(ref_file_path.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    
    // Initialize the parser.
    kseq_t* const parser = kseq_init(fileno(input));

    // Track the maximum sequence buffer size used.
    size_t max_buf_sz = 0;

    // Parse sequences one-by-one, and continue partial classification of the k-mers through them.
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


        // Single-threaded classification.
        // process_substring(seq, seq_len, 0, seq_len - k);


        // Multi-threaded classification.
        const uint16_t thread_count = params.thread_count();
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
                task.emplace_back(&CdBG::process_substring, this, seq, seq_len, left_end, right_end);
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

    std::cout << "Maximum buffer size used (in MB): " << max_buf_sz / (1024 * 1024) << "\n";


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done classifying the vertices. Time taken = " << elapsed_seconds << " seconds.\n";
}


void CdBG::process_substring(const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
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


size_t CdBG::process_contiguous_subseq(const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Directed_Kmer curr_kmer(cuttlefish::kmer_t(seq, kmer_idx));

    // The subsequence contains only an isolated k-mer,
    // i.e. there's no valid left or right neighboring k-mer to this k-mer.
    if((kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE) &&
        (kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE))
        while(!process_isolated_kmer(curr_kmer.canonical()));
    else    // At least one valid neighbor exists, either to the left or to the right, or on both sides.
    {
        // Process the leftmost k-mer of this contiguous subsequence.

        // No valid right neighbor exists for the k-mer.
        if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
        {
            // A valid left neighbor exists at it's not an isolated k-mer.
            while(!process_rightmost_kmer(curr_kmer.canonical(), curr_kmer.dir(), seq[kmer_idx - 1]));

            // The contiguous sequence ends at this k-mer.
            return kmer_idx + k;
        }

        // A valid right neighbor exists for the k-mer.
        Directed_Kmer next_kmer = curr_kmer;
        next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);
        
        // No valid left neighbor exists for the k-mer.
        if(kmer_idx == 0 || seq[kmer_idx - 1] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
            while(!process_leftmost_kmer(curr_kmer.canonical(), curr_kmer.dir(), next_kmer.canonical(), seq[kmer_idx + k]));
        // Both left and right valid neighbors exist for this k-mer.
        else
            while(!process_internal_kmer(curr_kmer.canonical(), curr_kmer.dir(), next_kmer.canonical(), seq[kmer_idx - 1], seq[kmer_idx + k]));
        

        // Process the internal k-mers of this contiguous subsequence.
        // Each of these k-mers have valid neighbors to their left and right.
        for(kmer_idx++; kmer_idx < right_end && seq[kmer_idx + k] != cuttlefish::PLACEHOLDER_NUCLEOTIDE; ++kmer_idx)
        {
            curr_kmer = next_kmer;
            next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

            while(!process_internal_kmer(curr_kmer.canonical(), curr_kmer.dir(), next_kmer.canonical(), seq[kmer_idx - 1], seq[kmer_idx + k]));
        }


        // Process the rightmost k-mer of this contiguous subsequence. This does not coincide with the leftmost k-mer.

        if(kmer_idx <= right_end)   // Required for the cases where the provided range has just length 1, so `start_idx = right_end`.
        {
            curr_kmer = next_kmer;
        
            // No valid right neighbor exists for the k-mer.
            if(kmer_idx + k == seq_len || seq[kmer_idx + k] == cuttlefish::PLACEHOLDER_NUCLEOTIDE)
                while(!process_rightmost_kmer(curr_kmer.canonical(), curr_kmer.dir(), seq[kmer_idx - 1]));
            // A valid right neighbor exists for the k-mer.
            else
            {
                next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

                while(!process_internal_kmer(curr_kmer.canonical(), curr_kmer.dir(), next_kmer.canonical(), seq[kmer_idx - 1], seq[kmer_idx + k]));
            }
        }
        else
            kmer_idx--; // `kmer_idx` has to be the index of the last valid k-mer encountered.
    }


    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


bool CdBG::is_self_loop(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_t& next_kmer_hat) const
{
    return kmer_hat == next_kmer_hat;
}


bool CdBG::process_leftmost_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t next_nucl)
{
    // Fetch the entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    State& state = hash_table_entry.get_state();
    const State old_state = state;


    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
        state = State(Vertex(cuttlefish::Vertex_Class::multi_in_multi_out));
    else if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::Vertex_Class::multi_in_single_out, next_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.vertex_class_ == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.exit_ == next_nucl)
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_single_out;
                else
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.vertex_class_ == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit_ != next_nucl)
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state = State(vertex);
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::Vertex_Class::single_in_multi_out, cuttlefish::kmer_t::complement(next_nucl)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.vertex_class_ == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter_ == cuttlefish::kmer_t::complement(next_nucl))
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::single_in_multi_out;
                else
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.vertex_class_ == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter_ != cuttlefish::kmer_t::complement(next_nucl))
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }


    // We can get away without updating the same value again, because -- (1) even if this k-mer's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    return state == old_state ? true : Vertices.update(hash_table_entry);
}


bool CdBG::process_rightmost_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::dir_t dir, const cuttlefish::nucleotide_t prev_nucl)
{
    // Fetch the entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    State& state = hash_table_entry.get_state();
    const State old_state = state;


    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::Vertex_Class::single_in_multi_out, prev_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.vertex_class_ == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter_ == prev_nucl)
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::single_in_multi_out;
                else
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.vertex_class_ == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;
                
                state = State(vertex);
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter_ != prev_nucl)
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::Vertex_Class::multi_in_single_out, cuttlefish::kmer_t::complement(prev_nucl)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.vertex_class_ == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.exit_ == cuttlefish::kmer_t::complement(prev_nucl))
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_single_out;
                else
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state= State(vertex);
            }
            else if(vertex.vertex_class_ == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit_ != cuttlefish::kmer_t::complement(prev_nucl))
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                state = State(vertex);
            }
        }
    }


    // We can get away without updating the same value again, because -- (1) even if this k-mer's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    return state == old_state ? true : Vertices.update(hash_table_entry);
}


bool CdBG::process_internal_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t prev_nucl, const cuttlefish::nucleotide_t next_nucl)
{
    // Fetch the hash table entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    State& state = hash_table_entry.get_state();
    const State old_state = state;


    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.


    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
        state = State(Vertex(cuttlefish::Vertex_Class::multi_in_multi_out));
    else if(dir == cuttlefish::FWD)
    {
        // The k-mer is encountered for the first time, and in the forward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::Vertex_Class::single_in_single_out, prev_nucl, next_nucl));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.vertex_class_ == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter_ == prev_nucl && vertex.exit_ == next_nucl)
                    return true;    // See note at the end on early return w/o hash table update.
                else if(vertex.enter_ != prev_nucl && vertex.exit_ != next_nucl)
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;
                else if(vertex.enter_ != prev_nucl)
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_single_out;
                else    // vertex.exit != next_nucl
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::single_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.vertex_class_ == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit_ != next_nucl)
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter_ != prev_nucl)
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }
    else
    {
        // The k-mer is encountered for the first time, and in the backward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::Vertex_Class::single_in_single_out, cuttlefish::kmer_t::complement(next_nucl), cuttlefish::kmer_t::complement(prev_nucl)));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.vertex_class_ == cuttlefish::Vertex_Class::single_in_single_out)
            {
                if(vertex.enter_ == cuttlefish::kmer_t::complement(next_nucl) && vertex.exit_ == cuttlefish::kmer_t::complement(prev_nucl))
                    return true;    // See note at the end on early return w/o hash table update.
                else if(vertex.enter_ != cuttlefish::kmer_t::complement(next_nucl) && vertex.exit_ != cuttlefish::kmer_t::complement(prev_nucl))
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;
                else if(vertex.enter_ != cuttlefish::kmer_t::complement(next_nucl))
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_single_out;
                else    // vertex.exit != complement(prev_nucl)
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::single_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.vertex_class_ == cuttlefish::Vertex_Class::multi_in_single_out)
            {
                if(vertex.exit_ != cuttlefish::kmer_t::complement(prev_nucl))
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.vertex_class == cuttlefish::Vertex_Class::single_in_multi_out
            {
                if(vertex.enter_ != cuttlefish::kmer_t::complement(next_nucl))
                {
                    vertex.vertex_class_ = cuttlefish::Vertex_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }


    // We can get away without updating the same value again, because -- (1) even if this k-mer's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    return state == old_state ? true : Vertices.update(hash_table_entry);
}


bool CdBG::process_isolated_kmer(const cuttlefish::kmer_t& kmer_hat)
{
    // Fetch the hash table entry for `kmer_hat`.
    Kmer_Hash_Entry_API hash_table_entry = Vertices[kmer_hat];
    State& state = hash_table_entry.get_state();


    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.
    
    // Classify the isolated k-mer as a complex node.
    state = State(Vertex(cuttlefish::Vertex_Class::multi_in_multi_out));
    return Vertices.update(hash_table_entry);
}


void CdBG::print_vertex_class_dist() const
{
    const std::string& kmc_db_path = params.kmc_db_path();

    Kmer_Container kmers(kmc_db_path);
    auto it_beg = kmers.begin();
    auto it_end = kmers.end();
    size_t C[4] = {0, 0, 0, 0};

    for(auto it = it_beg; it != it_end; ++it)
        C[uint8_t(Vertices[*it].decode().vertex_class_)]++;


    std::cout << "Single-In-Single-Out:\t" << C[0] << "\n";
    std::cout << "Multi-In-Single-Out:\t" << C[1] << "\n";
    std::cout << "Single-In-Multi-Out:\t" << C[2] << "\n";
    std::cout << "Multi-In-Multi-Out:\t" << C[3] << "\n";
}
