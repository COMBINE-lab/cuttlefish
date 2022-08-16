
#include "CdBG.hpp"
#include "DNA_Utility.hpp"
#include "Directed_Kmer.hpp"
#include "Ref_Parser.hpp"
#include "Thread_Pool.hpp"

#include <iomanip>
#include <chrono>


template <uint16_t k>
void CdBG<k>::classify_vertices()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const std::string& buckets_file_path = params.buckets_file_path();

    // The serialized hash table buckets (saved from some earlier execution) exists.
    if(!buckets_file_path.empty() && file_exists(buckets_file_path))
    {
        std::cout << "Found the hash table buckets at file " << buckets_file_path << "\n";
        std::cout << "Loading the buckets.\n";

        hash_table->load_hash_buckets(buckets_file_path);

        std::cout << "Loaded the buckets into memory.\n";
    }
    else    // No buckets file name provided, or does not exist. Build and save (if specified) one now.
    {
        // Open a parser for the FASTA / FASTQ file containing the reference.
        Ref_Parser parser(params.sequence_input());


        // Construct a thread pool.
        const uint16_t thread_count = params.thread_count();
        Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::classification);


        // Track the maximum sequence buffer size used and the total length of the references.
        size_t max_buf_sz = 0;
        uint64_t ref_len = 0;
        uint64_t seq_count = 0;

        // Parse sequences one-by-one, and continue partial classification of the k-mers through them.
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
            {
                short_refs.emplace_back(std::make_pair(parser.seq_name(), seq_len));
                continue;
            }


            // Single-threaded classification.
            // process_substring(seq, seq_len, 0, seq_len - k);


            // Multi-threaded classification.
            distribute_classification(seq, seq_len, thread_pool);
            thread_pool.wait_completion();
        }

        std::cerr << "\nProcessed " << seq_count << " sequences. Total reference length: " << ref_len << " bases.\n";
        std::cout << "Maximum input sequence buffer size used: " << max_buf_sz / (1024 * 1024) << " MB.\n";

        // Close the thread-pool.
        thread_pool.close();

        // Close the parser.
        parser.close();


        // Save the hash table buckets, if a file path is provided.
        if(params.save_buckets())
        {
            hash_table->save_hash_buckets(buckets_file_path);
            std::cout << "Saved the hash buckets at " << buckets_file_path << "\n";
        }
    }


    dbg_info.add_basic_info(*this);

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done computing the vertex-states. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void CdBG<k>::distribute_classification(const char* seq, const size_t seq_len, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();
    const size_t task_size = (seq_len - k + 1) / thread_count;
    const uint16_t partition_count = (task_size < PARTITION_SIZE_THRESHOLD ? 1 : thread_count);

    size_t left_end = 0;
    size_t right_end;

    for(uint16_t t_id = 0; t_id < partition_count; ++t_id)
    {
        right_end = (t_id == partition_count - 1 ? seq_len - k : left_end + task_size - 1);

        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_classification_task(idle_thread_id, seq, seq_len, left_end, right_end);

        left_end += task_size;
    }
}


template <uint16_t k> 
void CdBG<k>::process_substring(const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
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


template <uint16_t k> 
size_t CdBG<k>::process_contiguous_subseq(const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Directed_Kmer<k> curr_kmer(Kmer<k>(seq, kmer_idx));

    // The subsequence contains only an isolated k-mer,
    // i.e. there's no valid left or right neighboring k-mer to this k-mer.
    if((kmer_idx == 0 || DNA_Utility::is_placeholder(seq[kmer_idx - 1])) &&
        (kmer_idx + k == seq_len || DNA_Utility::is_placeholder(seq[kmer_idx + k])))
        while(!process_isolated_kmer(curr_kmer));
    else    // At least one valid neighbor exists, either to the left or to the right, or on both sides.
    {
        // Process the leftmost k-mer of this contiguous subsequence.

        // No valid right neighbor exists for the k-mer.
        if(kmer_idx + k == seq_len || DNA_Utility::is_placeholder(seq[kmer_idx + k]))
        {
            // A valid left neighbor exists at it's not an isolated k-mer.
            while(!process_rightmost_kmer(curr_kmer, seq[kmer_idx - 1]));

            // The contiguous sequence ends at this k-mer.
            return kmer_idx + k;
        }

        // A valid right neighbor exists for the k-mer.
        Directed_Kmer<k> next_kmer = curr_kmer;
        next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);
        
        // No valid left neighbor exists for the k-mer.
        if(kmer_idx == 0 || DNA_Utility::is_placeholder(seq[kmer_idx - 1]))
            while(!process_leftmost_kmer(curr_kmer, next_kmer, seq[kmer_idx + k]));
        // Both left and right valid neighbors exist for this k-mer.
        else
            while(!process_internal_kmer(curr_kmer, next_kmer, seq[kmer_idx - 1], seq[kmer_idx + k]));
        

        // Process the internal k-mers of this contiguous subsequence.
        // Each of these k-mers have valid neighbors to their left and right.
        for(kmer_idx++; kmer_idx < right_end && !DNA_Utility::is_placeholder(seq[kmer_idx + k]); ++kmer_idx)
        {
            curr_kmer = next_kmer;
            next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

            while(!process_internal_kmer(curr_kmer, next_kmer, seq[kmer_idx - 1], seq[kmer_idx + k]));
        }


        // Process the rightmost k-mer of this contiguous subsequence. This does not coincide with the leftmost k-mer.

        if(kmer_idx <= right_end)   // Required for the cases where the provided range has just length 1, so `start_idx = right_end`.
        {
            curr_kmer = next_kmer;
        
            // No valid right neighbor exists for the k-mer.
            if(kmer_idx + k == seq_len || DNA_Utility::is_placeholder(seq[kmer_idx + k]))
                while(!process_rightmost_kmer(curr_kmer, seq[kmer_idx - 1]));
            // A valid right neighbor exists for the k-mer.
            else
            {
                next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

                while(!process_internal_kmer(curr_kmer, next_kmer, seq[kmer_idx - 1], seq[kmer_idx + k]));
            }
        }
        else
            kmer_idx--; // `kmer_idx` has to be the index of the last valid k-mer encountered.
    }


    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


template <uint16_t k> 
bool CdBG<k>::is_self_loop(const Kmer<k>& kmer_hat, const Kmer<k>& next_kmer_hat) const
{
    return kmer_hat == next_kmer_hat;
}


template <uint16_t k>
bool CdBG<k>::process_loop(const Directed_Kmer<k>& kmer, const Directed_Kmer<k>& next_kmer, const char prev_char)
{
    // Note that, any loop that connects two different sides of a vertex makes it a
    // complex node. This is because, from whichever side you may try to include this
    // vertex into a non-trivial unitig, it violates the unitig-defining constraints on
    // that side due to that loop edge incident there. The only way that a vertex having
    // a loop incident to it may get included into a non-trivial unitig is that the loop
    // is totally contained at one of its sides. If the opposite has only one distinct
    // edge, then the vertex may get into a non-trivial unitig through that edge. In such
    // cases, it also becomes an endpoint for that unitig.

    // Gist: if a vertex has a loop that crosses its sides (due to an immediate direct repeat
    // in the underlying sequence), then the vertex is a complex node. Otherwise, the loop is
    // entirely contained at one of its sides (due to an immediate inverted repeat). In this
    // case, the vertex might get included into a maximal unitig through its other side. The
    // side with the loop blocks the vertex to extend the unitig any farther.


    // It is the leftmost k-mer in some sequence. So, its left side can't have incident edges,
    // and the loop (crossing or one-sided) blocks the other side for any unitigs-inclusion.
    // Or, the label encountered for this and the next k-mer is seen as the same. So there is
    // a direct repeat, that produces a crossing loop. In either case, this is a complex node.
    if(!prev_char || kmer.kmer() == next_kmer.kmer())
    {
        // Fetch the entry for `kmer_hat`.
        const Kmer<k>& kmer_hat = kmer.canonical();
        Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER> hash_table_entry = hash_table->at(kmer_hat);
        State& state = hash_table_entry.get_state();
        state = State(Vertex(cuttlefish::State_Class::multi_in_multi_out));

        return hash_table->update(hash_table_entry);
    }

    
    // The k-mer is internal, and the loop is one-sided. So it not possible to extend a maximal
    // unitig through that side, so this k-mer can equivalently be treated as a rightmost k-mer
    // (a sentinel) of some sequence.
    return process_rightmost_kmer(kmer, prev_char);
}


template <uint16_t k> 
bool CdBG<k>::process_leftmost_kmer(const Directed_Kmer<k>& kmer, const Directed_Kmer<k>& next_kmer, const char next_char)
{
    const Kmer<k>& kmer_hat = kmer.canonical();
    const cuttlefish::dir_t dir = kmer.dir();
    const Kmer<k>& next_kmer_hat = next_kmer.canonical();

    // Fetch the entry for `kmer_hat`.
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER> hash_table_entry = hash_table->at(kmer_hat);
    State& state = hash_table_entry.get_state();

    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.

    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
        return process_loop(kmer, next_kmer);


    const State old_state = state;
    const cuttlefish::base_t next_base = DNA_Utility::map_base(next_char);

    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::State_Class::multi_in_single_out, next_base));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.state_class_ == cuttlefish::State_Class::single_in_single_out)
            {
                if(vertex.back_ == next_base)
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_single_out;
                else
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.state_class_ == cuttlefish::State_Class::multi_in_single_out)
            {
                if(vertex.back_ != next_base)
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.state_class_ == cuttlefish::State_Class::single_in_multi_out
            {
                vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state = State(vertex);
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::State_Class::single_in_multi_out, DNA_Utility::complement(next_base)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.state_class_ == cuttlefish::State_Class::single_in_single_out)
            {
                if(vertex.front_ == DNA_Utility::complement(next_base))
                    vertex.state_class_ = cuttlefish::State_Class::single_in_multi_out;
                else
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.state_class_ == cuttlefish::State_Class::multi_in_single_out)
            {
                vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else    // vertex.state_class_ == cuttlefish::State_Class::single_in_multi_out
            {
                if(vertex.front_ != DNA_Utility::complement(next_base))
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }


    // We can get away without updating the same value again, because -- (1) even if this k-mer's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    return state == old_state ? true : hash_table->update(hash_table_entry);
}


template <uint16_t k> 
bool CdBG<k>::process_rightmost_kmer(const Directed_Kmer<k>& kmer, const char prev_char)
{
    const Kmer<k>& kmer_hat = kmer.canonical();
    const cuttlefish::dir_t dir = kmer.dir();

    // Fetch the entry for `kmer_hat`.
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER> hash_table_entry = hash_table->at(kmer_hat);
    State& state = hash_table_entry.get_state();

    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.


    const State old_state = state;
    const cuttlefish::base_t prev_base = DNA_Utility::map_base(prev_char);


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::State_Class::single_in_multi_out, prev_base));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.state_class_ == cuttlefish::State_Class::single_in_single_out)
            {
                if(vertex.front_ == prev_base)
                    vertex.state_class_ = cuttlefish::State_Class::single_in_multi_out;
                else
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.state_class_ == cuttlefish::State_Class::multi_in_single_out)
            {
                vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;
                
                state = State(vertex);
            }
            else    // vertex.state_class_ == cuttlefish::State_Class::single_in_multi_out
            {
                if(vertex.front_ != prev_base)
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::State_Class::multi_in_single_out, DNA_Utility::complement(prev_base)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.state_class_ == cuttlefish::State_Class::single_in_single_out)
            {
                if(vertex.back_ == DNA_Utility::complement(prev_base))
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_single_out;
                else
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state= State(vertex);
            }
            else if(vertex.state_class_ == cuttlefish::State_Class::multi_in_single_out)
            {
                if(vertex.back_ != DNA_Utility::complement(prev_base))
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.state_class_ == cuttlefish::State_Class::single_in_multi_out
            {
                vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                state = State(vertex);
            }
        }
    }


    // We can get away without updating the same value again, because -- (1) even if this k-mer's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    return state == old_state ? true : hash_table->update(hash_table_entry);
}


template <uint16_t k> 
bool CdBG<k>::process_internal_kmer(const Directed_Kmer<k>& kmer, const Directed_Kmer<k>& next_kmer, const char prev_char, const char next_char)
{
    const Kmer<k>& kmer_hat = kmer.canonical();
    const cuttlefish::dir_t dir = kmer.dir();
    const Kmer<k>& next_kmer_hat = next_kmer.canonical();

    // Fetch the hash table entry for `kmer_hat`.
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER> hash_table_entry = hash_table->at(kmer_hat);
    State& state = hash_table_entry.get_state();

    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.

    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
        return process_loop(kmer, next_kmer, prev_char);

    
    const State old_state = state;
    const cuttlefish::base_t prev_base = DNA_Utility::map_base(prev_char);
    const cuttlefish::base_t next_base = DNA_Utility::map_base(next_char);

    if(dir == cuttlefish::FWD)
    {
        // The k-mer is encountered for the first time, and in the forward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::State_Class::single_in_single_out, prev_base, next_base));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.state_class_ == cuttlefish::State_Class::single_in_single_out)
            {
                if(vertex.front_ == prev_base && vertex.back_ == next_base)
                    return true;    // See note at the end on early return w/o hash table update.
                else if(vertex.front_ != prev_base && vertex.back_ != next_base)
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;
                else if(vertex.front_ != prev_base)
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_single_out;
                else    // vertex.back_ != next_base
                    vertex.state_class_ = cuttlefish::State_Class::single_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.state_class_ == cuttlefish::State_Class::multi_in_single_out)
            {
                if(vertex.back_ != next_base)
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.state_class_ == cuttlefish::State_Class::single_in_multi_out
            {
                if(vertex.front_ != prev_base)
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }
    else
    {
        // The k-mer is encountered for the first time, and in the backward direction.
        if(!state.is_visited())
            state = State(Vertex(cuttlefish::State_Class::single_in_single_out, DNA_Utility::complement(next_base), DNA_Utility::complement(prev_base)));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = state.decode();

            if(vertex.state_class_ == cuttlefish::State_Class::single_in_single_out)
            {
                if(vertex.front_ == DNA_Utility::complement(next_base) && vertex.back_ == DNA_Utility::complement(prev_base))
                    return true;    // See note at the end on early return w/o hash table update.
                else if(vertex.front_ != DNA_Utility::complement(next_base) && vertex.back_ != DNA_Utility::complement(prev_base))
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;
                else if(vertex.front_ != DNA_Utility::complement(next_base))
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_single_out;
                else    // vertex.back_ != complement(prev_base)
                    vertex.state_class_ = cuttlefish::State_Class::single_in_multi_out;

                state = State(vertex);
            }
            else if(vertex.state_class_ == cuttlefish::State_Class::multi_in_single_out)
            {
                if(vertex.back_ != DNA_Utility::complement(prev_base))
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
            else    // vertex.state_class_ == cuttlefish::State_Class::single_in_multi_out
            {
                if(vertex.front_ != DNA_Utility::complement(next_base))
                {
                    vertex.state_class_ = cuttlefish::State_Class::multi_in_multi_out;

                    state = State(vertex);
                }
            }
        }
    }


    // We can get away without updating the same value again, because -- (1) even if this k-mer's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    return state == old_state ? true : hash_table->update(hash_table_entry);
}


template <uint16_t k> 
bool CdBG<k>::process_isolated_kmer(const Directed_Kmer<k>& kmer)
{
    const Kmer<k>& kmer_hat = kmer.canonical();

    // Fetch the hash table entry for `kmer_hat`.
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER> hash_table_entry = hash_table->at(kmer_hat);
    State& state = hash_table_entry.get_state();


    // The k-mer is already classified as a complex node.
    if(state.is_dead_end())
        return true;    // Early return is safe from here as this vertex is a dead-end for state transitions.
    
    // Classify the isolated k-mer as a complex node.
    state = State(Vertex(cuttlefish::State_Class::multi_in_multi_out));
    return hash_table->update(hash_table_entry);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
