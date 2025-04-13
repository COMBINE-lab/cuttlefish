
#include "CdBG.hpp"
#include "DNA_Utility.hpp"
#include "Annotated_Kmer.hpp"


template <uint16_t k>
void CdBG<k>::output_plain_off_substring(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
{
    size_t kmer_idx = left_end;
    while(kmer_idx <= right_end)
    {
        kmer_idx = search_valid_kmer(seq, kmer_idx, right_end);

        // No valid k-mer remains in the sequence.
        if(kmer_idx > right_end)
            break;

        // Process a maximal valid contiguous subsequence, and advance to the index following it.
        kmer_idx = output_maximal_unitigs_plain(thread_id, seq, seq_len, right_end, kmer_idx);
    }
}


template <uint16_t k>
size_t CdBG<k>::output_maximal_unitigs_plain(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx)
{
    size_t kmer_idx = start_idx;

    // assert(kmer_idx <= seq_len - k);

    Annotated_Kmer<k> curr_kmer(Kmer<k>(seq, kmer_idx), kmer_idx, *hash_table);

    // The subsequence contains only an isolated k-mer, i.e. there's no valid left or right
    // neighboring k-mer to this k-mer. So it's a maximal unitig by itself.
    if((kmer_idx == 0 || DNA_Utility::is_placeholder(seq[kmer_idx - 1])) &&
        (kmer_idx + k == seq_len || DNA_Utility::is_placeholder(seq[kmer_idx + k])))
        output_plain_unitig(thread_id, seq, curr_kmer, curr_kmer);
    else    // At least one valid neighbor exists, either to the left or to the right, or on both sides.
    {
        // No valid right neighbor exists for the k-mer.
        if(kmer_idx + k == seq_len || DNA_Utility::is_placeholder(seq[kmer_idx + k]))
        {
            // A valid left neighbor exists as it's not an isolated k-mer.
            Annotated_Kmer<k> prev_kmer(Kmer<k>(seq, kmer_idx - 1), kmer_idx, *hash_table);
            
            if(is_unipath_start(curr_kmer.state_class(), curr_kmer.dir(), prev_kmer.state_class(), prev_kmer.dir()))
                // A maximal unitig ends at the ending of a maximal valid subsequence.
                output_plain_unitig(thread_id, seq, curr_kmer, curr_kmer);

            // The contiguous sequence ends at this k-mer.
            return kmer_idx + k;
        }


        // A valid right neighbor exists for the k-mer.
        Annotated_Kmer<k> next_kmer = curr_kmer;
        next_kmer.roll_to_next_kmer(seq[kmer_idx + k], *hash_table);

        bool on_unipath = false;
        Annotated_Kmer<k> unipath_start_kmer;
        Annotated_Kmer<k> prev_kmer;

        // No valid left neighbor exists for the k-mer.
        if(kmer_idx == 0 || DNA_Utility::is_placeholder(seq[kmer_idx - 1]))
        {
            // A maximal unitig starts at the beginning of a maximal valid subsequence.
            on_unipath = true;
            unipath_start_kmer = curr_kmer;
        }
        // Both left and right valid neighbors exist for this k-mer.
        else
        {
            prev_kmer = Annotated_Kmer<k>(Kmer<k>(seq, kmer_idx - 1), kmer_idx, *hash_table);
            if(is_unipath_start(curr_kmer.state_class(), curr_kmer.dir(), prev_kmer.state_class(), prev_kmer.dir()))
            {
                on_unipath = true;
                unipath_start_kmer = curr_kmer;
            }
        }

        if(on_unipath && is_unipath_end(curr_kmer.state_class(), curr_kmer.dir(), next_kmer.state_class(), next_kmer.dir()))
        {
            output_plain_unitig(thread_id, seq, unipath_start_kmer, curr_kmer);
            on_unipath = false;
        }


        // Process the rest of the k-mers of this contiguous subsequence.
        for(kmer_idx++; on_unipath || kmer_idx <= right_end; ++kmer_idx)
        {
            prev_kmer = curr_kmer;
            curr_kmer = next_kmer;

            if(is_unipath_start(curr_kmer.state_class(), curr_kmer.dir(), prev_kmer.state_class(), prev_kmer.dir()))
            {
                on_unipath = true;
                unipath_start_kmer = curr_kmer;
            }


            // No valid right neighbor exists for the k-mer.
            if(kmer_idx + k == seq_len || DNA_Utility::is_placeholder(seq[kmer_idx + k]))
            {
                // A maximal unitig ends at the ending of a maximal valid subsequence.
                if(on_unipath)
                {
                    output_plain_unitig(thread_id, seq, unipath_start_kmer, curr_kmer);
                    on_unipath = false;
                }

                // The contiguous sequence ends at this k-mer.
                return kmer_idx + k;
            }
            else    // A valid right neighbor exists.
            {
                next_kmer.roll_to_next_kmer(seq[kmer_idx + k], *hash_table);
                
                if(on_unipath && is_unipath_end(curr_kmer.state_class(), curr_kmer.dir(), next_kmer.state_class(), next_kmer.dir()))
                {
                    output_plain_unitig(thread_id, seq, unipath_start_kmer, curr_kmer);
                    on_unipath = false;
                }
            }
        }
    }
    
    
    // Return the non-inclusive ending index of the processed contiguous subsequence.
    return kmer_idx + k;
}


template <uint16_t k>
void CdBG<k>::output_plain_unitig(const uint16_t thread_id, const char* const seq, const Annotated_Kmer<k>& start_kmer, const Annotated_Kmer<k>& end_kmer)
{
    // This is to avoid race conditions that may arise while multi-threading.
    // If two threads try to output the same unitig at the same time but
    // encounter it in the opposite orientations, then data races may arise.
    // For a particular unitig, always query the same well-defined canonical flanking
    // k-mer, irrespective of which direction the unitig may be traversed at.
    const Kmer<k> min_flanking_kmer = std::min(start_kmer.canonical(), end_kmer.canonical());
    const uint64_t bucket_id = hash_table->bucket_id(min_flanking_kmer);
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER> hash_table_entry = hash_table->at(bucket_id);
    State& state = hash_table_entry.get_state();

    if(state.is_outputted())
        return;
    

    state = state.outputted();

    // If the hash table update is successful, only then this thread may output this unitig.
    if(hash_table->update(hash_table_entry))
    {
        write_path(thread_id, seq, bucket_id, start_kmer.idx(), end_kmer.idx(), start_kmer.kmer() < end_kmer.rev_compl());
        
        unipaths_info_local[thread_id].add_maximal_unitig(end_kmer.idx() - start_kmer.idx() + 1);
    }
}


template <uint16_t k>
void CdBG<k>::write_path(const uint16_t thread_id, const char* const seq, const uint64_t unitig_id, const size_t start_kmer_idx, const size_t end_kmer_idx, const cuttlefish::dir_t dir) 
{
    std::string& buffer = output_buffer[thread_id];
    const size_t path_len = end_kmer_idx - start_kmer_idx + k;
    constexpr std::size_t header_len = 12;  // FASTA header len: '>' + <id> + '\n'


    ensure_buffer_space(buffer, path_len + header_len, output_[thread_id]);

    buffer += ">";
    buffer += fmt::format_int(unitig_id).c_str();
    buffer += "\n";

    if(dir == cuttlefish::FWD)
        for(size_t offset = 0; offset < path_len; ++offset)
            buffer += DNA_Utility::upper(seq[start_kmer_idx + offset]);
    else    // dir == cuttlefish::BWD
        for(size_t offset = 0; offset < path_len; ++offset)
            buffer += DNA_Utility::complement(seq[end_kmer_idx + k - 1 - offset]);

    // End the path.
    buffer += "\n";


    if(kmer_idx != nullptr)
    {
        const char* const unitig_seq = buffer.data() + (buffer.size() - (path_len + 1)); // The +1 length is to account for the ending line-break.
        kmer_idx->deposit(token[thread_id], unitig_seq, path_len);
    }

    
    // Mark buffer size increment.
    check_output_buffer(thread_id);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, CdBG)
