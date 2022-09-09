
#include "Validator.hpp"
#include "DNA_Utility.hpp"
#include "Ref_Parser.hpp"
#include "Kmer_Container.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <fstream>


template <uint16_t k>
void Validator<k>::validate_sequence_completion(bool& result)
{
    console->info("Testing validation of the completeness of coverage of the sequence by the produced unitigs.\n");

    const std::string& kmc_db_path = params.kmc_db_path();
    const uint16_t thread_count = params.thread_count();

    const Kmer_Container<k> kmer_container(kmc_db_path);
    const uint64_t kmer_count = kmer_container.size();

    console->info("Number of k-mers in the k-mer database: {}\n", kmer_count);


    // Allocate the unitig tables.
    console->info("Allocating the unitig tables.\n");

    unitig_id.resize(kmer_count);
    unitig_dir.resize(kmer_count);

    console->info("Done allocation of the unitig tables.\n");


    // Load the unitigs into memory and build the associated tables `unitig_id` and `unitig_dir`.
    build_unitig_tables();


    // Open a parser for the FASTA / FASTQ file containing the reference.
    Ref_Parser parser(params.reference_input());


    std::vector<std::thread> th(thread_count);  // Thread-pool (round-robin) to validate the sequences parallelly.
    uint32_t seqCount = 0;  // Number of sequences read.
    bool success = true;    // Whether the total validation succeeded.
    char** S = new char*[thread_count]; // Sequence buffers to pass on to the threads.
    bool* thread_result = new bool[thread_count];   // Validation result produced by the threads.
    
    // Parse sequences one-by-one, and continue spelling them using the resultant unitigs of the compaction algorithm.
    while(parser.read_next_seq())
    {
        const char* const seq = parser.seq();
        const size_t seq_len = parser.seq_len();

        console->info("Spelling out sequence {}, with length {}.\n", seqCount, seq_len);


        // ID of the thread to pass on the latest-read sequence `seq`.
        uint16_t thread_id = seqCount % thread_count;

        // A buffer exists for thread number `thread_id`, so it's already processing a sequence.
        if(seqCount >= thread_count)
        {
            if(!th[thread_id].joinable())
            {
                console->error("Early termination faced for a worker thread. Aborting.\n");
                std::exit(EXIT_FAILURE);
            }

            th[thread_id].join();
            delete S[thread_id];
            
            success = success && thread_result[thread_id];
            if(!success)
                break;
        }

        // Allocate a buffer for the thread number `thread_id` to pass on the latest-read sequence `seq`.
        S[thread_id] = new char[seq_len + 1];
        strcpy(S[thread_id], seq);

        th[thread_id] = std::thread(&Validator::walk_sequence, this, S[thread_id], seq_len, std::ref(thread_result[thread_id]));

        seqCount++;
    }


    // Finish up the threads that might still be running.
    for(uint16_t thread_id = 0; thread_id < thread_count; ++thread_id)
        if(th[thread_id].joinable())
        {
            th[thread_id].join();
            delete S[thread_id];

            success = success && thread_result[thread_id];
        }


    // Close the parser.
    parser.close();

    result = success;
}


template <uint16_t k>
void Validator<k>::build_unitig_tables()
{
    const std::string& cdbg_file_path = params.cdbg_file_path();

    std::ifstream input(cdbg_file_path.c_str(), std::ifstream::in);
    if(!input)
    {
        console->error("Error reading from the unitigs file {}. Aborting.\n", cdbg_file_path);
        std::exit(EXIT_FAILURE);
    }


    console->info("Loading the unitigs from the file {}\n", cdbg_file_path);

    size_t unitig_count = 0;
    std::string unitig;
    while(input >> unitig)
    {
        U.emplace_back(unitig);

        if(unitig.length() == k)
        {
            const Kmer<k> kmer(unitig);
            const uint64_t kmer_hash = mph->lookup(kmer.canonical());
            
            unitig_id[kmer_hash] = unitig_count;
            unitig_dir[kmer_hash] = Unitig_Dir::either;
        }
        else
        {
            const Kmer<k> last_kmer(unitig.substr(unitig.length() - k, k));
            const Kmer<k> kmer_fwd(unitig, 0);
            const Kmer<k> kmer_bwd = last_kmer.reverse_complement();

            const uint64_t kmer_fwd_hash = mph->lookup(kmer_fwd.canonical());
            const uint64_t kmer_bwd_hash = mph->lookup(kmer_bwd.canonical());
            
            unitig_id[kmer_fwd_hash] = unitig_count;
            unitig_dir[kmer_fwd_hash] = Unitig_Dir::fwd;

            unitig_id[kmer_bwd_hash] = unitig_count;
            unitig_dir[kmer_bwd_hash] = Unitig_Dir::bwd;
        }

        
        unitig_count++;
        
        // Track progress.
        if(unitig_count % PROGRESS_GRAIN_SIZE == 0)
            console->info("Loaded {}M unitigs.\n", unitig_count / 1000000);
    }

    console->info("Done loading a total of {} unitigs.\n", unitig_count);

    input.close();
}


template <uint16_t k>
void Validator<k>::walk_sequence(const char* const seq, const size_t seq_len, bool& result) const
{
    // Nothing to process for sequences with length shorter than `k`.
    if(seq_len < k)
        return;

    size_t kmer_idx = 0;
    while(kmer_idx <= seq_len - k)
    {
        kmer_idx = search_valid_kmer(seq, seq_len, kmer_idx);

        // No valid k-mer remains in this sequence anymore.
        if(kmer_idx > seq_len - k)
            break;

        // Walk a maximal valid contiguous subsequence, and advance to the index following it.
        kmer_idx = walk_first_unitig(seq, seq_len, kmer_idx);
        if(kmer_idx == std::numeric_limits<size_t>::max())
        {
            result = false;
            return;
        }
    }

    result = true;
}


template <uint16_t k>
size_t Validator<k>::walk_first_unitig(const char* const seq, const size_t seq_len, const size_t start_idx) const
{
    const Kmer<k> kmer(seq, start_idx);
    const uint64_t kmer_hash = mph->lookup(kmer.canonical());
    
    if(unitig_dir[kmer_hash] == Unitig_Dir::none)
    {
        console->error("Encountered k-mer(s) in sequence that are not flanking k-mers of any of the result unitigs,"
                        " yet unitig traversals were attempted from those. Aborting.\n");
        return std::numeric_limits<size_t>::max();
    }


    const size_t uid = unitig_id[kmer_hash];
    const std::string& unitig = U[uid];
    const Unitig_Dir dir = unitig_dir[kmer_hash];

    if(!walk_unitig(seq, seq_len, start_idx, unitig, dir))
    {
        console->error("Mismatching base(s) found during walking a resultant unitig. Aborting.\n");
        return std::numeric_limits<size_t>::max();
    }

    return start_idx + unitig.length() - k + 1;
}


template <uint16_t k>
bool Validator<k>::walk_unitig(const char* const seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const Unitig_Dir dir) const
{
    if(dir == Unitig_Dir::either)
        return walk_unitig(seq, seq_len, start_idx, unitig, true) || walk_unitig(seq, seq_len, start_idx, unitig, false);

    return walk_unitig(seq, seq_len, start_idx, unitig, dir == Unitig_Dir::fwd);
}


template <uint16_t k>
bool Validator<k>::walk_unitig(const char* const seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const bool in_forward) const
{
    if(in_forward)
    {
        for(size_t idx = 0; idx < unitig.length(); ++idx)
            if(start_idx + idx >= seq_len || DNA_Utility::upper(seq[start_idx + idx]) != unitig[idx])
                return false;

        return true;
    }


    const size_t len = unitig.length();
    for(size_t idx = 0; idx < len; ++idx)
        if(start_idx + idx >= seq_len || DNA_Utility::upper(seq[start_idx + idx]) != DNA_Utility::complement(unitig[len - 1 - idx]))
            return false;
            
    return true;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Validator)
