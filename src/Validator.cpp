
#include "Validator.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"
#include "Kmer_Hasher.hpp"
#include "BBHash/BooPHF.h"
#include "Directed_Kmer.hpp"
#include "kseq/kseq.h"

#include "spdlog/sinks/stdout_color_sinks.h"

#include <fstream>
#include <sys/stat.h>
#include <thread>


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


Validator::Validator(const Validation_Params& params, cuttlefish::logger_t console):
    params(params), k(params.k()), console(console)
{
    cuttlefish::kmer_t::set_k(k);
}


bool Validator::validate()
{
    build_mph_function();


    bool is_valid_sequence;
    std::thread sequence_validator(&Validator::validate_sequence_completion, this, std::ref(is_valid_sequence));

    bool is_valid_kmer_set;
    std::thread kmer_set_validator(&Validator::validate_kmer_set, this, std::ref(is_valid_kmer_set));


    if(!sequence_validator.joinable())
    {
        console->error("Early termination faced for a worker thread. Aborting.\n");
        std::exit(EXIT_FAILURE);
    }

    sequence_validator.join();
    console->info("{} validation of complete coverage of the sequence by the produced unitigs.\n",
                    is_valid_sequence ? "Passed" : "Failed");

    
    if(!kmer_set_validator.joinable())
    {
        console->error("Early termination faced for a worker thread. Aborting.\n");
        std::exit(EXIT_FAILURE);
    }

    kmer_set_validator.join();
    console->info("{} validation of the k-mer set.\n", is_valid_kmer_set ? "Passed" : "Failed");


    return is_valid_kmer_set && is_valid_sequence;
}


void Validator::build_mph_function()
{
    const std::string& kmc_db_path = params.kmc_db_path();
    const uint16_t thread_count = params.thread_count();
    const std::string& mph_file_path = params.mph_file_path();

    const Kmer_Container kmer_container(kmc_db_path);

    // The serialized BBHash file (saved from some earlier execution) exists.
    struct stat buffer;
    if(stat(mph_file_path.c_str(), &buffer) == 0)
    {
        console->info("Loading the MPH function from file {}\n", mph_file_path);
        
        std::ifstream input(mph_file_path.c_str(), std::ifstream::in);
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher>();
        mph->load(input);
        input.close();
        
        console->info("Loaded the MPH function into memory.\n");
    }
    else    // No BBHash file exists. Build and save one now.
    {
        // Build the MPHF.
        console->info("Building the MPH function from the k-mer database {}\n", kmer_container.container_location());

        auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> (kmer_container.size(), data_iterator, thread_count, gamma_factor);

        console->info("Built the MPH function in memory.\n");
        

        // Save the MPHF.
        console->info("Saving the MPH function in file {}\n", mph_file_path);

        std::ofstream output(mph_file_path.c_str(), std::ofstream::out);
        mph->save(output);
        output.close();

        console->info("Saved the MPH function in disk.\n");
    }
}


void Validator::validate_kmer_set(bool& result) const
{
    console->info("Testing validation of the uniqueness of the k-mers and completeness of the k-mer set in the produced unitigs.\n");

    const std::string& kmc_db_path = params.kmc_db_path();
    const std::string& cdbg_file_path = params.cdbg_file_path();

    const Kmer_Container kmer_container(kmc_db_path);
    const uint64_t kmer_count = kmer_container.size();

    console->info("Number of k-mers in the database: {}\n", kmer_count);

    std::vector<bool> is_present(kmer_count);
    uint64_t kmers_seen = 0;
    uint64_t kmers_repeated = 0;
    uint64_t unitigs_processed = 0;
    uint64_t kmers_invalid = 0;

    // Scan through the unitigs one-by-one.
    std::string unitig;
    std::ifstream input(cdbg_file_path.c_str(), std::ifstream::in);
    while(input >> unitig)
    {
        cuttlefish::kmer_t first_kmer(unitig, 0);
        Directed_Kmer kmer(first_kmer);

        // Scan through the k-mers one-by-one.
        for(size_t kmer_idx = 0; kmer_idx <= unitig.length() - k; ++kmer_idx)
        {
            uint64_t hash_val = mph->lookup(kmer.canonical());

            // Encountered a k-mer that is absent at the k-mer database and hashes outside of the valid range.
            if(hash_val >= kmer_count)
            {
                console->error("Invalid k-mer encountered.\n");
                kmers_invalid++;
            }
            // Encountered a k-mer for the first time that is either a k-mer present at the database, or is
            // absent there but hashes to the value of a present one (and that present one hasn't been seen yet).
            else if(!is_present[hash_val])
                is_present[hash_val] = true;
            // A repeated k-mer is seen.
            else
            {
                console->info("Repeated k-mer encountered.\n");
                kmers_repeated++;
            }

            if(kmer_idx < unitig.length() - k)
                kmer.roll_to_next_kmer(unitig[kmer_idx + k]);
        }


        kmers_seen += unitig.length() - k + 1;
        unitigs_processed++;

        if(unitigs_processed % PROGRESS_GRAIN_SIZE == 0)
            console->info("Validated {}M unitigs.\n", unitigs_processed / 1000000);
    }


    console->info("Total number of repeated k-mers: {}\n", kmers_repeated);
    console->info("Total number of invalid k-mers: {}\n", kmers_invalid);
    console->info("Total number of k-mers seen: {}\n", kmers_seen);
    console->info("Total number of k-mers expected: {}\n", kmer_count);

    input.close();

    result = (!kmers_repeated && !kmers_invalid && kmers_seen == kmer_count);
}


void Validator::validate_sequence_completion(bool& result)
{
    console->info("Testing validation of the completeness of coverage of the sequence by the produced unitigs.\n");

    const std::string& ref_file_path = params.ref_file_path();
    const std::string& kmc_db_path = params.kmc_db_path();
    const uint16_t thread_count = params.thread_count();

    const Kmer_Container kmer_container(kmc_db_path);
    const uint64_t kmer_count = kmer_container.size();

    console->info("Number of k-mers in the k-mer database: {}\n", kmer_count);


    // Allocate the unitig tables.
    console->info("Allocating the unitig tables.\n");

    unitig_id.resize(kmer_count);
    unitig_dir.resize(kmer_count);

    console->info("Done allocation of the unitig tables.\n");


    // Load the unitigs into memory and build the associated tables `unitig_id` and `unitig_dir`.
    build_unitig_tables();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* const input = fopen(ref_file_path.c_str(), "r");
    if(input == NULL)
    {
        console->error("Error opening the reference file {}. Aborting.\n", ref_file_path);
        std::exit(EXIT_FAILURE);
    }
    
    // Initialize the parser.
    kseq_t* const parser = kseq_init(fileno(input));


    std::vector<std::thread> th(thread_count);  // Thread-pool (round-robin) to validate the sequences parallelly.
    uint32_t seqCount = 0;  // Number of sequences read.
    bool success = true;    // Whether the total validation succeeded.
    char** S = new char*[thread_count]; // Sequence buffers to pass on to the threads.
    bool* thread_result = new bool[thread_count];   // Validation result produced by the threads.
    
    // Parse sequences one-by-one, and continue spelling them using the resultant unitigs of the compaction algorithm.
    while(kseq_read(parser) >= 0)
    {
        const char* const seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

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


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);

    result = success;
}


void Validator::build_unitig_tables()
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
            const cuttlefish::kmer_t kmer(unitig);
            const uint64_t kmer_hash = mph->lookup(kmer.canonical());
            
            unitig_id[kmer_hash] = unitig_count;
            unitig_dir[kmer_hash] = Unitig_Dir::either;
        }
        else
        {
            const cuttlefish::kmer_t last_kmer = cuttlefish::kmer_t(unitig.substr(unitig.length() - k, k));
            const cuttlefish::kmer_t kmer_fwd(unitig, 0);
            const cuttlefish::kmer_t kmer_bwd = last_kmer.reverse_complement();

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


void Validator::walk_sequence(const char* const seq, const size_t seq_len, bool& result) const
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


size_t Validator::search_valid_kmer(const char* const seq, const size_t seq_len, const size_t start_idx) const
{
    size_t valid_start_idx;
    uint16_t nucl_count;

    size_t idx = start_idx;
    while(idx <= seq_len - k)
    {
        // Go over the contiguous subsequence of 'N's.
        for(; idx <= seq_len - k && seq[idx] == cuttlefish::PLACEHOLDER_NUCLEOTIDE; idx++);

        // Go over the contiguous subsequence of non-'N's.
        if(idx <= seq_len - k)
        {
            valid_start_idx = idx;
            nucl_count = 0;

            for(; idx < seq_len && seq[idx] != cuttlefish::PLACEHOLDER_NUCLEOTIDE; ++idx)
                if(++nucl_count == k)
                    return valid_start_idx;
        }
    }


    return seq_len;
}


size_t Validator::walk_first_unitig(const char* const seq, const size_t seq_len, const size_t start_idx) const
{
    const cuttlefish::kmer_t kmer(seq, start_idx);
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
        console->error("Mismatching nucleotide(s) found during walking a resultant unitig. Aborting.\n");
        return std::numeric_limits<size_t>::max();
    }

    return start_idx + unitig.length() - k + 1;
}


bool Validator::walk_unitig(const char* const seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const Unitig_Dir dir) const
{
    if(dir == Unitig_Dir::either)
        return walk_unitig(seq, seq_len, start_idx, unitig, true) || walk_unitig(seq, seq_len, start_idx, unitig, false);

    return walk_unitig(seq, seq_len, start_idx, unitig, dir == Unitig_Dir::fwd);
}


bool Validator::walk_unitig(const char* const seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const bool in_forward) const
{
    if(in_forward)
    {
        for(size_t idx = 0; idx < unitig.length(); ++idx)
            if(start_idx + idx >= seq_len || seq[start_idx + idx] != unitig[idx])
                return false;

        return true;
    }


    const size_t len = unitig.length();
    for(size_t idx = 0; idx < len; ++idx)
        if(start_idx + idx >= seq_len || seq[start_idx + idx] != cuttlefish::kmer_t::complement(unitig[len - 1 - idx]))
            return false;
            
    return true;
}
