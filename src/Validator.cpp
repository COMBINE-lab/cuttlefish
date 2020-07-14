
#include "Validator.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"
#include "Kmer_Hasher.hpp"
#include "BBHash/BooPHF.h"
#include "Directed_Kmer.hpp"
#include "kseq/kseq.h"

#include <fstream>
#include <sys/stat.h>


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


Validator::Validator(const std::string& ref_file_name, const uint16_t k, const std::string& kmc_db_name, const std::string& cdbg_file_name):
    ref_file_name(ref_file_name), k(k), kmc_db_name(kmc_db_name), cdbg_file_name(cdbg_file_name), mph(NULL)
{
    Kmer::set_k(k);
}


bool Validator::validate(const std::string& bbhash_file_name, const uint16_t thread_count)
{
    build_mph_function(bbhash_file_name, thread_count);

    const bool valid_kmer_set = validate_kmer_set();
    std::cout << (valid_kmer_set ? "Passed" : "Failed") << " validation of the k-mer set.\n";

    const bool valid_sequence = validate_sequence_completion();
    std::cout << (valid_sequence ? "Passed" : "Failed") << " validation of complete coverage of the sequence"
                                                           " by the produced unitigs.\n";

    return valid_kmer_set && valid_sequence;
}


void Validator::build_mph_function(const std::string& bbhash_file_name, const uint16_t thread_count)
{
    const Kmer_Container kmer_container(kmc_db_name);

    // The serialized BBHash file (saved from some earlier execution) exists.
    struct stat buffer;
    if(stat(bbhash_file_name.c_str(), &buffer) == 0)
    {
        std::cout << "Loading the MPH function from file " << bbhash_file_name << "\n";
        
        std::ifstream input(bbhash_file_name.c_str(), std::ifstream::in);
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher>();
        mph->load(input);
        input.close();
        
        std::cout << "Loaded the MPH function into memory.\n";
    }
    else    // No BBHash file exists. Build and save one now.
    {
        // Build the MPHF.
        std::cout << "Building the MPH function from the k-mer database " << kmer_container.container_location() << "\n";

        auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> (kmer_container.size(), data_iterator, thread_count, gamma_factor);

        std::cout << "Built the MPH function in memory.\n";
        

        // Save the MPHF.
        std::cout << "Saving the MPH function in file " << bbhash_file_name << "\n";

        std::ofstream output(bbhash_file_name.c_str(), std::ofstream::out);
        mph->save(output);
        output.close();

        std::cout << "Saved the MPH function in disk.\n";
    }
}


bool Validator::validate_kmer_set() const
{
    std::cout << "Testing validation of the uniqueness of the k-mers and completeness of the k-mer set in the produced unitigs.\n";

    const Kmer_Container kmer_container(kmc_db_name);
    const uint64_t kmer_count = kmer_container.size();

    std::cout << "Number of k-mers in the database: " << kmer_count << "\n";

    std::vector<bool> is_present(kmer_count);
    uint64_t kmers_seen = 0;
    uint64_t kmers_repeated = 0;
    uint64_t unitigs_processed = 0;
    uint64_t kmers_invalid = 0;

    // Scan through the unitigs one-by-one.
    std::string unitig;
    std::ifstream input(cdbg_file_name.c_str(), std::ifstream::in);
    while(input >> unitig)
    {
        cuttlefish::kmer_t first_kmer(unitig, 0);
        Directed_Kmer kmer(first_kmer);

        // Scan through the k-mers one-by-one.
        for(size_t kmer_idx = 0; kmer_idx <= unitig.length() - k; ++kmer_idx)
        {
            uint64_t hash_val = mph->lookup(kmer.canonical);

            // Encountered a k-mer that is absent at the k-mer database and hashes outside of the valid range.
            if(hash_val >= kmer_count)
            {
                std::cout << "Invalid k-mer encountered.\n";
                kmers_invalid++;
            }
            // Encountered a k-mer for the first time that is either a k-mer present at the database, or is
            // absent there but hashes to the value of a present one (and that present one hasn't been seen yet).
            else if(!is_present[hash_val])
                is_present[hash_val] = true;
            // A repeated k-mer is seen.
            else
            {
                std::cout << "Repeated k-mer encountered.\n";
                kmers_repeated++;
            }

            if(kmer_idx < unitig.length() - k)
                kmer.roll_to_next_kmer(unitig[kmer_idx + k]);
        }


        kmers_seen += unitig.length() - k + 1;
        unitigs_processed++;

        if(unitigs_processed % progress_grain_size == 0)
            std::cout << "Processed " << unitigs_processed / 1000000 << "M unitigs.\n";
    }


    std::cout << "Total number of repeated k-mers: " << kmers_repeated << "\n";
    std::cout << "Total number of invalid k-mers: " << kmers_invalid << "\n";
    std::cout << "Total number of k-mers seen: " << kmers_seen << "\n";
    std::cout << "Total number of k-mers expected: " << kmer_count << "\n";

    input.close();

    return !kmers_repeated && !kmers_invalid && kmers_seen == kmer_count;
}


bool Validator::validate_sequence_completion()
{
    std::cout << "Testing validation of the completeness of coverage of the sequence by the produced unitigs.\n";

    const Kmer_Container kmer_container(kmc_db_name);
    const uint64_t kmer_count = kmer_container.size();

    std::cout << "Number of k-mers in the k-mer database: " << kmer_count << "\n";


    // Allocate the unitig tables.
    std::cout << "Allocating the unitig tables.\n";

    unitig_id.resize(kmer_count);
    unitig_dir.resize(kmer_count);

    std::cout << "Done allocation of the unitig tables.\n";


    // Load the unitigs into memory and build the associated tables `unitig_id` and `unitig_dir`.
    build_unitig_tables();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* input = fopen(ref_file_name.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening the reference file " << ref_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
    
    // Initialize the parser.
    kseq_t* parser = kseq_init(fileno(input));


    // Parse sequences one-by-one, and continue spelling them using the resultant unitigs
    // of the compaction algorithm.
    uint32_t seqCount = 0;
    bool success = true;
    while(success && kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::cout << "Spelling out sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        // Nothing to process for sequences with length shorter than `k`.
        if(seq_len < k)
            break;

        if(!walk_sequence(seq, seq_len))
            success = false;
    }


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);

    return success;
}


void Validator::build_unitig_tables()
{
    std::ifstream input(cdbg_file_name.c_str(), std::ifstream::in);
    if(!input)
    {
        std::cerr << "Error reading from the unitigs file " << cdbg_file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::cout << "Loading the unitigs from the file " << cdbg_file_name << ".\n";

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
        if(unitig_count % progress_grain_size == 0)
            std::cout << "Loaded " << unitig_count / 1000000 << "M unitigs.\n";
    }

    std::cout << "Done loading a total of " << unitig_count << " unitigs.\n";

    input.close();
}


bool Validator::walk_sequence(const char* seq, const size_t seq_len) const
{
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
            return false;
    }

    return true;
}


size_t Validator::search_valid_kmer(const char* seq, const size_t seq_len, const size_t start_idx) const
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


size_t Validator::walk_first_unitig(const char* seq, const size_t seq_len, const size_t start_idx) const
{
    const cuttlefish::kmer_t kmer(seq, start_idx);
    const uint64_t kmer_hash = mph->lookup(kmer.canonical());
    
    if(unitig_dir[kmer_hash] == Unitig_Dir::none)
    {
        std::cout << "Encountered k-mer(s) in sequence that are not flanking k-mers of any of the result unitigs,"
                     " yet unitig traversals were attempted from those. Aborting.\n";
        return std::numeric_limits<size_t>::max();
    }


    const size_t uid = unitig_id[kmer_hash];
    const std::string& unitig = U[uid];
    const Unitig_Dir dir = unitig_dir[kmer_hash];

    if(!walk_unitig(seq, seq_len, start_idx, unitig, dir))
    {
        std::cout << "Mismatching nucleotide(s) found during walking a resultant unitig. Aborting.\n";
        return std::numeric_limits<size_t>::max();
    }

    return start_idx + unitig.length() - k + 1;
}


bool Validator::walk_unitig(const char* seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const Unitig_Dir dir) const
{
    if(dir == Unitig_Dir::either)
        return walk_unitig(seq, seq_len, start_idx, unitig, true) || walk_unitig(seq, seq_len, start_idx, unitig, false);

    return walk_unitig(seq, seq_len, start_idx, unitig, dir == Unitig_Dir::fwd);
}


bool Validator::walk_unitig(const char* seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const bool in_forward) const
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
        if(start_idx + idx >= seq_len || seq[start_idx + idx] != Kmer::complement(unitig[len - 1 - idx]))
            return false;
            
    return true;
}
