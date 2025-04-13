
#ifndef INDEX_VALIDATOR_HPP
#define INDEX_VALIDATOR_HPP




#include "Kmer_Index.hpp"
#include "Minimizer_Utility.hpp"
#include "Ref_Parser.hpp"

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>


// =============================================================================

// TODO: better document and better log messages.


// A class to contain the validation algorithms for the Cuttlefish indexings.
template <uint16_t k, uint16_t l>
class Index_Validator
{
public:

    Index_Validator() = delete;

    // Validates the indexing algorithm by constructing an index naively for the
    // sequences at the file `seq_path` and validating the index constructed by
    // the algorithm against the naive index. Indexing is over `l`-minimizers.
    static bool validate(const std::string& file_path);

    // Validates the k-mer index stored at path `idx_path`, which is supposed to
    // be over the sequences stored at path `seq_path`. Indexing is over `l`-
    // minimizers.
    static bool validate(const std::string& seq_path, const std::string& idx_path);

    // Validates the k-mer index stored at path `idx_path`, which is supposed to
    // be over the sequences stored at path `seq_path`. Indexing is over `l`-
    // minimizers.
    static bool validate(const std::string& seq_path, const std::string& idx_path, uint16_t kmer_len, uint16_t min_len);
};


template <uint16_t k, uint16_t l>
inline bool Index_Validator<k, l>::validate(const std::string& file_path)
{
    // Build a naive and a cuttlefish index.
    Ref_Parser parser(file_path);

    std::string paths;  // The concatenated paths.
    std::vector<std::size_t> ends;  // Endpoints of the paths in the concatenated sequence.

    Kmer_Index<k> kmer_index(l, 1, true);   // The cuttlefish index.
    const auto token = kmer_index.get_token();

    std::unordered_map<cuttlefish::minimizer_t, std::vector<std::size_t>> M;    // The minimizer table.
    std::size_t inst_count = 0;

    while(parser.read_next_seq())
    {
        if(parser.seq_len() < k)
            continue;

        const char* const seq = parser.seq();
        const std::size_t len = parser.seq_len();
        std::size_t last_min_idx = len;

        kmer_index.deposit(token, seq, len);

        // Collect the minimizer instances.
        for(std::size_t kmer_idx = 0; kmer_idx + k <= len; ++kmer_idx)
        {
            Kmer<l> min_lmer(seq, kmer_idx);
            std::size_t min_idx = kmer_idx;
            uint64_t min_hash = Minimizer_Utility::hash(min_lmer.as_int());

            for(std::size_t i = kmer_idx + 1; i + l <= kmer_idx + k; ++i)
            {
                const Kmer<l> lmer(seq, i);
                const uint64_t lmer_hash = Minimizer_Utility::hash(lmer.as_int());

                if(min_hash < lmer_hash)
                    continue;

                if(min_hash > lmer_hash || min_lmer > lmer)
                    min_lmer = lmer, min_idx = i, min_hash = lmer_hash;
            }

            if(min_idx != last_min_idx)
            {
                const std::size_t abs_min_idx = paths.size() + min_idx;
                M[min_lmer.as_int()].emplace_back(abs_min_idx);
                inst_count++;
                last_min_idx = min_idx;
            }
        }

        paths += seq;
        ends.emplace_back(paths.size());
    }


    kmer_index.index();

    std::cout << "Constructed the naive and the Cuttlefish index.\n";

    std::cout <<    "\n\nCross-checking the indices.\n"
                    "===============================\n\n";

    // Check paths' validity.

    std::cout << "Path counts:\n\tNaive idx: " << ends.size() << ", Cuttlefish idx: " << kmer_index.path_count_ << ".\n";
    if(ends.size() != kmer_index.path_count_)
        return false;

    for(std::size_t i = 0; i < ends.size(); ++i)
        if(ends[i] != kmer_index.path_ends[i])
            return false;

    std::cout << "Path endpoint indices matched.\n";

    std::cout << "Path sequence lengths:\n\tNaive idx: " << paths.size() << ", Cuttlefish idx: " << kmer_index.sum_paths_len << ".\n";
    if(paths.size() != kmer_index.sum_paths_len)
        return false;

    for(std::size_t i = 0; i < paths.size(); ++i)
        if(paths[i] != DNA_Utility::map_char(DNA::Base(kmer_index.paths[i].get())))
            return false;

    std::cout << "Path sequences aligned.\n";


    // Check minimizer instance counts and their offsets in the paths.

    std::cout << "Unique minimizer count:\n\tNaive idx: " << M.size() << ", Cuttlefish idx: " << kmer_index.min_count_ << ".\n";
    std::cout << "Minimizer instance count:\n\tNaive idx: " << inst_count << ", Cuttlefish idx: " << kmer_index.num_instances_ << ".\n";
    if(M.size() != kmer_index.min_count_ || inst_count != kmer_index.num_instances_)
        return false;

    const auto& mi_count = kmer_index.min_inst_count;
    const auto& m_offset = *kmer_index.min_offset;
    for(const auto p : M)
    {
        const cuttlefish::minimizer_t min = p.first;
        const uint64_t hash = kmer_index.hash(min) - 1;
        if(hash >= M.size())
        {
            std::cout << "Alien minimizer encountered.\n";
            return false;
        }

        const std::size_t count = mi_count[hash] - (hash > 0 ? mi_count[hash - 1] : 0);
        if(count != p.second.size())
        {
            std::cout << "Instance count for minimizer: " << Kmer<l>(p.first) << " with hash " << hash << ":\n\t"
                            "Naive idx: " << p.second.size() << ", Cuttlefish idx: " << count << "\n";
            std::cout << "Instance counts don't match for some minimizers.\n";
            return false;
        }


        const std::size_t min_blk_off = (hash > 0 ? mi_count[hash - 1] : 0);
        for(std::size_t i = 0; i < count; ++i)
            if(m_offset[min_blk_off + i] != p.second[i])
            {
                std::cout << "Differing instance offset for minimizer: " << Kmer<l>(p.first) << " with hash " << hash << ":\n\t"
                                "Naive idx: " << p.second[i] << ", Cuttlefish idx: " << m_offset[min_blk_off + i] << "\n";
                std::cout << "Instance offsets don't match for some minimizers.\n";
                return false;
            }
    }

    std::cout << "Instance count and offsets of individual minimizers matched.\n";


    return true;
}


template <uint16_t k, uint16_t l>
inline bool Index_Validator<k, l>::validate(const std::string& seq_path, const std::string& idx_path)
{
    // Load the index.
    Kmer_Index<k> kmer_idx(idx_path);

    std::cout << "Loaded the Cuttlefish index.\n";

    if(kmer_idx.l() != l)
    {
        std::cout <<    "The minimizer length in the k-mer index is " << kmer_idx.l() <<
                        ", while the validation is requested for minimizer length " << l << ".\n";
        return false;
    }

    const std::size_t path_count = kmer_idx.path_count_;
    const auto& paths = kmer_idx.paths;
    const auto& p_end = kmer_idx.path_ends;


    // Load the original sequences.
    std::vector<std::string> seqs_original;
    seqs_original.reserve(path_count);

    Ref_Parser parser(seq_path);
    while(parser.read_next_seq())
    {
        if(parser.seq_len() < k)
            continue;

        seqs_original.emplace_back(parser.seq());
    }

    parser.close();

    std::cout << "Loaded the original paths.\n";
    if(seqs_original.size() != path_count)
    {
        std::cout << "Path counts:\ntSequence file: " << seqs_original.size() << ", Cuttlefish idx: " << path_count << ".\n";
        return false;
    }


    // Load the sequences from the index and collect the minimizer instances.
    std::vector<std::string> seqs_idx;
    seqs_idx.reserve(path_count);

    std::unordered_map<cuttlefish::minimizer_t, std::vector<std::size_t>> M;    // The minimizer table.
    std::size_t inst_count = 0;
    std::size_t sum_paths_len = 0;
    std::size_t last_idx = 0;
    std::string path;
    std::size_t kmer_id = 0;
    typename Kmer_Index<k>::Kmer_Alignment result;

    for(std::size_t path_id = 0; path_id < path_count; ++path_id)
    {
        path.clear();
        for(std::size_t i = last_idx; i < p_end[path_id]; ++i)
            path += DNA_Utility::map_char(DNA::Base(paths[i]));

        last_idx = p_end[path_id];
        seqs_idx.emplace_back(path);


        const char* const seq = path.c_str();
        const std::size_t len = path.length();

        std::size_t last_min_idx = len;
        for(std::size_t idx = 0; idx + k <= len; ++idx)
        {
            Kmer<l> min_lmer(seq, idx);
            std::size_t min_idx = idx;
            uint64_t min_hash = Minimizer_Utility::hash(min_lmer.as_int());

            for(std::size_t i = idx + 1; i + l <= idx + k; ++i)
            {
                const Kmer<l> lmer(seq, i);
                const uint64_t lmer_hash = Minimizer_Utility::hash(lmer.as_int());

                if(min_hash < lmer_hash)
                    continue;

                if(min_hash > lmer_hash || min_lmer > lmer)
                    min_lmer = lmer, min_idx = i, min_hash = lmer_hash;
            }

            if(min_idx != last_min_idx)
            {
                const std::size_t abs_min_idx = sum_paths_len + min_idx;
                M[min_lmer.as_int()].emplace_back(abs_min_idx);
                inst_count++;
                last_min_idx = min_idx;
            }

            // Align each k-mer to the paths in the index.
            const Kmer<k> kmer(seq, idx);
            if(!kmer_idx.align(kmer, sum_paths_len + idx))
            {
                std::cout << "Non-aligning true-positive k-mer: " << kmer << "\n";

                std::cout << "Some true-positive k-mer don't align to the index.\n";
                return false;
            }

            // Check if the k-mer's containing path ID is correct.
            if(!kmer_idx.query(kmer, result) || result.path_id() != path_id || result.kmer_id() != kmer_id || result.kmer_id_in_path() != idx)
            {
                std::cout << "Query failed for k-mer: " << kmer << "\n";

                std::cout << "Some k-mer queries failed for true-positive k-mers.\n";
                return false;
            }


            kmer_id++;
        }

        sum_paths_len += len;
    }

    std::cout << "Loaded the paths from the index and constructed the naive index.\n";
    std::cout << "All k-mers in the index aligned to the index itself.\n";


    // Check if the paths align exactly.
    std::sort(seqs_original.begin(), seqs_original.end());
    std::sort(seqs_idx.begin(), seqs_idx.end());
    for(std::size_t seq_id = 0; seq_id < path_count; ++seq_id)
        if(seqs_original[seq_id] != seqs_idx[seq_id])
        {
            std::cout << "Path sequences don't match for some paths.\n";
            return false;
        }

    std::cout << "Path sequences aligned.\n";


    std::cout << "Unique minimizer count:\n\tNaive idx: " << M.size() << ", Cuttlefish idx: " << kmer_idx.min_count_ << ".\n";
    std::cout << "Minimizer instance count:\n\tNaive idx: " << inst_count << ", Cuttlefish idx: " << kmer_idx.num_instances_ << ".\n";
    if(M.size() != kmer_idx.min_count_ || inst_count != kmer_idx.num_instances_)
        return false;

    const auto& mi_count = kmer_idx.min_inst_count;
    const auto& m_offset = *kmer_idx.min_offset;

    std::vector<std::size_t> offs;
    for(auto& p : M)
    {
        const cuttlefish::minimizer_t min = p.first;
        const uint64_t hash = kmer_idx.hash(min) - 1;
        if(hash >= M.size())
        {
            std::cout << "Alien minimizer encountered.\n";
            return false;
        }

        const std::size_t count = mi_count[hash] - (hash > 0 ? mi_count[hash - 1] : 0);
        if(count != p.second.size())
        {
            std::cout << "Instance count for minimizer: " << Kmer<l>(p.first) << " with hash " << hash << ":\n\t"
                            "Naive idx: " << p.second.size() << ", Cuttlefish idx: " << count << "\n";
            std::cout << "Instance counts don't match for some minimizers.\n";
            return false;
        }

        offs.clear();
        const std::size_t min_blk_off = (hash > 0 ? mi_count[hash - 1] : 0);
        for(std::size_t i = 0; i < count; ++i)
            offs.emplace_back(m_offset[min_blk_off + i]);

        std::sort(p.second.begin(), p.second.end());
        std::sort(offs.begin(), offs.end());

        if(offs != p.second)
        {
            std::cout << "Differing instance offsets for minimizer: " << Kmer<l>(min) << " with hash " << hash << ":\n\t";

            std::cout << "Naive index:";
            std::for_each(p.second.begin(), p.second.end(), [](const auto x) { std::cout << " " << x; } );

            std::cout << "\n\tCuttlefish index:";
            std::for_each(offs.begin(), offs.end(), [](const auto x) { std::cout << " " << x; } );

            std::cout << "Instance offsets don't match for some minimizers.\n";
            return false;
        }
    }

    std::cout << "Instance count and offsets of individual minimizers matched.\n";


    return true;
}


template <uint16_t k, uint16_t l>
bool Index_Validator<k, l>::validate(const std::string& seq_path, const std::string& idx_path, const uint16_t kmer_len, const uint16_t min_len)
{
    if constexpr(l == 0 || k == 1)
        return false;
    else
    {
        if(l < min_len)
        {
            std::cerr << "The provided minimizer length is not valid. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        if(l > min_len)
            return Index_Validator<k, l - 1>::validate(seq_path, idx_path, kmer_len, min_len);


        if(k < kmer_len)
        {
            std::cerr << "The provided k-mer length is not valid. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        if(k > kmer_len)
            return Index_Validator<k - 2, l>::validate(seq_path, idx_path, kmer_len, min_len);


        return validate(seq_path, idx_path);
    }
}



#endif
