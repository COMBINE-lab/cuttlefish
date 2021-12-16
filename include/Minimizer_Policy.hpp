
#ifndef MINIMIZER_POLICY_HPP
#define MINIMIZER_POLICY_HPP



#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>


// Forward declarations.
template <uint16_t k> class Kmer_SPMC_Iterator;
class Spin_Lock;


// A class to manage l-minimizer related policies and functions for k-mers.
template <uint16_t k, uint8_t l>
class Minimizer_Policy
{
private:

    // Maximum supported length for minimizers. Note that, memory usage for the associated data structures
    // may grow exponentially, as `4 ^ l` different minimizers are possible.
    static constexpr uint8_t MAX_LEN = 16;
    
    static constexpr uint32_t NUM_LMERS = 0b1U << (2 * l);   // Number of different possible `l`-mers.
    std::string kmer_db_path;   // Path to the underlying k-mer database.
    std::vector<uint32_t> order;    // `order[i]` denotes the order of the minimizer `i` in the policy.


    // Sets the lexicographic ordering for the l-minimizers of the k-mers.
    void set_lexicographic_ordering();

    // Sets a random ordering for the l-minimizers of the k-mers.
    void set_random_ordering();

    // Sets the frequency-based ordering for the l-minimizers of the k-mers, using up-to `thread_count`
    // number of processor threads.
    void set_frequency_ordering(uint16_t thread_count);

    // Counts the l-mers provided to the consumer thread with ID `thread_id` by the k-mer parser `parser`.
    // The count results are stored into the vector `count` — `count[i]` is the frequency of the l-mer `i`.
    // The spin-lock `lock` is used for thread-safe access to `count`.
    void count_lmers(Kmer_SPMC_Iterator<k>& parser, uint16_t thread_id, std::vector<uint64_t>& count, Spin_Lock& lock);

    // Counts the l-mer minimizers of the k-mers provided to the consumer thread with ID `thread_id` by the
    // k-mer parser `parser`. The count results are stored into the vector `count` — `count[i]` is the
    // frequency of the minimizer `i`. The spin-lock `lock` is used for thread-safe access to count.
    void count_minimizers(Kmer_SPMC_Iterator<k>& parser, uint16_t thread_id, std::vector<uint64_t>& count, Spin_Lock& lock);

public:

    // Minimizer-ordering policies supported.
    enum class Policy: uint8_t
    {
        lexicographic,
        random,
        frequency,
    };

    // Constructs an l-minimizer policy object for the k-mers in the database with path `kmer_db_path`,
    // ordering policy set as per `policy`. Construction may use up-to `thread_count` number of processor
    // threads depending upon `policy`.
    Minimizer_Policy(const std::string& kmer_db_path, Policy policy = Policy::lexicographic, uint16_t thread_count = 1);

    // Prints some statistics over the l-minimizers of the underlying k-mer set.
    void print_minimizer_stats(uint16_t thread_count);
};



#endif
