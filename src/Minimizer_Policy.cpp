
#include "Minimizer_Policy.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Spin_Lock.hpp"
#include "globals.hpp"

#include <memory>
#include <utility>
#include <numeric>
#include <random>
#include <thread>
#include <functional>
#include <algorithm>


template <uint16_t k, uint8_t l>
Minimizer_Policy<k, l>::Minimizer_Policy(const std::string& kmer_db_path, const Policy policy, const uint16_t thread_count):
    kmer_db_path(kmer_db_path),
    order(NUM_LMERS)
{
    switch(policy)
    {
        default:
        case Policy::lexicographic:
            set_lexicographic_ordering();
            break;

        case Policy::random:
            set_random_ordering();
            break;
        
        case Policy::frequency:
            set_frequency_ordering(thread_count);
            break;
    }
}


template <uint16_t k, uint8_t l>
void Minimizer_Policy<k, l>::set_lexicographic_ordering()
{
    std::iota(order.begin(), order.end(), 0U);
}


template <uint16_t k, uint8_t l>
void Minimizer_Policy<k, l>::set_random_ordering()
{
    set_lexicographic_ordering();
    
    std::shuffle(order.begin(), order.end(), std::mt19937(std::random_device()()));
}


template <uint16_t k, uint8_t l>
void Minimizer_Policy<k, l>::set_frequency_ordering(const uint16_t thread_count)
{
    const Kmer_Container<k> kmer_container(kmer_db_path);
    Kmer_SPMC_Iterator<k> parser(&kmer_container, thread_count);


    parser.launch_production();

    std::vector<uint64_t> count(NUM_LMERS);
    Spin_Lock lock;
    std::vector<std::unique_ptr<std::thread>> T(thread_count);

    for(uint16_t thread_id = 0; thread_id < thread_count; ++thread_id)
        T[thread_id].reset(
            new std::thread(&Minimizer_Policy::count_lmers, this, std::ref(parser), thread_id, std::ref(count), std::ref(lock))
        );

    parser.seize_production();

    for(uint16_t thread_id = 0; thread_id < thread_count; ++thread_id)
        T[thread_id]->join();
    

    const std::size_t max_freq_lmer = std::max_element(count.begin(), count.end()) - count.begin();
    std::cout << "Most frequent l-mer:   " << max_freq_lmer << ".\n";
    std::cout << "Associated frequency:  " << count[max_freq_lmer] << ".\n";

    const std::size_t min_freq_lmer = std::min_element(count.begin(), count.end()) - count.begin();
    std::cout << "Least frequent l-mer:  " << min_freq_lmer << ".\n";
    std::cout << "Associated frequency:  " << count[min_freq_lmer] << ".\n";


    std::vector<std::pair<uint64_t, uint32_t>> freq_lmer_pair;
    freq_lmer_pair.reserve(NUM_LMERS);

    for(uint32_t lmer = 0; lmer < NUM_LMERS; ++lmer)
        freq_lmer_pair.emplace_back(count[lmer], lmer);

    std::sort(freq_lmer_pair.begin(), freq_lmer_pair.end());

    
    for(uint32_t idx = 0; idx < NUM_LMERS; ++idx)
        order[freq_lmer_pair[idx].second] = idx;
}


template <uint16_t k, uint8_t l>
void Minimizer_Policy<k, l>::count_lmers(Kmer_SPMC_Iterator<k>& parser, const uint16_t thread_id, std::vector<uint64_t>& count, Spin_Lock& lock)
{
    std::vector<uint64_t> local_count(NUM_LMERS);
    Kmer<k> kmer;

    while(parser.tasks_expected(thread_id))
        if(parser.value_at(thread_id, kmer))
            kmer.template count_lmers<l>(local_count);

    lock.lock();
    std::transform(count.begin(), count.end(), local_count.begin(), count.begin(), std::plus<uint64_t>());
    lock.unlock();
}


template <uint16_t k, uint8_t l>
void Minimizer_Policy<k, l>::print_minimizer_stats(const uint16_t thread_count)
{
    const Kmer_Container<k> kmer_container(kmer_db_path);
    Kmer_SPMC_Iterator<k> parser(&kmer_container, thread_count);

    parser.launch_production();

    std::vector<uint64_t> count(NUM_LMERS);
    Spin_Lock lock;
    std::vector<std::unique_ptr<std::thread>> T(thread_count);

    for(uint16_t thread_id = 0; thread_id < thread_count; ++thread_id)
        T[thread_id].reset(
                new std::thread(&Minimizer_Policy<k, l>::count_minimizers, this, std::ref(parser), thread_id, std::ref(count), std::ref(lock))
        );

    parser.seize_production();

    for(uint16_t thread_id = 0; thread_id < thread_count; ++thread_id)
        T[thread_id]->join();

    
    const std::size_t max_freq_minmzr = std::max_element(count.begin(), count.end()) - count.begin();
    std::cout << "Most frequent l-minmizer:     " << max_freq_minmzr << ".\n";
    std::cout << "Associated frequency:         " << count[max_freq_minmzr] << ".\n";

    const std::size_t min_freq_minmzr = std::min_element(count.begin(), count.end()) - count.begin();
    std::cout << "Least frequent l-minimizer:   " << min_freq_minmzr << ".\n";
    std::cout << "Associated frequency:         " << count[min_freq_minmzr] << ".\n";
}


template <uint16_t k, uint8_t l>
void Minimizer_Policy<k, l>::count_minimizers(Kmer_SPMC_Iterator<k>& parser, const uint16_t thread_id, std::vector<uint64_t>& count, Spin_Lock& lock)
{
    std::vector<uint64_t> local_count(NUM_LMERS);
    Kmer<k> kmer;

    while(parser.tasks_expected(thread_id))
        if(parser.value_at(thread_id, kmer))
            local_count[kmer.template minimizer<l>(order)]++;

    lock.lock();
    std::transform(count.begin(), count.end(), local_count.begin(), count.begin(), std::plus<uint64_t>());
    lock.unlock();
}
