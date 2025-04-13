
#include "Kmer_Index.hpp"
#include "Sparse_Lock.hpp"
#include "Minimizer_Iterator.hpp"
#include "CdBG.hpp"
#include "Read_CdBG.hpp"
#include "Build_Params.hpp"
#include "Input_Defaults.hpp"
#include "utility.hpp"
#include "key-value-collator/Key_Value_Iterator.hpp"
#include "elias_fano/sequence.hpp"

#include <string>
#include <functional>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <thread>
#include <chrono>


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const uint16_t l, const uint16_t producer_count, const bool retain, const std::string& output_pref, const std::string& working_dir, const std::optional<Build_Params> params):
    // TODO: reserve space for `paths`, preferably from additional k-mer count field
    l_(l),
    output_pref(output_pref),
    working_dir(working_dir + "/"),
    producer_count(producer_count),
    worker_count(producer_count),
    params(params),
    retain(retain),
    num_instances_(0),
    min_count_(0),
    max_inst_count_(0),
    producer_path_buf(producer_count),
    producer_path_end_buf(producer_count),
    producer_min_inst_buf(producer_count, nullptr),
    min_collator(min_instance_path_pref(), 2 * producer_count),
    overflow_min_count_(0),
    overflow_kmer_count_(0),
    serialize_stream(index_file_path(output_pref), std::ios::out | std::ios::binary),
    curr_token(0)
{
    assert(l <= 32);

    save_config();
}


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const std::string& idx_path):
    l_(minimizer_len(index_file_path(idx_path))),
    output_pref(idx_path),
    working_dir(dirname(idx_path) + "/"),
    producer_count(0),
    worker_count(0),
    retain(true),
    num_instances_(0),
    min_count_(0),
    max_inst_count_(0),  // TODO: think if this might have repercussions, as we don't have it saved (though, inferrable).
    overflow_min_count_(0),
    overflow_kmer_count_(0)
{
    const uint16_t idx_k = kmer_len(index_file_path(idx_path));
    if(idx_k != k)
    {
        std::cerr << "The k-mer length of the requested index is " << idx_k << ", which doesn't match with this execution. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::ifstream deserialize_stream(index_file_path(idx_path), std::ios::in | std::ios::binary);

    deserialize_stream.seekg(sizeof(k) + sizeof(l_), std::ios::beg);    // Fast forward over the `k` and the `l` values.

    paths.deserialize(deserialize_stream);
    sum_paths_len_ = paths.size();

    elias_fano::essentials::load(path_ends, std::ref(deserialize_stream));
    path_count_ = path_ends.size();

    min_mphf.emplace(deserialize_stream);
    // min_count = min_mphf->nbKeys();

    elias_fano::essentials::load(min_inst_count, std::ref(deserialize_stream));
    min_count_ = min_inst_count.size();

    min_offset.emplace();
    min_offset->deserialize(deserialize_stream);
    num_instances_ = min_offset->size();

    kmer_mphf.emplace(deserialize_stream);
    // overflow_kmer_count_ = kmer_mphf->nbKeys();

    overflow_kmer_map.emplace();
    overflow_kmer_map->deserialize(deserialize_stream);
    overflow_kmer_count_ = overflow_kmer_map->size();

    deserialize_stream.close();
    if(!deserialize_stream)
    {
        std::cerr << "Error reading from the k-mer index file at " << index_file_path(idx_path) << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    min_collator.close_deposit_stream();
}


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const Build_Params& params):
    Kmer_Index(params.min_len(), params.thread_count(), true, params.output_prefix(), params.working_dir_path(), params)
{}


template <uint16_t k>
void Kmer_Index<k>::construct()
{
    assert(params);

    if(params->is_read_graph() || params->is_ref_graph())
    {
        Read_CdBG<k> cdbg(*params, this);
        cdbg.construct();
    }
    else
    {
        CdBG<k> cdbg(*params, this);
        cdbg.construct();
    }

    index();
}


template <uint16_t k>
void Kmer_Index<k>::save_config() const
{
    // TODO: maybe also output configuration to the CF json file?
    constexpr uint16_t K = k;
    serialize_stream.write(reinterpret_cast<const char*>(&K), sizeof(K));
    serialize_stream.write(reinterpret_cast<const char*>(&l_), sizeof(l_));
}


template <uint16_t k>
const typename Kmer_Index<k>::Producer_Token Kmer_Index<k>::get_token()
{
    lock.lock();
    const std::size_t token = curr_token++;
    assert(token < producer_count);
    lock.unlock();

    return Producer_Token(token);
}


template <uint16_t k>
void Kmer_Index<k>::index()
{
    typedef std::chrono::high_resolution_clock::time_point time_point_t;
    constexpr auto now = std::chrono::high_resolution_clock::now;
    constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };

    // Dump the remaining in-memory content to disk.
    close_deposit_stream();
    std::cout << "Closed the sequence deposit stream.\n";
    const time_point_t t_deposit = now();

    // Collate the minimizer-instances.
    min_collator.collate(worker_count, true);
    num_instances_ = min_collator.pair_count();
    min_count_ = min_collator.unique_key_count();
    std::cout << "Minimizer-instance count: " << num_instances_ << ".\n";
    std::cout << "Unique minimizer count:   " << min_count_ << ".\n";
    const time_point_t t_collate = now();
    std::cout << "Collated the minimizer-instances. Time taken = " << duration(t_collate - t_deposit) << " seconds.\n";

    // Construct an MPHF over the unique minimizers.
    construct_minimizer_mphf();
    const uint64_t total_bits = min_mphf->totalBitSize();
    std::cout <<    "For the minimizer-MPHF:\n"
                    "\tTotal size: " << total_bits / (8U * 1024U * 1024U) << " MB."
                    " Bits per k-mer: " << static_cast<double>(total_bits) / min_count_ << ".\n";
    const time_point_t t_mphf = now();
    std::cout << "Constructed the minimizer-MPHF. Time taken = " << duration(t_mphf - t_collate) << " seconds.\n";

    // Count the instances per minimizer.
    count_minimizer_instances();
    const time_point_t t_count = now();
    std::cout << "Counted the minimizer instances. Time taken = " << duration(t_count - t_mphf) << " seconds.\n";

    // Get the offsets of each minimizer.
    get_minimizer_offsets();
    const time_point_t t_off = now();
    std::cout << "Gathered the minimizer instances' offsets. Time taken = " << duration(t_off - t_count) << " seconds.\n";

    // Construct the minimizer-overflow index.
    construct_overflow_index();
    const time_point_t t_overflow = now();
    std::cout << "Constructed the minimizer-overflow index. Time taken = " << duration(t_overflow - t_off) << " seconds.\n";

    // Save the index and remove temporary files.
    close_output();
}


template <uint16_t k>
void Kmer_Index<k>::close_deposit_stream()
{
    // Flush the remaining buffer content, and release memory of all the path-sequence buffers of the producers.
    for(std::size_t id = 0; id < producer_count; ++id)
    {
        if(!producer_path_buf[id].empty())
            flush(id);

        force_free(producer_path_buf[id]);
        force_free(producer_path_end_buf[id]);
    }

    min_collator.close_deposit_stream();


    sum_paths_len_ = paths.size();

    paths.serialize(serialize_stream, true);

    // Release memory of the concatenated paths if required.
    if(!retain)
        force_free(paths);


    // Elias-Fano encode the path-endpoints.
    path_count_ = path_ends_vec.size();
    path_ends.encode(path_ends_vec.cbegin(), path_count_, path_ends_vec.back());
    assert(path_ends_vec.back() == sum_paths_len_);

    force_free(path_ends_vec);

    elias_fano::essentials::save(path_ends, std::ref(serialize_stream));
    if(!retain)
        force_free(path_ends);
}


template <uint16_t k>
void Kmer_Index<k>::construct_minimizer_mphf()
{
    const auto min_collation_range = boomphf::range(min_collator.begin(), min_collator.end());
    min_mphf.emplace(min_count_, min_collation_range, working_dir, worker_count, gamma);

    min_mphf->save(serialize_stream);
}


template <uint16_t k>
void Kmer_Index<k>::count_minimizer_instances()
{
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(num_instances_ + 1)));
    assert(bits_per_entry > 0);
    min_inst_count_bv.emplace(bits_per_entry, min_count_ + 1);
    min_inst_count_bv->clear_mem();

    auto min_inst_iter = min_collator.begin();
    constexpr std::size_t min_buf_sz = 5U * 1024U * 1024U / sizeof(min_inst_t); // Minimizer-instance chunk size, in elements: 5 MB total.

    std::vector<std::thread> worker;
    Sparse_Lock<Spin_Lock> idx_lock(min_count_ + 1, idx_lock_count);

    for(uint16_t t = 0; t < worker_count; ++t)
        worker.emplace_back(
            [this, &min_inst_iter, &idx_lock]()
            {
                min_inst_t* const min_chunk = allocate<min_inst_t>(min_buf_sz);
                std::size_t buf_elem_count;
                uint64_t h;
                auto& mi_count = *min_inst_count_bv;

                while((buf_elem_count = min_inst_iter.read(min_chunk, min_buf_sz)) != 0)
                    for(std::size_t idx = 0; idx < buf_elem_count; ++idx)
                    {
                        h = hash(min_chunk[idx].first);

                        idx_lock.lock(h);
                        mi_count[h] = mi_count[h] + 1;
                        idx_lock.unlock(h);
                    }

                deallocate(min_chunk);
            }
            );

    for(uint16_t t = 0; t < worker_count; ++t)
    {
        if(!worker[t].joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[t].join();
    }

    max_inst_count_ = min_collator.mode_frequency();
    std::cout << "Maximum instance count of a minimizer: " << max_inst_count_ << ".\n";


    // Transform the instance counts to cumulative counts, i.e. prefix-sums, ordered per the minimizers' hash.
    auto& mi_count = *min_inst_count_bv;
    uint64_t pref_sum = (mi_count[0] = 0);
    for(std::size_t i = 1; i <= min_count_; ++i)
        mi_count[i] = (pref_sum += mi_count[i]);
}


template <uint16_t k>
void Kmer_Index<k>::get_minimizer_offsets()
{
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(sum_paths_len_))); // The path-sequence is 0-indexed—no need for +1.
    assert(bits_per_entry > 0);
    min_offset.emplace(bits_per_entry, num_instances_);

    auto min_inst_iter = min_collator.begin();
    constexpr std::size_t min_buf_sz = 5U * 1024U * 1024U / sizeof(min_inst_t); // Minimizer instance chunk size, in elements: 5 MB total.

    std::vector<std::thread> worker;
    Sparse_Lock<Spin_Lock> idx_lock(min_count_, idx_lock_count);

    for(uint16_t t = 0; t < worker_count; ++t)
        worker.emplace_back(
            [this, &min_inst_iter, &idx_lock]()
            {
                min_inst_t* const min_chunk = allocate<min_inst_t>(min_buf_sz);
                std::size_t buf_elem_count;
                uint64_t h;
                std::size_t offset;
                auto& mi_count = *min_inst_count_bv;
                auto& m_offset = *min_offset;

                while((buf_elem_count = min_inst_iter.read(min_chunk, min_buf_sz)) != 0)
                    for(std::size_t idx = 0; idx < buf_elem_count; ++idx)
                    {
                        h = hash(min_chunk[idx].first) - 1;

                        // Tracking the next offsets to put the minimizers' offsets, within `min_instance_count` in-place,
                        // transforms it by sliding it one index to the left. That is, now, `min_instance_count[MPH(min)]`
                        // has the prefix-sum, opposed to the index `MPH(min) + 1` as earlier.
                        idx_lock.lock(h);
                        offset = mi_count[h];
                        mi_count[h] = offset + 1;
                        idx_lock.unlock(h);

                        m_offset[offset] = min_chunk[idx].second;
                    }

                deallocate(min_chunk);
            }
        );

    for(uint16_t t = 0; t < worker_count; ++t)
    {
        if(!worker[t].joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[t].join();
    }


    // After the in-place transformation (see note above), `min_instance_count` need not have that extra last entry anymore.
    min_inst_count_bv->resize(min_count_);
    min_inst_count.encode(min_inst_count_bv->cbegin(), min_count_, min_inst_count_bv->back());
    assert(min_inst_count_bv->back() == num_instances_);

    min_inst_count_bv.reset();
    elias_fano::essentials::save(min_inst_count, std::ref(serialize_stream));

    min_offset->serialize(serialize_stream);


    // Release the portion of the index still in memory if required.
    if(!retain)
    {
        min_mphf.reset();
        force_free(min_inst_count);
        min_offset.reset();
    }
}


template <uint16_t k>
void Kmer_Index<k>::construct_overflow_index()
{
    collect_overflown_kmers();

    construct_overflow_kmer_mphf();

    map_overflown_kmers();


    // Release the overflow index if required.
    if(!retain)
    {
        kmer_mphf.reset();
        overflow_kmer_map.reset();
    }
}


template <uint16_t k>
void Kmer_Index<k>::collect_overflown_kmers()
{
    std::ofstream kmer_op(overflow_kmers_path());
    std::ofstream inst_idx_op(overflow_min_insts_path());

    const std::size_t task_size = min_count_ / worker_count;
    const uint16_t partition_count = (task_size < 1 ? 1 : worker_count);
    std::vector<std::size_t> num_min(partition_count);
    std::vector<std::size_t> num_kmer(partition_count);
    std::vector<std::thread> worker;
    worker.reserve(partition_count);

    std::size_t min_id_low = 0;
    void (Kmer_Index<k>::* const task)(std::size_t, std::size_t, std::ofstream&, std::ofstream&, std::size_t&, std::size_t&) const = &Kmer_Index::collect_overflown_kmers;
    for(uint16_t t_id = 0; t_id < partition_count; ++t_id)
    {
        const std::size_t min_id_high = (t_id == partition_count - 1 ? min_count_ : min_id_low + task_size);
        worker.emplace_back(task, this,
                            min_id_low, min_id_high, std::ref(kmer_op), std::ref(inst_idx_op), std::ref(num_min[t_id]), std::ref(num_kmer[t_id]));

        min_id_low += task_size;
    }

    for(uint16_t t_id = 0; t_id < partition_count; ++t_id)
    {
        if(!worker[t_id].joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[t_id].join();
    }

    kmer_op.close();
    inst_idx_op.close();


    overflow_min_count_ = std::accumulate(num_min.begin(), num_min.end(), uint64_t(0));
    overflow_kmer_count_ = std::accumulate(num_kmer.begin(), num_kmer.end(), uint64_t(0));

    std::cout << "Minimizer overflow threshold:     " << overflow_threshold << ".\n";
    std::cout << "Number of overflowing minimizers: " << overflow_min_count_ << ".\n";
    std::cout << "Number of corresponding k-mers:   " << overflow_kmer_count_ << ".\n";
}


template <uint16_t k>
void Kmer_Index<k>::collect_overflown_kmers(const std::size_t low, const std::size_t high, std::ofstream& kmer_op, std::ofstream& inst_idx_op, std::size_t& num_min, std::size_t& num_kmer) const
{
    const auto& m_offset = *min_offset;

    std::vector<Kmer<k>> kmers; // Buffer for the k-mers from the overflowing minimizers.
    std::vector<std::size_t> inst_idx;  // Buffer for the indices of the minimizer instances into their corresponding blocks, of the overflown k-mers.
    Kmer<k> kmer;   // Each k-mer from the overflowing minimizers.

    constexpr std::size_t buf_elem_th = buf_sz_th / (sizeof(Kmer<k>) + sizeof(std::size_t));    // Maximum number of k-mers to keep in the buffer.
    static_assert(buf_elem_th > 0);

    num_min = 0;
    num_kmer = 0;
    const auto locked_dump =    [this, &kmers, &kmer_op, &inst_idx, &inst_idx_op]()
                                {
                                    lock.lock();
                                    dump(kmers, kmer_op), dump(inst_idx, inst_idx_op);
                                    lock.unlock();
                                };

    kmers.reserve(buf_elem_th), inst_idx.reserve(buf_elem_th);
    for(std::size_t id = low; id < high; ++id)  // Go over each minimizer.
    {
        const std::size_t idx_begin = (id > 0 ? min_inst_count[id - 1] : 0);    // Offset of the block of indices of the minimizer's instances.
        const std::size_t idx_end = min_inst_count[id]; // End-offset of the instance block (exclusive).
        const std::size_t inst_count = idx_end - idx_begin; // Number of instances of this minimizer.

        if(inst_count >= overflow_threshold)
        {
            for(std::size_t i = idx_begin; i < idx_end; ++i)    // Go over each instance of this minimizer.
            {
                const std::size_t min_inst_idx = m_offset[i];   // Index of the minimizer instance in the paths-sequence.

                // Extract each k-mer from the paths that have this minimizer instance as their minimizer.

                const int64_t l = path_ends.prev_leq(min_inst_idx); // ID of the path preceding the path containing this instance.
                const int64_t r = l + 1;    // ID of the path containing this instance.
                const std::size_t l_end = (l < 0 ? 0 : path_ends[l]);   // Index of the left end of the path containing this instance.
                const std::size_t r_end = path_ends[r]; // Index of the right end (exclusive) of the path containing this instance.

                const std::size_t window_start = (l_end + k - l_ <= min_inst_idx ? (min_inst_idx + l_ - k) : l_end);    // (Inclusive) left-end of the k-mer sliding window.
                const std::size_t window_end = (min_inst_idx + k <= r_end ? min_inst_idx + k : r_end);  // (Exclusive) right-end of the k-mer sliding window.
                const std::size_t window_len = window_end - window_start;

                Minimizer_Iterator min_iter(paths.cbegin() + window_start, window_len, k, l_);  // To iterate over each minimizer in the k-mer sliding window.
                cuttlefish::minimizer_t min;    // The current minimizer.
                std::size_t min_idx;    // Index of the current minimizer in the window.
                std::size_t kmer_idx = 0;   // Index of the current k-mer in the window.
                do
                {
                    min_iter.value_at(min, min_idx);
                    if(window_start + min_idx == min_inst_idx)  // The current k-mer has this minimizer instance as its minimizer.
                    {
                        // Retrieve the k-mer and this minimizer-instance's offset into the instance block.
                        get_kmer(window_start + kmer_idx, kmer);
                        kmers.emplace_back(kmer), inst_idx.emplace_back(i - idx_begin);

                        if(kmers.size() == buf_elem_th)
                            locked_dump();

                        num_kmer++;
                    }

                    kmer_idx++;
                }
                while(++min_iter);
            }

            num_min++;
        }
    }

    if(!kmers.empty())
        locked_dump();
}


template <uint16_t k>
void Kmer_Index<k>::construct_overflow_kmer_mphf()
{
    const boomphf::file_binary<Kmer<k>> kmer_file(overflow_kmers_path().c_str());
    const auto data_iterator = boomphf::range(kmer_file.begin(), kmer_file.end());

    kmer_mphf.emplace(overflow_kmer_count_, data_iterator, working_dir, worker_count, gamma);
    kmer_mphf->save(serialize_stream);

    std::cout <<    "Constructed an MPHF over the overflown k-mers.\n"
                    "\tBits / k-mer: " << (static_cast<double>(kmer_mphf->totalBitSize()) / overflow_kmer_count_) << ".\n";
}


template <uint16_t k>
void Kmer_Index<k>::map_overflown_kmers()
{
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(max_inst_count_)));   // The minimizer-instance blocks are 0-indexed–no need for +1.
    assert(bits_per_entry > 0);
    overflow_kmer_map.emplace(bits_per_entry, overflow_kmer_count_);

    std::FILE* const kmer_file = std::fopen(overflow_kmers_path().c_str(), "rb");
    std::FILE* const inst_id_file = std::fopen(overflow_min_insts_path().c_str(), "rb");
    constexpr std::size_t buf_elem_count = buf_sz_th / (sizeof(Kmer<k>) + sizeof(std::size_t)); // Maximum number of (k-mer, inst-ID) pair in memory per worker.

    std::vector<std::thread> worker;
    worker.reserve(worker_count);

    for(uint16_t t_id = 0; t_id < worker_count; ++t_id)
        worker.emplace_back(
            [this, kmer_file, inst_id_file]()
            {
                auto& kmer_map = *overflow_kmer_map;
                Kmer<k>* const kmer_buf = allocate<Kmer<k>>(buf_elem_count);
                std::size_t* const inst_id_buf = allocate<std::size_t>(buf_elem_count);

                while(true)
                {
                    lock.lock();

                    const std::size_t kmers_read = std::fread(static_cast<void*>(kmer_buf), sizeof(Kmer<k>), buf_elem_count, kmer_file);
                    const std::size_t inst_ids_read = std::fread(static_cast<void*>(inst_id_buf), sizeof(std::size_t), buf_elem_count, inst_id_file);
                    if(kmers_read != inst_ids_read)
                    {
                        std::cerr << "Error reading from the overflow data files. Aborting.\n";
                        std::exit(EXIT_FAILURE);
                    }

                    lock.unlock();

                    if(!kmers_read)
                        break;

                    for(std::size_t i = 0; i < kmers_read; ++i)
                        kmer_map[kmer_mphf->lookup(kmer_buf[i])] = inst_id_buf[i];
                }

                deallocate(kmer_buf);
                deallocate(inst_id_buf);
            }
        );

    for(auto& w : worker)
    {
        if(!w.joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        w.join();
    }

    std::fclose(kmer_file);
    std::fclose(inst_id_file);

    overflow_kmer_map->serialize(serialize_stream);

    std::cout << "Mapped the overflown k-mers to their associated minimizer-instances.\n";
}


template <uint16_t k>
void Kmer_Index<k>::close_output() const
{
    if(!remove_file(overflow_kmers_path()) || !remove_file(overflow_min_insts_path()))
    {
        std::cerr << "Error removing temporary index files. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    serialize_stream.close();
    if(!serialize_stream)
    {
        std::cerr << "Error writing to the k-mer index file at " << index_file_path(output_pref) << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
