
#include "Graph_Partitioner.hpp"
#include "Subgraphs_Manager.hpp"
#include "Minimizer_Iterator.hpp"
#include "DNA_Utility.hpp"
#include "Spin_Lock.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "RabbitFX/io/Reference.h"
#include "RabbitFX/io/Formater.h"
#include "RabbitFX/io/Globals.h"
#include "parlay/parallel.h"

#include <functional>
#include <memory>
#include <utility>
#include <thread>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k, bool Is_FASTQ_, bool Colored_>
Graph_Partitioner<k, Is_FASTQ_, Colored_>::Graph_Partitioner(Subgraphs_Manager<k, Colored_>& subgraphs, const Data_Logistics& logistics, const uint16_t l):
      subgraphs(subgraphs)
    //, seqs(logistics.input_paths_collection())
    , l_(l)
    , sup_km1_mer_len_th(2 * (k - 1) - l_)
    , chunk_pool_sz(parlay::num_workers() * (Is_FASTQ_ ? 2 : 4))    // TODO: make a more informed choice.
    , chunk_pool(chunk_pool_sz)
    , chunk_q(chunk_pool_sz)
    , parsed_chunk_w(parlay::num_workers())
    , reader_c(parlay::num_workers() < 32 ? 2 : 4)
    , stat_w(parlay::num_workers())
{
    auto& input_paths = logistics.input_paths_collection();
    //std::cerr << "num input paths = " << std::distance(input_paths.begin(), input_paths.end()) << "\n";
    seqs = std::deque<std::string>(input_paths.begin(), input_paths.end());
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::partition()
{
    // const bool large_src = true;    // Whether dealing with large individual source-files or not.   TODO: fix.

    double t_part = 0;
    double t_collate = 0;

    // Number of consumers when partitioning is done in the producer-consumer model.
    // const auto consumer_c = (parlay::num_workers() > (reader_c - 1) ? parlay::num_workers() - (reader_c - 1) : 1);

/*
    if constexpr(!Colored_)
    {
        std::thread reader([&](){ read_chunks(); });

        parlay::parallel_for(0, consumer_c,
            [&](auto)
            {
                process_uncolored_chunks();
            }, 1);

        reader.join();
    }
    else if(large_src)
    {
        std::thread reader([&](){ read_chunks(); });

        std::cerr << "\n";

        uint64_t bytes_processed = 0;
        while(!m_pushed_all_data)
        {
            const auto t_0 = timer::now();
            bytes_consumed = 0;

            source_id_t min_source = std::numeric_limits<source_id_t>::max();
            source_id_t max_source = 0;
            Spin_Lock source_lock;
            auto last_cp_count = last_checkpoint.load();
            if (bytes_processed > last_checkpoint.load()) {
                last_checkpoint.compare_exchange_strong(last_cp_count, bytes_processed);
            }

            m_do_reading = true;
            parlay::parallel_for(0, consumer_c,
                [&](auto)
                {
                    source_id_t smallest_src = std::numeric_limits<source_id_t>::max();
                    source_id_t largest_src = 0;
                    process_colored_chunks(smallest_src, largest_src);
                    source_lock.lock();
                    min_source = std::min(smallest_src, min_source);
                    max_source = std::max(largest_src, max_source);
                    source_lock.unlock();

                }, 1);

            const auto t_1 = timer::now();
            t_part += timer::duration(t_1 - t_0);

            // Collate and flush all buckets.
            if(max_source > 0)
                subgraphs.collate_super_kmer_buffers(min_source, max_source);
            bytes_processed += bytes_consumed;
            const auto t_2 = timer::now();
            t_collate += timer::duration(t_2 - t_1);

            std::cerr << "\rProcessed " << (bytes_processed / (1024 * 1024)) << "MB of uncompressed input data. Source range [" << min_source << ", " << max_source << "]";
        }

        std::cerr << "\n";
        reader.join();
    }
    else
*/
    {
        std::vector<std::pair<std::size_t, std::size_t>> sz_src;    // Size of the sources and their IDs.
        sz_src.reserve(seqs.size());
        for(std::size_t s = 0; s < seqs.size(); ++s)
            sz_src.emplace_back(file_size(seqs[s]), s);

        // Process sources in decreasing order of size.
        std::sort(sz_src.begin(), sz_src.end(), std::greater<>());


        std::size_t src_idx = 0;
        Spin_Lock src_lock;
        bytes_consumed = 0;
        constexpr std::size_t log_step = 512 * 1024 * 1024; // 512 MB.
        std::size_t log_checkp = log_step;
        parlay::parallel_for(0, parlay::num_workers(),
        [&](auto)
        {
            const auto w = parlay::worker_id();
            typedef typename RabbitFX_DS_type<Is_FASTQ_>::reader_t reader_t;
            std::unique_ptr<reader_t> reader;

            while(true)
            {
                src_lock.lock();

                std::size_t src;    // Next source to partition.
                const bool done = (src_idx == sz_src.size());
                if(!done)
                    src = sz_src[src_idx++].second;

                if(bytes_consumed >= log_checkp)
                {
                    std::cerr << "\rPartitioned " << (bytes_consumed / (1024 * 1024)) << " MB of uncompressed data.";
                    log_checkp = bytes_consumed + log_step;
                }

                src_lock.unlock();

                if(done)
                    break;

                if(!reader)
                    reader = std::make_unique<reader_t>(seqs[src], chunk_pool);
                else
                    reader->set_new_file(seqs[src]);

                chunk_t* chunk;
                while(true)
                {
                    if constexpr(Is_FASTQ_)
                        chunk = reader->readNextChunk();
                    else
                        chunk = reader->readNextChunkList();

                    if(chunk == NULL)
                        break;

                    bytes_consumed += process_chunk(chunk, src + 1);
                }

                if constexpr(Colored_)
                    subgraphs.flush_worker_if_req(w);
            }
        }, 1);
    }

    std::cerr << "\rPartitioned " << (bytes_consumed / (1024 * 1024)) << " MB of uncompressed data.\n";


    Worker_Stats stat;
    std::for_each(stat_w.cbegin(), stat_w.cend(), [&](const auto& s){ stat += s.unwrap(); });
    std::cerr << "Number of processed chunks: " << stat.chunk_count << ".\n";
    std::cerr << "Total size of chunks: " << stat.chunk_bytes << ".\n";
    std::cerr << "Number of records: " << stat.record_count << ".\n";
    std::cerr << "Number of super (k - 1)-mers: " << stat.weak_super_kmer_count << ".\n";
    std::cerr << "Total length of the weak super k-mers:  " << stat.weak_super_kmers_len << ".\n";
    std::cerr << "Total length of the super (k - 1)-mers: " << stat.super_km1_mers_len << ".\n";
    std::cerr << "Total work in parse: " << stat.parse_time << "s.\n";
    std::cerr << "Total work in processing records: " << stat.process_time << "s.\n";
    std::cerr << "Max work in processing records: " <<
        [&](){ double t = 0; std::for_each(stat_w.cbegin(), stat_w.cend(), [&](const auto& s){ t = std::max(t, s.unwrap().process_time); }); return t; }() << "s.\n";
    if constexpr(Colored_)
    {
        std::cerr << "Time taken in processing colored chunks: " << t_part << "s.\n";
        std::cerr << "Time taken in collating colored chunks:  " << t_collate << "s.\n";
    }
}

template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::read_chunks()
{
    const auto t_s = timer::now();
    std::atomic<uint64_t> chunk_count{0};   // Number of read chunks.
    std::atomic<uint64_t> bytes_pushed{0};
    std::atomic<uint64_t> last_update{0};
    //std::atomic<uint64_t> last_checkpoint{0};
    rabbit::int64 source_id = 1;    // Source (i.e. file) ID of a chunk.
    size_t n_files = seqs.size();

    // mutex to control obtaining the next file to read from the queue
    Spin_Lock fp_mutex;
    std::vector<std::thread> reader_threads;

    // TODO: make the number of readers dynamic and based on the workload between 
    // the reading and parsing / super-kmerization
    for (size_t t = 0; t < reader_c; ++t) {
        reader_threads.push_back(std::thread([&]{
            using reader_t = typename RabbitFX_DS_type<Is_FASTQ_>::reader_t;
            std::unique_ptr<reader_t> reader;
            // keep reading (we will break when we have read all of the 
            // data from all of the files).
            while(true) {
                // the `m_do_reading` variable controls if we want to pause 
                // reading of new files (to give the super-kmerization threads
                // a chance to process data that has already been read and avoid
                // the read queue growing too large).
                if (!Colored_ || m_do_reading) {
                    // get the next input file's name
                    std::string next_file;
                    fp_mutex.lock();
                    bool done = seqs.empty();
                    size_t current_source_id = source_id;
                    if (!done) {
                        source_id++;
                        max_read_source_id++;
                        next_file = seqs.front();
                        seqs.pop_front();
                    }
                    fp_mutex.unlock();
                    // if there are no more files to read, then this
                    // thread is done
                    if (done) { m_do_reading = false; break; }

                    // if we don't have a reader for this thread yet, then 
                    // make one, otherwise assign the existing reader the 
                    // next file.
                    if (!reader) { 
                        reader = std::make_unique<reader_t>(next_file, chunk_pool);
                    } else {
                        reader->set_new_file(next_file);
                    }

                    // lambda to handle pushing a chunk (either FASTQ or FASTA)
                    // to the data chunk queue.
                    const auto push_chunk = [&](auto const chunk)
                    {
                        if(chunk == NULL)
                            return false;

                        // for FASTA chunks
                        if constexpr(!Is_FASTQ_) {
                            // NOTE: The metadata in the RabbitFX::io::FastaChunk
                            // is not filled in when it is read (i.e. start, end, nseqs).
                            // So, to get the total number of bytes here, we actually have to 
                            // walk the list of data chunks produced.

                            if (chunk->chunk == nullptr) {
                                return false;
                            }
                            // get the first RabitFX::io::core::DataChunk
                            auto* dchunk = chunk->chunk;
                            uint64_t nb = 0;

                            // while we've not reached the end of the list, accumulate
                            // the total number of bytes.
                            while (dchunk != nullptr) {
                                nb += dchunk->size;
                                dchunk = dchunk->next;
                            }

                            // keeping track of the number of bytes pushed into
                            // the queue.
                            bytes_pushed += nb;

                            // if we've processed > bytes_per_batch since the last check, 
                            // then print it out and pause reading of further files until
                            // these files are done.
                            // NOTE: multiple threads could enter this condition before we set
                            // m_do_reading = false; may not be a problem, but think about this.
                            if (bytes_pushed > (last_checkpoint + bytes_per_batch)) {
                                last_checkpoint.store(bytes_pushed.load());
                                //std::cerr << "::Pushed " << (bytes_pushed/ (1024 * 1024)) << "MB of uncompressed input data.\n";
                                m_do_reading = false;
                            }
                        } else {
                            bytes_pushed += chunk->size;
                        }
                        // regardless of the chunk type, push it.
                        chunk_q.Push(current_source_id, chunk);
                        chunk_count++;

                        if constexpr(!Colored_) {
                            if (bytes_pushed >= (last_update + 512 * 1024 * 1024)) {
                                const auto b = bytes_pushed.load();
                                last_update.store(b);
                                std::cerr << "\rPushed " << b / (1024 * 1024) << "MB of parsed data onto chunk queue.";
                            }
                        }
                        return true;
                    };

                    // for the current file, read all of the chunks and 
                    // push them onto the queue.
                    bool chunks_remain = true;
                    while(chunks_remain) {
                        if constexpr(Is_FASTQ_) {
                            chunks_remain = push_chunk(reader->readNextChunk());
                        } else {
                            chunks_remain = push_chunk(reader->readNextChunkList());
                        }
                    }
                    //m_do_reading = true;

                    // increment the source id
                    //source_id++;
                    // keep track of the maximum source
                    // id we've read so far.
                    //max_read_source_id++;
                }
            }
        }));
    }

    // wait for all of our reader threads to 
    // finish.
    for (auto& r: reader_threads) {
        r.join();
    }

    // we're done with the queue, so signal that nothing further 
    // will be pushed.
    m_do_reading = false;
    chunk_q.SetCompleted();
    m_pushed_all_data = true;
    const auto t_e = timer::now();
    std::cerr << "\n\nRead " << chunk_count << " chunks in total from " << n_files << " files. Time elapsed: " << timer::duration(t_e - t_s) << " s.\n";
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::process_uncolored_chunks()
{
    assert(!Colored_);

    rabbit::int64 source_id = 0;    // Source (i.e. file) ID of a chunk.
    chunk_t* chunk = nullptr; // Current chunk.
    rabbit::int64 last_source = 0;  // Source ID of the last chunk processed.
    (void)last_source;

    // while things can still be put into the queue
    while(!chunk_q.IsCompleted()) {
        // try and get a chunk out
        // bool saw_max_id = false;
        while(chunk_q.Pop(source_id, chunk)) {
            // assert(source_id >= last_source);

            bytes_consumed += process_chunk(chunk, source_id);

            last_source = source_id;
            if (last_source >= max_read_source_id - 1) { 
                m_do_reading = true; 
            } else {
                //std::cerr << "last_source = " << last_source << ", max_read_source_id = " << max_read_source_id << "\n";
            }
        }
    }
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
bool Graph_Partitioner<k, Is_FASTQ_, Colored_>::process_colored_chunks(source_id_t& min_source_id, source_id_t& max_source_id)
{
    assert(Colored_);

    rabbit::int64 source_id;    // Source (i.e. file) ID of a chunk.
    chunk_t* chunk; // Current chunk.
    rabbit::int64 last_source = 0;  // Source ID of the last chunk processed.
    (void)last_source;

    // while there may still be things pushing into the queue 
    // or the queue still has things in it
    while(m_do_reading || !chunk_q.IsEmpty())
    {
        // if we can get something out
        if(chunk_q.TryPop(source_id, chunk)) {
            //assert(source_id >= last_source);

            bytes_consumed += process_chunk(chunk, source_id);

            last_source = source_id;
            if (last_source < min_source_id) { min_source_id = last_source; }
            if (last_source > max_source_id) { max_source_id = last_source; }
        } else if (m_pushed_all_data) {
            return false;
        }
    }

    //if (saw_last_source) { m_do_reading = true; }
    return !chunk_q.IsCompleted();
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
uint64_t Graph_Partitioner<k, Is_FASTQ_, Colored_>::process_chunk(chunk_t* chunk, const source_id_t source_id)
{
    auto& w_stat = stat_w[parlay::worker_id()].unwrap();

    uint64_t sup_km1_mer_count = 0; // Number of parsed super (k - 1)-mers.
    std::size_t sup_km1_mers_len = 0;   // Total length of the super (k - 1)-mers, in bases.
    std::size_t weak_sup_kmers_len = 0; // Total length of the (weak) super k-mers, in bases.
    uint64_t chunk_bytes = 0;   // Count of bytes in the chunk.

    // Minimizer_Iterator<const char*, k - 1, true> min_it(l_, min_seed);   // `l`-minimizer iterator for `(k - 1)`-mers.
    Min_Iterator<k - 1> min_it(l_); // `l`-minimizer iterator for `(k - 1)`-mers.

    const auto t_0 = timer::now();
    auto& parsed_chunk = parsed_chunk_w[parlay::worker_id()].unwrap();
    parsed_chunk.clear();
    if constexpr(Is_FASTQ_)
        chunk_bytes = chunk->size,
        w_stat.record_count += rabbit::fq::chunkFormat(chunk, parsed_chunk);
    else
    {
        w_stat.record_count += rabbit::fa::chunkListFormat(*chunk, parsed_chunk);

        auto ptr = chunk->chunk;
        do
        {
            chunk_bytes += ptr->size;
            chunk_pool.Release(ptr);
            ptr = ptr->next;
        }
        while(ptr != NULL);

        rabbit::fa::FastaFileReader::release_chunk_list(chunk);
    }

    const auto t_1 = timer::now();
    w_stat.parse_time += timer::duration(t_1 - t_0);

    for(auto& record : parsed_chunk)
    {
        const char* seq;
        std::size_t seq_len;

        if constexpr(Is_FASTQ_)
            seq = reinterpret_cast<const char*>(record.base + record.pseq),
            seq_len = record.lseq;
        else
        {
            record.seq.resize(record.seq.length() + 32);    // Deal with UB in `Super_Kmer_Chunk::add_encoded_label`
            seq = record.seq.c_str(),
            seq_len = record.length;
        }

        std::size_t last_frag_end = 0;  // Ending index (exclusive) of the last sequence fragment.
        while(true)
        {
            std::size_t frag_beg;   // Offset of the next fragment in the sequence.
            std::size_t frag_len;   // Length of the fragment.

            // Skip placeholder bases.
            for(frag_beg = last_frag_end; frag_beg + (k + 1) <= seq_len; ++frag_beg)
                if(DNA_Utility::is_DNA_base(seq[frag_beg]))
                    break;

            if(frag_beg + (k + 1) > seq_len)    // No more fragment remains with complete (k + 1)-mers, i.e. edges.
                break;

            const auto frag = seq + frag_beg;   // Current fragment.

            // Check whether the first (k + 1)-mer of the fragment has any placeholder bases.
            for(frag_len = 1; frag_len < (k + 1); ++frag_len)
                if(!DNA_Utility::is_DNA_base(frag[frag_len]))
                    break;

            if(frag_len < (k + 1))
            {
                last_frag_end = frag_beg + frag_len;
                continue;
            }


            // minimizer_t cur_min;    // Minimizer of the current super (k - 1)-mer in the iteration.
            // minimizer_t next_min;   // Minimizer of the next super (k - 1)-mer in the iteration.
            uint64_t cur_h; // 64-bit hash of the current super (k - 1)-mer's minimizer.
            uint64_t next_h;    // 64-bit hash of the next super (k - 1)-mer's minimizer.
            std::size_t prev_g; // Subgraph ID of the previous super (k - 1)-mer's minimizer.
            std::size_t cur_g;  // Subgraph ID of the current super (k - 1)-mer's minimizer.
            std::size_t next_g; // Subgraph ID of the next super (k - 1)-mer's minimizer.
            // std::size_t cur_min_off, next_min_off;  // Relative offsets of the minimizers in the fragment. Used for assertion checks.
            std::size_t cur_sup_km1_mer_off = 0;    // Relative offset of the current super (k - 1)-mer in the fragment.
            std::size_t km1_mer_idx = 0;    // Index of the current (k - 1)-mer in the current super (k - 1)-mer.
            frag_len = k - 1;

            // min_it.reset(frag, seq_len - frag_beg); // The fragment length is an estimate; upper-bound to be exact.
            // min_it.value_at(cur_min, cur_min_off, cur_h);
            min_it.reset(frag);
            cur_h = min_it.hash();
            cur_g = subgraphs.graph_ID(cur_h);
            prev_g = subgraphs.graph_count();   // To deal with false-positive `-Wmaybe-uninitialized` later on.

            while(DNA_Utility::is_DNA_base(frag[frag_len]))
            {
                const auto len = km1_mer_idx + (k - 1); // Length of the current super (k - 1)-mer.

                min_it.advance(frag[frag_len]);
                km1_mer_idx++, frag_len++;

                // min_it.value_at(next_min, next_min_off, next_h);
                next_h = min_it.hash();
                next_g = subgraphs.graph_ID(next_h);
/*                  assert(next_min_off >= cur_sup_km1_mer_off + km1_mer_idx);

                if(next_min_off != cur_min_off)
                    // Either the last minimizer just fell out of the (k - 1)-mer window or the new minimizer sits at the last l-mer.
                    assert( cur_min_off == cur_sup_km1_mer_off + km1_mer_idx - 1 ||
                            next_min_off == cur_sup_km1_mer_off + km1_mer_idx + (k - 1) - l_);
*/
                // Either encountered a discontinuity vertex—the k-mer whose suffix is the current (k - 1)-mer, or the super (k - 1)-mer extends too long.
                if(next_g != cur_g || len == sup_km1_mer_len_th)
                {
                    if(next_g != cur_g)
                        // The `(k - 1)`-mers of this discontinuity k-mer have minimizers mapping to different subgraphs.
                        assert(subgraphs.G().is_discontinuity(frag + cur_sup_km1_mer_off + km1_mer_idx - 1));

                    const auto next_sup_km1_mer_off = cur_sup_km1_mer_off + km1_mer_idx;
                    assert(next_sup_km1_mer_off == frag_len - (k - 1));
                    sup_km1_mers_len += len;
                    sup_km1_mer_count++;

                    const bool l_joined = (cur_sup_km1_mer_off > 0);    // Whether this weak super k-mer is joined to the one to its left.
                    const bool r_joined = true; // Whether this weak super k-mer is joined to the one to its right.
                    const bool l_disc = (l_joined && prev_g != cur_g);  // Whether it's left-discontinuous.
                    const bool r_disc = (next_g != cur_g);  // Whether it's right-discontinuous.
                    // const bool l_cont = (l_joined && prev_g == cur_g);  // Whether it's left-continuous.
                    // const bool r_cont = (next_g == cur_g);  // Whether it's right-continuous.
                    const auto len_weak = l_joined + len + r_joined;    // Length of the weak super k-mer.
                    assert(len_weak >= k);
                    // TODO: the following add, being to different subgraphs' different worker-buffers, causes lots of cache misses.
                    if constexpr(!Colored_)
                        subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, l_disc, r_disc);
                    else
                        subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, source_id, l_disc, r_disc);
                    weak_sup_kmers_len += len_weak;

                    cur_sup_km1_mer_off = next_sup_km1_mer_off;
                    prev_g = cur_g;
                    cur_g = next_g;
                    km1_mer_idx = 0;
                }

                // cur_min = next_min;
                // cur_min_off = next_min_off; // The minimizer-instance offsets are tracked only for assertion checks.
                cur_h = next_h;
            }

            const auto len = frag_len - cur_sup_km1_mer_off;
            sup_km1_mers_len += len;
            sup_km1_mer_count++;

            const bool l_joined = (cur_sup_km1_mer_off > 0);
            const bool r_joined = false;
            const bool l_disc = (l_joined && prev_g != cur_g);
            const bool r_disc = false;
            // const bool l_cont = (l_joined && prev_g == cur_g);
            // const bool r_cont = false;
            const auto len_weak = l_joined + len + r_joined;
            assert(len_weak >= k);
            // TODO: the following add, being to different subgraphs' different worker-buffers, causes lots of cache misses.
            if constexpr(!Colored_)
                subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, l_disc, r_disc);
            else
                subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, source_id, l_disc, r_disc);
            weak_sup_kmers_len += len_weak;

            last_frag_end = frag_beg + frag_len;
        }
    }

    const auto t_2 = timer::now();
    w_stat.process_time += timer::duration(t_2 - t_1);

    if constexpr(Is_FASTQ_) // In case of FASTA, chunks have been released earlier during bytes-counting.
        chunk_pool.Release(chunk);

    w_stat.chunk_count++;
    w_stat.chunk_bytes += chunk_bytes;
    w_stat.weak_super_kmer_count += sup_km1_mer_count;
    w_stat.weak_super_kmers_len += weak_sup_kmers_len;
    w_stat.super_km1_mers_len += sup_km1_mers_len;

    return chunk_bytes;
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::Worker_Stats::operator+=(const Worker_Stats& rhs)
{
    chunk_count += rhs.chunk_count;
    chunk_bytes += rhs.chunk_bytes;
    record_count += rhs.record_count;
    weak_super_kmer_count += rhs.weak_super_kmer_count;
    weak_super_kmers_len += rhs.weak_super_kmers_len;
    super_km1_mers_len += rhs.super_km1_mers_len;
    parse_time += rhs.parse_time;
    process_time += rhs.process_time;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL_L2, cuttlefish::Graph_Partitioner)
