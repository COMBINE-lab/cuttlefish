
#ifndef GRAPH_PARTITIONER_HPP
#define GRAPH_PARTITIONER_HPP



#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "RabbitFX/io/FastxChunk.h"
#include "RabbitFX/io/DataQueue.h"
#include "RabbitFX/io/FastxStream.h"
#include "RabbitFX/io/Reference.h"

#include <cstdint>
#include <cstddef>
#include <atomic>
#include <string>
#include <vector>
#include <deque>

namespace cuttlefish
{

template <uint16_t k, bool Colored_> class Subgraphs_Manager;


// Data types for the `rabbitfx` parser. `Is_FASTQ` denotes whether the parsing
// is over FASTQ data or not (i.e. FASTA).
template <bool Is_FASTQ_>
struct RabbitFX_DS_type
{};


template <>
struct RabbitFX_DS_type<true>
{
    typedef rabbit::fq::FastqDataChunk chunk_t; // Type of chunks containing read sequences.
    typedef rabbit::fq::FastqDataPool chunk_pool_t; // Type of memory pools for chunks.
    typedef rabbit::core::TDataQueue<chunk_t> chunk_q_t;    // Type of queue of read chunks.
    typedef rabbit::fq::FastqFileReader reader_t;   // Type of file-reader.
    typedef neoReference ref_t; // Type of parsed data.
};


template <>
struct RabbitFX_DS_type<false>
{
    typedef rabbit::fa::FastaChunk chunk_t; // Type of chunks containing read sequences.
    typedef rabbit::fa::FastaDataPool chunk_pool_t; // Type of memory pools for chunks.
    typedef rabbit::core::TDataQueue<chunk_t> chunk_q_t;    // Type of queue of read chunks.
    typedef rabbit::fa::FastaFileReader reader_t;   // Type of file-reader.
    typedef Reference ref_t;    // Type of parsed data.
};


// =============================================================================
// Class to partition de Bruijn graphs into subgraphs based on minimizers of the
// `k`-mers. Effectively, splits input sequences to maximal weak super k-mers
// and distributes those to appropriate subgraphs. `Is_FASTQ_` denotes whether
// the input is FASTQ or not (i.e. FASTA). `Colored_` denotes whether the
// vertices in the graph have associated colors.
// TODO: only supports FASTQ for now; extend to FASTA also.
template <uint16_t k, bool Is_FASTQ_, bool Colored_>
class Graph_Partitioner
{

private:

    Subgraphs_Manager<k, Colored_>& subgraphs;  // Subgraphs of the de Bruijn graph.

    typedef typename RabbitFX_DS_type<Is_FASTQ_>::chunk_t chunk_t;  // Type of chunks containing read sequences.
    typedef typename RabbitFX_DS_type<Is_FASTQ_>::chunk_pool_t chunk_pool_t;    // Type of memory pools for chunks.
    typedef typename RabbitFX_DS_type<Is_FASTQ_>::chunk_q_t chunk_q_t;  // Type of queue of read chunks.
    typedef typename RabbitFX_DS_type<Is_FASTQ_>::ref_t parsed_rec_t;   // Type of parsed records.

    std::deque<std::string> seqs;    // Input sequence collection.
    std::atomic_bool m_do_reading{true}; // Signal if it's ok to continue reading input, or if we should wait
    std::atomic_bool m_pushed_all_data{false}; // Signal if it's ok to continue reading input, or if we should wait
    std::atomic<uint64_t> last_checkpoint{0};
    std::atomic<int64_t> max_read_source_id{1};

    const uint16_t l_;  // Size of minimizers for the super k-mers.
    const std::size_t sup_km1_mer_len_th;   // Length threshold of super (k - 1)-mers.

    const std::size_t chunk_pool_sz;    // Maximum number of chunks in the chunk memory pool.
    chunk_pool_t chunk_pool;    // Memory pool for chunks of sequences.
    chunk_q_t chunk_q;  // Read chunks.

    std::vector<Padded<std::vector<parsed_rec_t>>> parsed_chunk_w;  // Parsed record collection per worker.

    std::atomic_uint64_t bytes_consumed;    // Counts of input bytes consumed across all workers in one batch in the colored-case.
    constexpr static uint64_t bytes_per_batch = 1024 * 1024 * 1024lu;   // 1GB per input batch, at least.

    const std::size_t reader_c; // Number of working doing input-reads.

    struct Worker_Stats
    {
        uint64_t chunk_count = 0;   // Number of chunks processed from the input.
        uint64_t chunk_bytes = 0;   // Total size of chunks in bytes.
        uint64_t record_count = 0;  // Number of records in the sequences.

        uint64_t weak_super_kmer_count = 0; // Number of weak super k-mers in the sequences.
        uint64_t weak_super_kmers_len = 0;  // Total length of the weak super k-mers in the sequences.
        uint64_t super_km1_mers_len = 0;    // Total length of the super (k - 1)-mers in the sequences.

        double parse_time = 0;  // Total time taken in parsing read chunks.
        double process_time = 0;    // Total time taken in processing parsed records.

        Worker_Stats()
        {}

        void operator+=(const Worker_Stats& rhs);
    };

    std::vector<Padded<Worker_Stats>> stat_w;   // Sequence-processing statistics per worker.

    // Reads the provided sequences into chunks from the chunk memory pool and
    // puts the read chunks into the read-queue`.
    void read_chunks();

    // Processes the read chunks from the read-queue and returns the processed
    // chunks to the chunk memory pool for uncolored graphs.
    void process_uncolored_chunks();

    // Processes the read chunks from the read-queue and returns the processed
    // chunks to the chunk memory pool for colored graphs until total bytes
    // processed across all workers reaches a certain threshold. Returns
    // `false` if no data remain anymore in the queue after this processing.
    bool process_colored_chunks(source_id_t& min_source, source_id_t& max_source);

    // Processes the chunk `chunk` with source-ID `source_id` and releases it
    // to the chunk-pool. The parsed sequences are stored in `parsed_chunk`.
    // Returns the count of bytes in chunk.
    uint64_t process_chunk(chunk_t* chunk, source_id_t source_id);

public:

    // Constructs a de Bruijn graph partitioner with `l`-minimizers for the
    // sequences from the data logistics manager `logistics`. The graph is
    // partitioned into the subgraph-manager `subgraphs`.
    Graph_Partitioner(Subgraphs_Manager<k, Colored_>& subgraphs, const Data_Logistics& logistics, uint16_t l);

    // Partitions the passed sequences into maximal weak super k-mers and
    // deposits those to corresponding subgraphs.
    void partition();
};

}



#endif
