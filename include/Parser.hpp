
#ifndef PARSER_HPP
#define PARSER_HPP



#include "utility.hpp"
#include "RabbitFX/io/FastxChunk.h"

#include <cstddef>
#include <string>
#include <iostream>
#include <atomic>


namespace cuttlefish
{

// A class to parse FASTX files.
class Parser
{

private:

    typedef rabbit::fq::FastqDataChunk chunk_t; // Type of chunks containing sequences.
    typedef rabbit::fq::FastqDataPool fq_chunk_pool_t;  // Type of memory pool for chunks.
    typedef rabbit::core::TDataQueue<chunk_t> fq_chunk_queue_t; // Type of queue of chunks ready for consumption.

    std::string file_path;  // File to parse.

    std::size_t consumer_count; // Number of concurrent consumers of the read sequence data.

    std::atomic_uint64_t record_count;  // Number of FASTX records in the input.

    // TODO: document.
    struct timing_info
    {
        double q_wait_time = 0;
        double chunk_format_time = 0;
        double min_it_init_time = 0;
        double min_it_iter_time = 0;


        friend std::ostream& operator<<(std::ostream& os, const timing_info& t)
        {
            os << "Queue-wait time:    " << t.q_wait_time << "s.\n";
            os << "Chunk format time:  " << t.chunk_format_time << "s.\n";
            os << "Iterator init time: " << t.min_it_init_time << "s.\n";
            os << "Iterator iter time: " << t.min_it_iter_time << "s.\n";

            return os;
        }

        void operator+=(const timing_info& rhs)
        {
            q_wait_time += rhs.q_wait_time;
            chunk_format_time += rhs.chunk_format_time;
            min_it_init_time += rhs.min_it_init_time;
            min_it_iter_time += rhs.min_it_iter_time;
        }
    };

    std::vector<Padded<timing_info>> T;


    // Proof-of-concept production method of parsed sequences.
    void produce(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q);

    // Proof-of-concept consumption method for parsed sequences.
    void consume_count_bases(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q, std::vector<std::atomic_uint64_t>& count);

    // Proof-of-concept consumption method for parsed sequences.
    void consume_split_super_kmers(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q, std::atomic_uint64_t& count, timing_info& t);

public:

    // Constructs a parser for the file at path `file_path` for `parser_count`
    // parsers.
    Parser(const std::string& file_path, std::size_t parser_count = 1);

    // Proof-of-concept parse method.
    void parse();
};

}



#endif
