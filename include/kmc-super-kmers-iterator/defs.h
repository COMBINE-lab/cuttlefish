#pragma once
#include <cinttypes>
#include "libs/refresh/memory_chunk/memory_chunk.h"
#include "libs/refresh/parallel_queues/parallel-queues.h"

constexpr size_t READER_PART_BUFF_SIZE = 1ull << 23;

class reader_data_t {
	std::unique_ptr<uint8_t[]> ptr;
public:
	reader_data_t(reader_data_t&& rhs) = default;
	reader_data_t& operator=(reader_data_t&& rhs) = default;

	reader_data_t(const reader_data_t& rhs) = delete;
	reader_data_t& operator=(const reader_data_t& rhs) = delete;

	refresh::memory_chunk<uint8_t> data;
	reader_data_t(size_t size = READER_PART_BUFF_SIZE) :
		ptr(new uint8_t[size]),
		data(ptr.get(), size) {
	}

};

using bin_reader_super_kmers_packs_maker_queue_t = refresh::parallel_queue<reader_data_t>;

class super_kmers_packer_data_t {
	std::unique_ptr<uint8_t[]> ptr;
public:
	size_t n_super_kmers; 
	using data_t = refresh::memory_chunk<uint8_t>;
	data_t data;

	super_kmers_packer_data_t() = default;

	super_kmers_packer_data_t(size_t n_super_kmers, size_t size) :
		ptr(new uint8_t[size]),
		n_super_kmers(n_super_kmers),
		data(ptr.get(), size) {

	}

	super_kmers_packer_data_t(const super_kmers_packer_data_t& rhs) = delete;
	super_kmers_packer_data_t& operator=(const super_kmers_packer_data_t& rhs) = delete;

	super_kmers_packer_data_t(super_kmers_packer_data_t&& rhs) = default;
	super_kmers_packer_data_t& operator=(super_kmers_packer_data_t&& rhs) = default;

};

using super_kmers_packs_maker_super_kmer_iterator_queue_t = refresh::parallel_queue<super_kmers_packer_data_t>;

//in uint64_t, so 8 means max k = 256 for example
constexpr size_t MAX_KMER_SIZE = 8;