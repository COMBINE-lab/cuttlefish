/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.1
  Date   : 2022-01-04
*/

#ifndef _SUPER_KMERS_SERIALIZATION_HELPER_H
#define _SUPER_KMERS_SERIALIZATION_HELPER_H

#include "defs.h"
#include "kmer.h"
#include "libs/refresh/bitmemory/bitmemory.h"
#include "libs/refresh/memory_chunk/memory_chunk.h"

template<unsigned SIZE>
class PackedSuperKmer {
	CKmer<2 * SIZE> data; //super-k-mer must fill in 2 * SIZE
public:
	//idx: 0 -> most significant
	void set_8bytes(uint64_t idx, uint64_t v) {
		data[2 * SIZE - 1 - idx] = v;
	}
	void set_tail(uint64_t idx, uint64_t v, uint64_t n_bits)
	{
		data[2 * SIZE - 1 - idx] = v << (64 - n_bits);
	}
	void clear()
	{
		data.clear();
	}

	//for super-k-mer dump
	const unsigned long long* raw_data() const {
		return data.data;
	}
};

struct SuperKmerSerializationHelper {
public:

	template<unsigned SIZE>
	static void load_super_kmer(refresh::bit_memory_reader<refresh::memory_chunk>& in,
		uint32_t n,
		uint32_t mmer_pos,
		uint32_t signature_len,
		uint32_t mmer,
		PackedSuperKmer<SIZE>& super_kmer) {
		//load part before m-mer
		auto tot_full_uint64 = mmer_pos / 32;
		size_t i = 0;
		for (; i < tot_full_uint64; ++i)
			super_kmer.set_8bytes(i, in.read_8bytes());

		auto tail_symbols = mmer_pos - 32 * tot_full_uint64;
		uint64_t x = in.read_bits(2 * tail_symbols);
		x <<= 64 - 2 * tail_symbols;
		//load m-mer
		uint64_t free_symbols_after_mmer_added;
		if (signature_len > 32 - tail_symbols) {
			auto mmer_prefix_len = 32 - tail_symbols;
			auto mmer_suffix_len = signature_len - mmer_prefix_len;
			uint64_t mmer_prefix = (uint64_t)mmer >> (2 * mmer_suffix_len);
			x += mmer_prefix;
			super_kmer.set_8bytes(i, x);
			++i;
			x = (uint64_t)mmer << (64 - 2 * mmer_suffix_len);
			free_symbols_after_mmer_added = 32 - mmer_suffix_len;
		}
		else {
			free_symbols_after_mmer_added = 32 - tail_symbols - signature_len;
			x += (uint64_t)mmer << (2 * free_symbols_after_mmer_added);
		}

		n -= mmer_pos + signature_len;
		tot_full_uint64 = n / 32;
		auto head_symbols = n - 32 * tot_full_uint64;

		if (head_symbols == free_symbols_after_mmer_added) {
			//head_symbols may be 0 and it should be OK
			x += in.read_bits(2 * head_symbols);
			super_kmer.set_8bytes(i++, x);

			for (size_t j = 0; j < tot_full_uint64; ++j)
				super_kmer.set_8bytes(i++, in.read_8bytes());
		}
		else if (head_symbols < free_symbols_after_mmer_added) {
			auto part_size_bits = 2 * (free_symbols_after_mmer_added - head_symbols);
			auto tmp = in.read_bits(2 * head_symbols);
			x += tmp << part_size_bits;
			for (size_t j = 0; j < tot_full_uint64; ++j) {
				tmp = in.read_8bytes();
				x += tmp >> (64 - part_size_bits);
				super_kmer.set_8bytes(i++, x);
				x = tmp << part_size_bits;
			}
			super_kmer.set_8bytes(i++, x);
		}
		else {
			auto tmp = in.read_bits(2 * head_symbols);
			auto part_size_bits = 2 * (head_symbols - free_symbols_after_mmer_added);
			x += tmp >> part_size_bits;
			super_kmer.set_8bytes(i++, x);

			x = tmp << (64 - part_size_bits);
			for (size_t j = 0; j < tot_full_uint64; ++j) {
				tmp = in.read_8bytes();
				x += tmp >> part_size_bits;
				super_kmer.set_8bytes(i++, x);
				x = tmp << (64 - part_size_bits);
			}
			super_kmer.set_8bytes(i++, x);
		}
	}
};

#endif
