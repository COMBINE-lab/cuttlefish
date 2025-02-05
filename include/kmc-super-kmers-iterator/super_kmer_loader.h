#pragma once
#include <cinttypes>
#include <string>
#include <vector>
#include <iostream>
#include "super_kmers_serialization_helper.h"
#include "libs/refresh/serialization/serialization.h"

class SuperKmerLoader
{
	uint32_t bin_id;
	uint32_t n_bins;
	uint32_t kmer_len;
	uint32_t signature_len;

	std::vector<uint32_t> id_to_mmer;
	uint32_t bits_for_sig_pos;
	uint32_t bits_for_sig_id;
	uint64_t bits_for_sig_pos_mask;
	uint64_t bits_for_sig_id_mask;

	bool encode_plain;

	uint32_t bits_for_sample_id;
	uint32_t bits_for_n_additional_symbols;

	uint32_t bits_required_base;

	uint32_t bits_required_to_represent(uint32_t val)
	{
		uint32_t res{};
		while (val)
		{
			++res;
			val >>= 1;
		}
		return res;
	}

public:

	SuperKmerLoader(std::istream& in) {

		std::string start_token_pat = "SERIALIZER_START";
		std::string end_token_pat = "SERIALIZER_END";

		std::string start_token(start_token_pat.size(), ' ');
		std::string end_token(end_token_pat.size(), ' ');

		in.read((char*)start_token.data(), start_token.size());

		if (start_token != start_token_pat) {
			std::cerr << "Error: wrong serialize start token\n";
			exit(1);
		}

		refresh::load_little_endian(bin_id, in);
		refresh::load_little_endian(n_bins, in);
		refresh::load_little_endian(kmer_len, in);
		refresh::load_little_endian(signature_len, in);



		refresh::load_little_endian(id_to_mmer, in);

		refresh::load_little_endian(bits_for_sig_pos, in);
		refresh::load_little_endian(bits_for_sig_id, in);
		refresh::load_little_endian(bits_for_sig_pos_mask, in);
		refresh::load_little_endian(bits_for_sig_id_mask, in);

		refresh::load_little_endian(encode_plain, in);

		refresh::load_little_endian(bits_for_sample_id, in);
		refresh::load_little_endian(bits_for_n_additional_symbols, in);

		refresh::load_little_endian(bits_required_base, in);

		in.read((char*)end_token.data(), end_token.size());

		if (end_token != end_token_pat) {
			std::cerr << "Error: wrong serialize end token\n";
			exit(1);
		}
	}

	SuperKmerLoader(const SuperKmerLoader&) = delete;
	SuperKmerLoader& operator=(const SuperKmerLoader&) = delete;

	SuperKmerLoader(SuperKmerLoader&&) noexcept = default;
	SuperKmerLoader& operator=(SuperKmerLoader&&) noexcept = default;

	template<unsigned SIZE>
	void LoadSuperKmer(refresh::bit_memory_reader<refresh::memory_chunk>& in, PackedSuperKmer<SIZE>& super_kmer, uint32_t& sample_id, uint32_t& additional_symbols)
	{
		super_kmer.clear();

		if (!encode_plain)
		{
			uint64_t x = in.read_bits(bits_for_n_additional_symbols + bits_for_sig_id + bits_for_sig_pos);

			uint32_t mmer_pos = x & bits_for_sig_pos_mask;
			x >>= bits_for_sig_pos;
			uint32_t sig_id = x & bits_for_sig_id_mask;
			x >>= bits_for_sig_id;
			additional_symbols = x;

			uint32_t mmer = id_to_mmer[sig_id];

			SuperKmerSerializationHelper::load_super_kmer(in, additional_symbols + kmer_len, mmer_pos, signature_len, mmer, super_kmer);

			sample_id = in.read_bits(bits_for_sample_id);
		}
		else
		{

			additional_symbols = in.read_bits(bits_for_n_additional_symbols);
			auto tot_symbols = kmer_len + additional_symbols;
			auto tot_full_uint64 = tot_symbols / 32;

			size_t i = 0;
			for (; i < tot_full_uint64; ++i)
				super_kmer.set_8bytes(i, in.read_8bytes());

			auto tail_bits = 2 * (tot_symbols - 32 * tot_full_uint64);
			if (tail_bits)
				super_kmer.set_tail(i, in.read_bits(tail_bits), tail_bits);

			sample_id = in.read_bits(bits_for_sample_id);

		}
	}
	uint64_t LoadMetadata(refresh::bit_memory_reader<refresh::memory_chunk>& in, size_t n_bits)
	{
		return in.read_bits(n_bits);
	}
};

