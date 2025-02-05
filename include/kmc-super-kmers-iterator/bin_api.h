#pragma once
#include "libs/refresh/serialization/serialization.h"
#include <fstream>
#include <list>

struct BinsGlobalConfig {
	size_t k{};
	bool is_zstd_compr{};
	size_t n_bins{};
	bool super_kmers_with_overlaps{};
public:
	void Serialize(std::ofstream& out) const {
		refresh::serialize_little_endian(k, out);
		refresh::serialize_little_endian(is_zstd_compr, out);
		refresh::serialize_little_endian(n_bins, out);
		refresh::serialize_little_endian(super_kmers_with_overlaps, out);
	}

	void Load(std::ifstream& in) {
		refresh::load_little_endian(k, in);
		refresh::load_little_endian(is_zstd_compr, in);
		refresh::load_little_endian(n_bins, in);
		refresh::load_little_endian(super_kmers_with_overlaps, in);
	}
};

struct simple_pack_data_t {
	uint64_t end_pos{};
	uint64_t n_super_kmers{};

	simple_pack_data_t(uint64_t end_pos, uint64_t n_super_kmers) :end_pos(end_pos), n_super_kmers(n_super_kmers) {}

	void Serialize(std::ostream& out) const {
		refresh::serialize_little_endian(end_pos, out);
		refresh::serialize_little_endian(n_super_kmers, out);
	}

	void Load(std::istream& in) {
		refresh::load_little_endian(end_pos, in);
		refresh::load_little_endian(n_super_kmers, in);
	}
};

struct compressed_packs_elem_t {
	size_t start_uncompressed;
	size_t start_compressed;
	compressed_packs_elem_t(size_t start_uncompressed, size_t start_compressed) :
		start_uncompressed(start_uncompressed),
		start_compressed(start_compressed)
	{

	}
	void Serialize(std::ostream& out) const {
		refresh::serialize_little_endian(start_uncompressed, out);
		refresh::serialize_little_endian(start_compressed, out);
	}

	void Load(std::istream& in) {
		refresh::load_little_endian(start_uncompressed, in);
		refresh::load_little_endian(start_compressed, in);
	}
};

struct expand_packs_t {
	std::list<simple_pack_data_t> data;
	void Serialize(std::ostream& out) const {
		uint64_t size = data.size();
		refresh::serialize_little_endian(size, out);
		for (auto& x : data)
			x.Serialize(out);
	}
	void Load(std::istream& in) {
		uint64_t size;
		refresh::load_little_endian(size, in);
		simple_pack_data_t p(0, 0);
		for (size_t i = 0; i < size; ++i) {
			p.Load(in);
			data.emplace_back(p);
		}
	}
};

struct compressed_packs_t {
	std::list<compressed_packs_elem_t> data;
	void Serialize(std::ostream& out) const {
		uint64_t size = data.size();
		refresh::serialize_little_endian(size, out);
		for (auto& x : data)
			x.Serialize(out);
	}
	void Load(std::istream& in) {
		uint64_t size;
		refresh::load_little_endian(size, in);
		compressed_packs_elem_t p(0, 0);
		for (size_t i = 0; i < size; ++i) {
			p.Load(in);
			data.emplace_back(p);
		}
	}
};

struct bin_meta_t {
	expand_packs_t expand_packs;
	compressed_packs_t compressed_packs;
	void Serialize(std::ostream& out) const {
		expand_packs.Serialize(out);
		compressed_packs.Serialize(out);
	}
	void Load(std::istream& in) {
		expand_packs.Load(in);
		compressed_packs.Load(in);
	}
};
