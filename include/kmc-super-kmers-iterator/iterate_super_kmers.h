#pragma once
#include "defs.h"
#include <vector>
#include <thread>
#include <string>
#include "super_kmer_loader.h"
#include "bin_api.h"
#include "bin_reader.h"
#include "super_kmers_packs_maker.h"

template<typename CALLBACK_T>
class IIterateSuperKmers {
public:
	virtual void Run(CALLBACK_T&& callback, uint32_t k, super_kmers_packs_maker_super_kmer_iterator_queue_t& in_queue, SuperKmerLoader* super_kmer_serializer, bool super_kmers_with_overlaps) = 0;
	~IIterateSuperKmers() = default;
};

template<typename CALLBACK_T, unsigned KMER_SIZE>
class IterateSuperKmersImpl : public IIterateSuperKmers<CALLBACK_T>
{
public:
	void Run(CALLBACK_T&& callback, uint32_t k, super_kmers_packs_maker_super_kmer_iterator_queue_t& in_queue, SuperKmerLoader* super_kmer_serializer, bool super_kmers_with_overlaps) override {
		while (in_queue.pop_and_consume([&](super_kmers_packer_data_t&& data) {
			//std::cerr << "pop\n";
			uint32_t additional_symbols;
			uint64_t pack_super_kmers = 0;
			uint32_t sample_id;
			uint64_t offset = 0;

			PackedSuperKmer<KMER_SIZE> super_kmer;

			for (size_t i = 0; i < data.n_super_kmers; i += pack_super_kmers) {
				auto ptr = reinterpret_cast<uint64_t*>(data.data.data() + offset);

				uint64_t size = *ptr++; //of pack, in bytes
				pack_super_kmers = *ptr++;

				offset += size;
				assert(size % sizeof(uint64_t) == 0);
				size /= sizeof(uint64_t); //now in uint64

				refresh::memory_chunk<uint64_t> _in(ptr, size);
				_in.resize(size);

				refresh::bit_memory_reader<refresh::memory_chunk> in(std::move(_in));
				//std::cerr << "pack_super_kmers: " << pack_super_kmers << "\n";
				for (size_t j = 0; j < pack_super_kmers; ++j)
				{
					super_kmer_serializer->LoadSuperKmer(in, super_kmer, sample_id, additional_symbols);

					bool left_flag = false;
					bool right_flag = false;
					if (super_kmers_with_overlaps)
					{
						left_flag = super_kmer_serializer->LoadMetadata(in, 1);
						right_flag = super_kmer_serializer->LoadMetadata(in, 1);
					}
					callback(super_kmer.raw_data(), k + additional_symbols, left_flag, right_flag);
				}
			}
		}));
	}
};

template<typename CALLBACK_T, unsigned KMER_SIZE>
struct IIterateSuperKmersDispatch {
	static std::unique_ptr<IIterateSuperKmers<CALLBACK_T>> Get(uint32_t k) {
		if (k > 32 * KMER_SIZE) {
			std::cerr << "Please increase MAX_KMER_SIZE in defs.h\n";
			exit(1);
		}

		if (k <= 32 * KMER_SIZE && k > 32 * (KMER_SIZE - 1))
			return std::make_unique<IterateSuperKmersImpl<CALLBACK_T, KMER_SIZE>>();
		return IIterateSuperKmersDispatch<CALLBACK_T, KMER_SIZE - 1>::Get(k);
	}
};

template<typename CALLBACK_T>
struct IIterateSuperKmersDispatch<CALLBACK_T, 0> {
	static std::unique_ptr<IIterateSuperKmers<CALLBACK_T>> Get(uint32_t k) {
		(void)k;
		std::cerr << "Error: This should never happen\n";
		exit(1);
	}
};

class IterateSuperKmers {
	std::vector<std::thread> threads;

	bin_reader_super_kmers_packs_maker_queue_t bin_reader_super_kmers_packs_maker_queue;
	super_kmers_packs_maker_super_kmer_iterator_queue_t super_kmers_packs_maker_super_kmer_iterator_queue;

	BinsGlobalConfig bins_global_config;
	bin_meta_t bin_meta;

	std::unique_ptr<SuperKmerLoader> super_kmers_serializer;
	
	std::string get_bin_path(const std::string& bins_path, size_t bin_id) {
		std::string s_tmp = std::to_string(bin_id);
		while (s_tmp.length() < 5)
			s_tmp = std::string("0") + s_tmp;

		return bins_path + "kmc_" + s_tmp + ".bin";
	}

	void read_global_config(const std::string& path) {
		std::ifstream in(path, std::ios::binary);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		bins_global_config.Load(in);
	}

	void read_bin_meta(const std::string& path) {
		std::ifstream in(path, std::ios::binary);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		bin_meta.Load(in);

		super_kmers_serializer = std::make_unique<SuperKmerLoader>(in);

	}
public:
	//bins_path is a tmp path of kmc run
	IterateSuperKmers(std::string bins_path, size_t bin_id, size_t queue_size) :
		bin_reader_super_kmers_packs_maker_queue(2), //mkokot_TODO: size OK?
		super_kmers_packs_maker_super_kmer_iterator_queue(queue_size + 1)
	{	
		if (bins_path.back() != '/' && bins_path.back() != '\\')
			bins_path.push_back('/');

		read_global_config(bins_path + "bins.global");

		auto bin_path = get_bin_path(bins_path, bin_id);
		auto bin_path_meta = bin_path + ".meta";

		read_bin_meta(bin_path_meta);

		threads.emplace_back([&, bin_path]() {
			BinReader bin_reader(bin_path, bins_global_config.is_zstd_compr, bin_reader_super_kmers_packs_maker_queue, bin_meta.compressed_packs.data);
			bin_reader.Process();
		});

		threads.emplace_back([&]() {
			SuperKmersPacksMaker super_kmers_pack_maker(bin_reader_super_kmers_packs_maker_queue, super_kmers_packs_maker_super_kmer_iterator_queue, bin_meta.expand_packs.data);
			super_kmers_pack_maker.Process();
		});
	}

	size_t GetSuperKmerDataLen() {
		auto k = bins_global_config.k;
		auto kmer_words = (k + 31) / 32;
		auto super_kmer_words = kmer_words * 2;
		return super_kmer_words;
	}
	template<typename SUPER_KMER_CALLBACK_T>
	void AddConsumer(SUPER_KMER_CALLBACK_T&& super_kmer_callback) {
		threads.emplace_back([this] (SUPER_KMER_CALLBACK_T&& super_kmer_callback) {
			auto impl = IIterateSuperKmersDispatch<SUPER_KMER_CALLBACK_T, MAX_KMER_SIZE>::Get(bins_global_config.k);
			impl->Run(std::move(super_kmer_callback), bins_global_config.k, super_kmers_packs_maker_super_kmer_iterator_queue, super_kmers_serializer.get(), bins_global_config.super_kmers_with_overlaps);
		}, std::move(super_kmer_callback));
	}
	void WaitForAll() {
		for (auto& thread : threads)
			thread.join();
	}
};
