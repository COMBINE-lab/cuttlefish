#include "iterate_super_kmers.h"
#include <string>
#include <sstream>

class super_kmer_decoder {
	const unsigned long long* data;
	size_t len;
	size_t pos = 0;
	size_t off = 62;
	size_t words = 0;

	uint8_t get_2bits(const uint32_t p) const {
		return (data[p >> 6] >> (p & 63)) & 3;
	}
public:
	super_kmer_decoder(const unsigned long long* ptr, size_t len, size_t words):
		data(ptr), len(len), words(words)
	{
		
	}
	template<typename SYMBOL_CALLBACLK_T>
	void iterate_symbols(const SYMBOL_CALLBACLK_T& callback) {
		auto p = words * 64 - 2;
		for (uint32_t i = 0; i < len; ++i) {
			callback("ACGT"[get_2bits(p)]);
			p -= 2;
		}
	}
};

int main(int argc, char**argv) {

	if (argc < 4) {
		std::cerr << "Usage: " << argv[0] << " <bins_path> <bin_id> <n_threads>\n";
		exit(1);
	}

	std::string bins_path = argv[1];
	size_t bin_id;
	std::istringstream iss(argv[2]);
	iss >> bin_id;

	size_t n_threads;
	iss = std::istringstream(argv[3]);
	iss >> n_threads;

	IterateSuperKmers iterate(bins_path, bin_id, n_threads);

	//in uint64_t
	size_t data_len = iterate.GetSuperKmerDataLen();

	std::mutex mtx;

	std::vector<std::vector<std::string>> to_flush(n_threads);
	
	size_t size_when_flush = 1000;

	for (size_t tid = 0 ; tid < n_threads ;++tid) {
		iterate.AddConsumer([&data_len, &my_to_flush = to_flush[tid], size_when_flush, &mtx](const unsigned long long* ptr, size_t super_kmer_len_symbols, bool left_flag, bool right_flag) {

			my_to_flush.emplace_back();
			std::string& super_kmer = my_to_flush.back();
			super_kmer.reserve(super_kmer_len_symbols);

			super_kmer_decoder decoder(ptr, super_kmer_len_symbols, data_len);

			decoder.iterate_symbols([&super_kmer](char symb) {
				super_kmer.push_back(symb);
				});

			if (left_flag)
				super_kmer.push_back('l');

			if (right_flag)
				super_kmer.push_back('r');

			if (my_to_flush.size() >= size_when_flush) {
				std::lock_guard<std::mutex> lck(mtx);
				for (const auto& x : my_to_flush)
					std::cout << ">\n" << x << "\n";
				my_to_flush.clear();
			}
		});
	}

	iterate.WaitForAll();

	//flush remaining
	for (size_t thread_id = 0 ; thread_id < n_threads ;++thread_id) {
		for (const auto& x : to_flush[thread_id])
			std::cout << ">\n" << x << "\n";
	}
}
