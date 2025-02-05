#pragma once
#include "defs.h"

#include "libs/refresh/compression/zstd_wrapper.h"

#include <string>
#include <stdexcept>
#include <list>

//mkokot_TODO: decompression could be parallelized like this is done in kb_sorter
//but its a litte tricky so I left it simple now
//for now reading zstd files is not best, because I'm just reading packs of a given compressed size
//probably some reading buffer should be used for better performance...
class BinReader {

	class CompressedPacksWalker {
		std::list<compressed_packs_elem_t> compressed_packs;
		public:
			CompressedPacksWalker(std::list<compressed_packs_elem_t>&& compressed_packs) : compressed_packs(std::move(compressed_packs)) {

			}

			bool Next(size_t& size_compressed, size_t& size_uncompressed, size_t& start_compressed) {
				//shoule be possible only if there is no data in the bin
				if (compressed_packs.size() == 1) {
					assert(compressed_packs.front().start_compressed == 0);
					assert(compressed_packs.front().start_uncompressed == 0);
					compressed_packs.pop_front();
				}

				if (compressed_packs.empty())
					return false;

				start_compressed = compressed_packs.front().start_compressed;
				size_t start_uncompressed = compressed_packs.front().start_uncompressed;

				compressed_packs.pop_front();

				size_compressed = compressed_packs.front().start_compressed - start_compressed;
				size_uncompressed = compressed_packs.front().start_uncompressed - start_uncompressed;

				//remove guard
				if (compressed_packs.size() == 1)
					compressed_packs.pop_front();

				return true;
			}
	};

	FILE* file{};
	bool is_zstd_compressed;

	bin_reader_super_kmers_packs_maker_queue_t& out_q;
	CompressedPacksWalker compressed_walker;
public:
	BinReader(const std::string& path, bool is_zstd_compressed, bin_reader_super_kmers_packs_maker_queue_t& out_q, std::list<compressed_packs_elem_t>& compressed_packs):
		is_zstd_compressed(is_zstd_compressed),
		out_q(out_q),
		compressed_walker(std::move(compressed_packs)) {

		file = fopen(path.c_str(), "rb");
		if (!file) {
			std::cerr << "Error: Cannot open file " << path << "\n";
			exit(1);
		}
	}

	void Process() {
		if (is_zstd_compressed) {

			size_t size_compressed;
			size_t size_uncompressed;
			size_t start_compressed;
			reader_data_t in;
			while (compressed_walker.Next(size_compressed, size_uncompressed, start_compressed)) {
				fseek(file, start_compressed, SEEK_SET);

				if (in.data.capacity() < size_compressed)
					in = reader_data_t(size_compressed);

				in.data.resize(size_compressed);
				if (size_compressed != fread(in.data.data(), 1, size_compressed, file)) {
					std::cerr <<"Error: something went wrong reading from zstd compressed bin\n";
					exit(1);
				}

				reader_data_t out(size_uncompressed);
				out.data.resize(size_uncompressed);

				refresh::zstd_in_memory zstd;
				if (out.data.size() != zstd.decompress(in.data.data(), in.data.size(), out.data.data(), out.data.size())) {
					std::cerr << "Error: something went wrong decompressing zstd pack\n";
					exit(1);
				}

				out_q.push(std::move(out));
			}
			out_q.mark_completed();
		}
		else {
			while (true) {
				reader_data_t _out;
				auto& out = _out.data;
				auto readed = fread(out.data(), 1, out.capacity(), file);
				out.resize(readed);
				if (!readed) {
					fclose(file);
					out_q.mark_completed();
					return;
				}
				out_q.push(std::move(_out));
			}
		}
	}
};
