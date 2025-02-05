#ifndef _LZ4_WRAPPER_H
#define _LZ4_WRAPPER_H

#include <cstdint>
#include <cstdio>
#include <string>

#include "../ext_lib/lz4/lz4.h"
#include "../ext_lib/lz4/lz4hc.h"

namespace refresh
{
	// **********************************************************************************
	class lz4_in_memory
	{
		int compression_level;
		bool low_memory;

		void* lz4_state;
		size_t lz4_state_size;

		void prepare_state(int level)
		{
			size_t requested_size;

			if (level < 0)
				requested_size = LZ4_sizeofState();
			else
				requested_size = LZ4_sizeofStateHC();

			if (requested_size > lz4_state_size)
			{
				if(lz4_state)
					free(lz4_state);
				lz4_state_size = requested_size;
				lz4_state = malloc(lz4_state_size);
			}
		}

	public:
		lz4_in_memory(int compression_level = -1, bool low_memory = false) :
			compression_level(compression_level),
			low_memory(low_memory),
			lz4_state(nullptr),
			lz4_state_size(0)
		{}

		lz4_in_memory(lz4_in_memory&& rhs)
		{
			compression_level = rhs.compression_level;
			low_memory = rhs.low_memory;

			lz4_state = rhs.lz4_state;
			rhs.lz4_state = nullptr;
			lz4_state_size = rhs.lz4_state_size;
			rhs.lz4_state_size = 0;
		};

		lz4_in_memory& operator=(lz4_in_memory&& rhs)
		{
			compression_level = rhs.compression_level;
			low_memory = rhs.low_memory;

			if (lz4_state)
				free(lz4_state);

			lz4_state = rhs.lz4_state;
			rhs.lz4_state = nullptr;
			lz4_state_size = rhs.lz4_state_size;
			rhs.lz4_state_size = 0;
		}

		~lz4_in_memory()
		{
			if (lz4_state)
				free(lz4_state);
		}

		void set_compression_level(size_t _compression_level)
		{
			compression_level = _compression_level;
		}

		static int get_min_compression_level()
		{
			return -65537;
		}

		static int get_max_compression_level()
		{
			return LZ4HC_CLEVEL_MAX;
		}

		static size_t get_memory_usage(int level)
		{
			if (level < 0)
				return LZ4_sizeofState();
			else
				return LZ4_sizeofStateHC();
		}

		size_t get_overhead(size_t to_compress_size)	const
		{
			return LZ4_compressBound(to_compress_size) - to_compress_size;
		}

		size_t compress(const void* src, const size_t src_size, void* dest, size_t dest_size, int level = 0)
		{
			if (LZ4_compressBound(src_size) > dest_size)
				return 0;

			if (level == 0)
				level = compression_level;

			if (low_memory)
			{
				if (level < 0)
					return LZ4_compress_fast((char*)src, (char*)dest, src_size, dest_size, -level);
				else
					return LZ4_compress_HC((char*)src, (char*)dest, src_size, dest_size, level);
			}
			else
			{
				prepare_state(level);

				if (level < 0)
					return LZ4_compress_fast_extState(lz4_state, (char*)src, (char*)dest, src_size, dest_size, -level);
				else
					return LZ4_compress_HC_extStateHC(lz4_state, (char*)src, (char*)dest, src_size, dest_size, level);
			}
		}

		size_t decompress(const void* src, const size_t src_size, void* dest, size_t dest_size)
		{
			return LZ4_decompress_safe((char*)src, (char*)dest, src_size, dest_size);
		}
	};	
}

#endif
