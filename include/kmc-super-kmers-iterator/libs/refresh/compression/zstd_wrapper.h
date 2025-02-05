#ifndef _ZSTD_WRAPPER_H
#define _ZSTD_WRAPPER_H

#include <cstdint>
#include <cstdio>
#include <string>

#include "ext_lib/zstd/zstd.h"

namespace refresh
{
	// **********************************************************************************
	class zstd_in_memory
	{
		int compression_level;
		bool low_memory;

		ZSTD_CCtx* zstd_cctx;
		ZSTD_DCtx* zstd_dctx;
		ZSTD_CStream* zstd_cstream;

		enum class mode_t {none, buffer, stream};

		char* stream_data;
		size_t stream_capacity;
		size_t stream_in_use;

		mode_t working_mode;

	public:
		zstd_in_memory(int compression_level = 9, bool low_memory = false) :
			compression_level(compression_level),
			low_memory(low_memory),
			zstd_cctx(nullptr),
			zstd_dctx(nullptr),
			zstd_cstream(nullptr),
			stream_data(nullptr),
			stream_capacity(0),
			stream_in_use(0),
			working_mode(mode_t::none)
		{}

		zstd_in_memory(zstd_in_memory&& rhs) noexcept
		{
			compression_level = rhs.compression_level;
			low_memory = rhs.low_memory;

			working_mode = rhs.working_mode;
			rhs.working_mode = mode_t::none;

			zstd_cctx = rhs.zstd_cctx;
			rhs.zstd_cctx = nullptr;

			zstd_dctx = rhs.zstd_dctx;
			rhs.zstd_dctx = nullptr;

			zstd_cstream = rhs.zstd_cstream;
			rhs.zstd_cstream = nullptr;

			stream_data = rhs.stream_data;
			rhs.stream_data = nullptr;
			
			stream_capacity = rhs.stream_capacity;
			rhs.stream_capacity = 0;
			
			stream_in_use = rhs.stream_in_use;
			rhs.stream_in_use = 0;
		};

		zstd_in_memory& operator=(zstd_in_memory&& rhs) noexcept
		{
			compression_level = rhs.compression_level;
			low_memory = rhs.low_memory;
			working_mode = rhs.working_mode;

			rhs.working_mode = mode_t::none;

			if (zstd_cctx)
				ZSTD_freeCCtx(zstd_cctx);

			zstd_cctx = rhs.zstd_cctx;
			rhs.zstd_cctx = nullptr;

			if (zstd_dctx)
				ZSTD_freeDCtx(zstd_dctx);

			zstd_dctx = rhs.zstd_dctx;
			rhs.zstd_dctx = nullptr;

			zstd_cstream = rhs.zstd_cstream;
			rhs.zstd_cstream = nullptr;

			stream_data = rhs.stream_data;
			rhs.stream_data = nullptr;

			stream_capacity = rhs.stream_capacity;
			rhs.stream_capacity = 0;

			stream_in_use = rhs.stream_in_use;
			rhs.stream_in_use = 0;

			return *this;
		}

		~zstd_in_memory() noexcept
		{
			if (zstd_cctx)
				ZSTD_freeCCtx(zstd_cctx);

			if (zstd_dctx)
				ZSTD_freeDCtx(zstd_dctx);

			if (zstd_cstream)
				ZSTD_freeCStream(zstd_cstream);
		}

		void set_compression_level(size_t _compression_level)
		{
			compression_level = _compression_level;
		}

		static int get_min_compression_level()
		{
			return 1;
		}

		static int get_max_compression_level()
		{
			return 19;
		}

/*		static size_t get_memory_usage(int level)
		{
			// Internal function - static linking only
			return ZSTD_estimateCCtxSize(level);
		}
*/
		static size_t get_overhead(size_t to_compress_size)
		{
			return ZSTD_compressBound(to_compress_size) - to_compress_size;
		}

		size_t compress(const void* src, const size_t src_size, void* dest, size_t dest_size, int level = 0)
		{
			if (working_mode == mode_t::none)
				working_mode = mode_t::buffer;
			else if (working_mode != mode_t::buffer)
				return 0;

			if (ZSTD_compressBound(src_size) > dest_size)
				return 0;

			if (level == 0)
				level = compression_level;

			if (!zstd_cctx)
				zstd_cctx = ZSTD_createCCtx();

			auto r = ZSTD_compressCCtx(zstd_cctx, dest, dest_size, src, src_size, level);

			if (low_memory)
			{
				ZSTD_freeCCtx(zstd_cctx);
				zstd_cctx = nullptr;
			}

			return r;
		}
		
		bool init_stream(void *dest, size_t dest_size, int level = 0)
		{
			if (working_mode != mode_t::none)
				return false;

			working_mode = mode_t::stream;

			if (level != 0)
				compression_level = level;

			working_mode = mode_t::stream;
			stream_data = (char*) dest;
			stream_capacity = dest_size;
			stream_in_use = 0;

			if (!zstd_cstream)
				zstd_cstream = ZSTD_createCStream();

			ZSTD_initCStream(zstd_cstream, compression_level);

			return true;
		}

		bool add_to_stream(const void* src, const size_t src_size)
		{
			if (working_mode != mode_t::stream)
				return false;

			ZSTD_inBuffer zstd_in_buffer;
			ZSTD_outBuffer zstd_out_buffer;

			zstd_in_buffer.src = src;
			zstd_in_buffer.pos = 0;
			zstd_in_buffer.size = src_size;

			zstd_out_buffer.dst = stream_data;
			zstd_out_buffer.pos = stream_in_use;
			zstd_out_buffer.size = stream_capacity;

			ZSTD_compressStream2(zstd_cstream, &zstd_out_buffer, &zstd_in_buffer, ZSTD_e_continue);
			
			stream_in_use = zstd_out_buffer.pos;

			if (zstd_in_buffer.pos < zstd_in_buffer.size)
				return false;		// Not enough space in out buffer

			return true;
		}

		size_t close_stream()
		{
			if (working_mode != mode_t::stream)
				return 0;

			ZSTD_inBuffer zstd_in_buffer;
			ZSTD_outBuffer zstd_out_buffer;

			char aux;
			zstd_in_buffer.src = &aux;
			zstd_in_buffer.pos = 0;
			zstd_in_buffer.size = 0;

			zstd_out_buffer.dst = stream_data;
			zstd_out_buffer.pos = stream_in_use;
			zstd_out_buffer.size = stream_capacity;

			bool all_sent = ZSTD_compressStream2(zstd_cstream, &zstd_out_buffer, &zstd_in_buffer, ZSTD_e_end) == 0;

			stream_in_use = zstd_out_buffer.pos;

			if (low_memory)
			{
				ZSTD_freeCStream(zstd_cstream);
				zstd_cstream = nullptr;
			}

			working_mode = mode_t::none;

			if (!all_sent)
				return 0;		// Not enough space in output buffer

			return stream_in_use;
		}

		size_t decompress(const void* src, const size_t src_size, void* dest, size_t dest_size)
		{
			if (!zstd_dctx)
				zstd_dctx = ZSTD_createDCtx();

			auto r = ZSTD_decompressDCtx(zstd_dctx, dest, dest_size, src, src_size);

			if (low_memory)
			{
				ZSTD_freeDCtx(zstd_dctx);
				zstd_dctx = nullptr;
			}

			return r;
		}
	};

	// **********************************************************************************
	class zstd_file
	{
		enum class working_mode_t {none, reading, writing};

		working_mode_t working_mode;
		FILE* fio;
		int compression_level;
		size_t io_buffer_size;
		size_t zstd_buffer_size;

		char *mem_buf[2];

		ZSTD_CStream *zstd_cstream;
		ZSTD_DStream *zstd_dstream;

		ZSTD_inBuffer zstd_read_buffer;
		bool read_eof;

		void flush()
		{
			ZSTD_inBuffer zstd_in_buffer;
			ZSTD_outBuffer zstd_out_buffer;

			zstd_in_buffer.src = mem_buf[0];
			zstd_in_buffer.pos = 0;
			zstd_in_buffer.size = 0;

			zstd_out_buffer.dst = mem_buf[1];
			zstd_out_buffer.pos = 0;
			zstd_out_buffer.size = zstd_buffer_size;

			for(bool need_work = true; need_work;)
			{
				zstd_out_buffer.pos = 0;
				need_work = ZSTD_compressStream2(zstd_cstream, &zstd_out_buffer, &zstd_in_buffer, ZSTD_e_end) != 0;
				fwrite(zstd_out_buffer.dst, 1, zstd_out_buffer.pos, fio);
			}
		}

		void close_reading()
		{
			if (!fio)
				return;
			fclose(fio);
			fio = nullptr;

			ZSTD_freeDStream(zstd_dstream);
			zstd_dstream = nullptr;

			working_mode = working_mode_t::none;
		}

		void close_writing()
		{
			if (!fio)
				return;

			flush();

			fclose(fio);
			fio = nullptr;

			ZSTD_freeCStream(zstd_cstream);
			zstd_cstream = nullptr;

			working_mode = working_mode_t::none;
		}

		void allocate_buffers()
		{
			release_buffers();

			mem_buf[0] = new char[zstd_buffer_size];
			mem_buf[1] = new char[zstd_buffer_size];
		}

		void release_buffers()
		{
			if (mem_buf[0])
			{
				delete[] mem_buf[0];
				mem_buf[0] = nullptr;
			}

			if (mem_buf[1])
			{
				delete[] mem_buf[1];
				mem_buf[1] = nullptr;
			}
		}

	public:
		zstd_file(int compression_level = 9, size_t io_buffer_size = 1 << 20, size_t zstd_buffer_size = 1 << 20) : 
			working_mode(working_mode_t::none), 
			fio(nullptr), 
			compression_level(compression_level),
			io_buffer_size(io_buffer_size),
			zstd_buffer_size(zstd_buffer_size)
		{
			mem_buf[0] = nullptr;
			mem_buf[1] = nullptr;

			zstd_cstream = nullptr;
			zstd_dstream = nullptr;
		}

		zstd_file(zstd_file&& rhs)
		{
			working_mode = rhs.working_mode;
			rhs.working_mode = working_mode_t::none;

			fio = rhs.fio;
			rhs.fio = nullptr;

			compression_level = rhs.compression_level;

			io_buffer_size = rhs.io_buffer_size;
			zstd_buffer_size = rhs.zstd_buffer_size;

			mem_buf[0] = rhs.mem_buf[0];
			rhs.mem_buf[0] = nullptr;

			mem_buf[1] = rhs.mem_buf[1];
			rhs.mem_buf[1] = nullptr;

			zstd_cstream = rhs.zstd_cstream;
			rhs.zstd_cstream = nullptr;

			zstd_dstream = rhs.zstd_dstream;
			rhs.zstd_dstream = nullptr;

			zstd_read_buffer = rhs.zstd_read_buffer;
			rhs.zstd_read_buffer.src = nullptr;

			read_eof = rhs.read_eof;
		};

		zstd_file& operator=(zstd_file&& rhs)
		{
			if (fio)
				close();

			working_mode = rhs.working_mode;
			rhs.working_mode = working_mode_t::none;

			fio = rhs.fio;
			rhs.fio = nullptr;

			compression_level = rhs.compression_level;

			io_buffer_size = rhs.io_buffer_size;
			zstd_buffer_size = rhs.zstd_buffer_size;

			if (mem_buf[0])
				delete[] mem_buf[0];

			mem_buf[0] = rhs.mem_buf[0];
			rhs.mem_buf[0] = nullptr;

			if (mem_buf[1])
				delete[] mem_buf[1];

			mem_buf[1] = rhs.mem_buf[1];
			rhs.mem_buf[1] = nullptr;

			if (zstd_cstream)
				ZSTD_freeCStream(zstd_cstream);

			zstd_cstream = rhs.zstd_cstream;
			rhs.zstd_cstream = nullptr;

			if (zstd_dstream)
				ZSTD_freeDStream(zstd_dstream);

			zstd_dstream = rhs.zstd_dstream;
			rhs.zstd_dstream = nullptr;

			zstd_read_buffer = rhs.zstd_read_buffer;
			rhs.zstd_read_buffer.src = nullptr;

			read_eof = rhs.read_eof;

			return *this;
		};

		zstd_file(const zstd_file& rhs) = delete;

		~zstd_file()
		{
			close();

			release_buffers();
		}

		void close()
		{
			if (working_mode == working_mode_t::reading)
				close_reading();
			else if (working_mode == working_mode_t::writing)
				close_writing();

			release_buffers();
		}

		void set_compression_level(size_t _compression_level)
		{
			if(working_mode == working_mode_t::none)
				compression_level = _compression_level;
		}

		static int get_min_compression_level()
		{
			return 1;
		}

		static int get_max_compression_level()
		{
			return 19;
		}

		void set_io_buffer_size(size_t _io_buffer_size)
		{
			if (working_mode == working_mode_t::none)
				io_buffer_size = _io_buffer_size;
		}

		void set_zstd_buffer_size(size_t _zstd_buffer_size)
		{
			if (working_mode == working_mode_t::none)
			{
				zstd_buffer_size = _zstd_buffer_size;

				allocate_buffers();
			}
		}

		bool eof()
		{
			return read_eof;
		}

		bool is_opened()
		{
			return working_mode != working_mode_t::none;
		}

		bool is_opened_for_reading()
		{
			return working_mode == working_mode_t::reading;
		}
		
		bool is_opened_for_writing()
		{
			return working_mode == working_mode_t::writing;
		}

		bool open_reading(const std::string& file_name)
		{
			if (working_mode != working_mode_t::none)
				return false;

			fio = fopen(file_name.c_str(), "rb");
			if (!fio)
				return false;

			setvbuf(fio, nullptr, _IOFBF, io_buffer_size);

			working_mode = working_mode_t::reading;

			allocate_buffers();

			zstd_dstream = ZSTD_createDStream();

			zstd_read_buffer.src = mem_buf[0];
			zstd_read_buffer.pos = 0;
			zstd_read_buffer.size = 0;

			read_eof = false;

			ZSTD_initDStream(zstd_dstream);

			return true;
		}

		bool open_writing(const std::string& file_name)
		{
			if (working_mode != working_mode_t::none)
				return false;

			fio = fopen(file_name.c_str(), "wb");
			if (!fio)
				return false;

			setvbuf(fio, nullptr, _IOFBF, io_buffer_size);

			working_mode = working_mode_t::writing;

			allocate_buffers();

			zstd_cstream = ZSTD_createCStream();

			ZSTD_initCStream(zstd_cstream, compression_level);

			return true;
		}

		bool put(char c)
		{
			return write(&c, 1);
		}

		bool write(char* p, size_t size)
		{
			if (working_mode != working_mode_t::writing)
				return false;

			ZSTD_inBuffer zstd_in_buffer;
			ZSTD_outBuffer zstd_out_buffer;

			zstd_in_buffer.src = p;
			zstd_in_buffer.pos = 0;
			zstd_in_buffer.size = size;

			zstd_out_buffer.dst = mem_buf[1];
			zstd_out_buffer.pos = 0;
			zstd_out_buffer.size = zstd_buffer_size;

			while (zstd_in_buffer.pos < zstd_in_buffer.size)
			{
				ZSTD_compressStream2(zstd_cstream, &zstd_out_buffer, &zstd_in_buffer, ZSTD_e_continue);

				fwrite(zstd_out_buffer.dst, 1, zstd_out_buffer.pos, fio);
				zstd_out_buffer.pos = 0;
			}

			return true;
		}

		size_t get(char& c)
		{
			return read(&c, 1);
		}

		size_t read(char* p, size_t size)
		{
			if (working_mode != working_mode_t::reading || read_eof)
				return 0;

			size_t readed = 0;

			ZSTD_outBuffer zstd_out_buffer;

			while (readed < size && !read_eof)
			{
				if (zstd_read_buffer.pos == zstd_read_buffer.size)
				{
					zstd_read_buffer.pos = 0;
					zstd_read_buffer.size = fread((void*)zstd_read_buffer.src, 1, zstd_buffer_size, fio);
				}

				zstd_out_buffer.dst = p + readed;
				zstd_out_buffer.pos = 0;
				zstd_out_buffer.size = size - readed;

				read_eof = ZSTD_decompressStream(zstd_dstream, &zstd_out_buffer, &zstd_read_buffer) == 0;

				readed += zstd_out_buffer.pos;
			}

			return readed;
		}
	};
}

#endif
