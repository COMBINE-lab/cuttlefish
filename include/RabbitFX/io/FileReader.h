#include "Buffer.h"
#include "FastxChunk.h"
#include "utils.h"
#include <cstddef>
#include <iostream>
#include <string>
#include <cstring>
#include <memory>
#include <cassert>
#include "Reference.h"

#if defined(USE_IGZIP)
#include "igzip_lib.h"
#endif

#include <zlib.h>  //support gziped files, functional but inefficient

#if defined(_WIN32)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)  // D_SCL_SECURE
#pragma warning(disable : 4244)  // conversion uint64 to uint32
//# pragma warning(disable : 4267)
#define FOPEN fopen
#define FDOPEN fdopen
#define FSEEK _fseeki64
#define FTELL _ftelli64
#define FCLOSE fclose
#elif __APPLE__  // Apple by default suport 64 bit file operations (Darwin 10.5+)
#define FOPEN fopen
#define FDOPEN fdopen
#define FSEEK fseek
#define FTELL ftell
#define FCLOSE fclose
#else
#if !defined(_LARGEFILE_SOURCE)
#define _LARGEFILE_SOURCE
#if !defined(_LARGEFILE64_SOURCE)
#define _LARGEFILE64_SOURCE
#endif
#endif
#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#undef _FILE_OFFSET_BITS
#endif
#if !defined(_FILE_OFFSET_BITS)
#define _FILE_OFFSET_BITS 64
#endif
#define FOPEN fopen64
#define FDOPEN fdopen
#define FSEEK fseeko64
#define FTELL ftello64
#define FCLOSE fclose
#endif


// Forward declarations for `rapidgzip`. Multi-defn errors arise from `rapidgzip` if its headers are included here.
namespace rapidgzip
{

struct ChunkData;
template<typename T_ChunkData> class ParallelGzipReader;

}

namespace rabbit{
	class FileReader{
	private:
		static const uint32 IGZIP_IN_BUF_SIZE = 1 << 22; // 4M gziped file onece fetch
		static const uint32 GZIP_HEADER_BYTES_REQ = 1<<16;
	public:
		FileReader(const std::string &fileName_, bool isZipped, std::size_t worker_count = 1);

		FileReader(int fd, bool isZipped = false);

#if	defined(USE_IGZIP)
		int64 igzip_read(FILE* zipFile, byte *memory_, size_t size_);
#endif

		int64 Read(byte *memory_, uint64 size_);

		/// True means no need to call Read() function
		bool FinishRead();

		bool Eof() const {
			return eof;
			//if(eof) return eof;
			//return feof(mFile);
		}
		void setEof(){
			eof = true;
		}

		~FileReader();

	private:
		FILE* mFile = NULL;
		gzFile mZipFile = NULL;
		//igzip usage
		unsigned char *mIgInbuf = NULL;

#if defined(USE_IGZIP)	
		isal_gzip_header mIgzipHeader;
		inflate_state mStream;

		bool par_deflate;
		std::unique_ptr<rapidgzip::ParallelGzipReader<rapidgzip::ChunkData>> par_gzip_reader;
#endif	
		bool isZipped = false;
		bool eof = false;
	};

} //namespace rabbit
