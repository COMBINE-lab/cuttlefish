
#include "FileReader.h"

#include "rapidgzip/ParallelGzipReader.hpp"


namespace rabbit
{

FileReader::FileReader(const std::string &fileName_, bool isZipped, const std::size_t worker_count){
    if(ends_with(fileName_, ".gz") || isZipped) {

        this->isZipped = true;
        par_deflate = (worker_count > 1);

#if defined(USE_IGZIP)
        if(par_deflate)
        {
            auto file_reader = std::make_unique<StandardFileReader>(fileName_);
            par_gzip_reader = std::make_unique<rapidgzip::ParallelGzipReader<>>(std::move(file_reader), worker_count);
            par_gzip_reader->setCRC32Enabled(false);

            return;
        }

        // printf("using igzip!!\n");
        mFile = fopen(fileName_.c_str(), "rb");
        if (mFile == NULL){
            throw RioException(
                ("Can not open file to read: " + fileName_).c_str());  
        }
        mIgInbuf = new unsigned char[IGZIP_IN_BUF_SIZE];	// TODO: add some `init` / `reset` method to reuse this mem.
        isal_gzip_header_init(&mIgzipHeader);
        isal_inflate_init(&mStream);
        mStream.crc_flag = ISAL_GZIP_NO_HDR_VER;

        mStream.next_in = mIgInbuf;
        mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, mFile);

        int ret =0;
        ret = isal_read_gzip_header(&mStream, &mIgzipHeader);
        if( ret != ISAL_DECOMP_OK ){
            std::cerr << "error invalid gzip header found: " << fileName_ << std::endl;
            if(mFile != NULL){
                fclose(mFile);
            }
            exit(-1);
        }
#else
        mZipFile = gzopen(fileName_.c_str(), "r");
        if (mZipFile == NULL) {
            throw RioException(("Can not open file to read: "+fileName_).c_str());  
        }
        gzrewind(mZipFile);
#endif			
    }else {
        // mFile = FOPEN(fileName_.c_str(), "rb");
        if (fileName_ != "") {
            mFile = FOPEN(fileName_.c_str(), "rb");
            if (mFile == NULL)
                throw RioException(("Can not open file to read: " + fileName_).c_str());  
        }
        if (mFile == NULL) {
            throw RioException(
                ("Can not open file to read: " + fileName_).c_str());  
        }
    }
}


/*
FileReader::FileReader(int fd, bool isZipped){
    if (isZipped) {
#if defined(USE_IGZIP)
        mFile = fdopen(fd, "r");
        if (mFile == NULL)
        {
            std::cerr << "Can not open file id to read: " + std::to_string(fd) << std::endl;
        }
        mIgInbuf = new unsigned char[IGZIP_IN_BUF_SIZE];
        isal_gzip_header_init(&mIgzipHeader);
        isal_inflate_init(&mStream);
        mStream.crc_flag = ISAL_GZIP_NO_HDR_VER;

        mStream.next_in = mIgInbuf;
        mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, mFile);

        int ret = 0;
        ret = isal_read_gzip_header(&mStream, &mIgzipHeader);
        if (ret != ISAL_DECOMP_OK)
        {
            std::cerr << "error! invalid gzip header found!" << std::endl;
            if (mFile != NULL)
            {
                fclose(mFile);
            }
            exit(-1);
        }
#else
        mZipFile = gzdopen(fd, "r");
        gzrewind(mZipFile);
#endif
        this->isZipped = true;
    } else {
        mFile = FDOPEN(fd, "rb");
        if (fd != -1) {
            mFile = FDOPEN(fd, "rb");
            if (mFile == NULL)
                throw RioException("Can not open file to read: ");  
        }
        if (mFile == NULL) {
            throw RioException("Can not open file to read: ");  
        }
    }
}
*/


#if	defined(USE_IGZIP)
__attribute__((no_sanitize("memory")))
int64 FileReader::igzip_read(FILE* zipFile, byte *memory_, size_t size_){
    uint64_t offset = 0;
    int ret = 0;
    do{
        if(mStream.avail_in == 0){
            mStream.next_in = mIgInbuf;
            mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, zipFile);
        }
        mStream.next_out = memory_ + offset;
        mStream.avail_out = size_ - offset;
        if(isal_inflate(&mStream) != ISAL_DECOMP_OK){
            std::cerr << "decompress error" << std::endl;
            return -1;
        }
        offset = (mStream.next_out - memory_);

        if(mStream.block_state == ISAL_BLOCK_FINISH) {
            if(feof(mFile) && mStream.avail_in == 0){
                return offset;
            }
            // a new block begins
            if (mStream.avail_in == 0) {
                isal_inflate_reset(&mStream);
                mStream.next_in = mIgInbuf;
                mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, mFile);
                //mGzipInputUsedBytes += mStream.avail_in;
            } else if (mStream.avail_in >= GZIP_HEADER_BYTES_REQ){
                unsigned char* old_next_in = mStream.next_in;
                size_t old_avail_in = mStream.avail_in;
                isal_inflate_reset(&mStream);
                mStream.avail_in = old_avail_in;
                mStream.next_in = old_next_in;
            } else {
                size_t old_avail_in = mStream.avail_in;
                memmove(mIgInbuf, mStream.next_in, mStream.avail_in);
                size_t readed = 0;
                if(!feof(mFile)){
                    readed = fread(mIgInbuf + mStream.avail_in, 1, IGZIP_IN_BUF_SIZE - mStream.avail_in, mFile);
                }
                isal_inflate_reset(&mStream);
                mStream.next_in = mIgInbuf;
                mStream.avail_in = old_avail_in + readed;;
            }
            ret = isal_read_gzip_header(&mStream, &mIgzipHeader);
            if (ret != ISAL_DECOMP_OK) {
                std::cerr << "igzip: invalid gzip header found: " << mStream.avail_in <<  " : " << mStream.avail_out << std::endl;
                exit(-1);
            }
        }
    }while(mStream.avail_out > 0);
    assert(offset <= size_);
    return offset;
}
#endif


int64 FileReader::Read(byte *memory_, uint64 size_) {
    if (isZipped) {
#if defined(USE_IGZIP)			
        if(par_deflate)
        {
            const auto r = par_gzip_reader->read(reinterpret_cast<char*>(memory_), size_);
            if(r < size_)
                setEof();

            return r;
        }

        int64 n = igzip_read(mFile, memory_, size_);
#else
        int64 n = gzread(mZipFile, memory_, size_);
#endif
        if (n == -1) std::cerr << "Error to read gzip file" << std::endl;
        return n;
    } else {
        int64 n = fread(memory_, 1, size_, mFile);
        return n;
    }
}


FileReader::~FileReader(){
    if(mIgInbuf != NULL) delete[] mIgInbuf;
    if(mFile != NULL){
        fclose(mFile);
        mFile = NULL;
    }
    if(mZipFile != NULL){
        gzclose(mZipFile);
    }
}


bool FileReader::FinishRead(){
    if(isZipped){
#if defined(USE_IGZIP)			
        if(par_deflate)
            return Eof();

        return feof(mFile) && (mStream.avail_in == 0);
#else
        return gzeof(mZipFile);
#endif			
    }else{
        return feof(mFile);
    }
}
}
