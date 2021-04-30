
#ifndef ASYNC_LOGGER_WRAPPER_HPP
#define ASYNC_LOGGER_WRAPPER_HPP



#include "spdlog/async_logger.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>


// A class wrapping the `spdlog` library's asynchronous logger.
class Async_Logger_Wrapper
{
private:

    // `spdlog`'s queue size, i.e. the maximum number of log message units it can contain before a flush to sink.
    static constexpr std::size_t QUEUE_CAP = 1024;

    // Number of backing worker threads for `spdlog`, i.e. the threads that actually make the writes to sink.
    static constexpr uint16_t NUM_THREADS = 1;

    // `spdlog` thread pool for performing the outputting task to sink. Unless multiple distinct sinks are to be
    // present (e.g. in the writing algorithm for the GFA-variants in reference dBG compaction), only one thread
    // pool is needed, and it does not make much sense in having multiple backing threads in that pool. And through
    // having a dedicated thread pool for the actual (aynchronous) sink-flushes, the disk-write happens in parallel
    // to the algorithm operation.
    std::shared_ptr<spdlog::details::thread_pool> tp;

    // Output logger.
    std::shared_ptr<spdlog::logger> logger;


public:

    // Initializes the `spdlog` logger wrapper that writes to a file with path `output_file_path`.
    void init_logger(const std::string& output_file_path);

    // Log the passed null-terminated message `str`.
    void write(const char* str) const;
};


inline void Async_Logger_Wrapper::write(const char* const str) const
{
    logger->info(str);
}



#endif
