
#include "Async_Logger_Wrapper.hpp"

#include "spdlog/sinks/basic_file_sink.h"


void Async_Logger_Wrapper::init_logger(const std::string& output_file_path)
{
    // Instantiate an `spdlog` thread pool for background output operations. The logger constructed with this
    // pool may contain up-to `QUEUE_CAP` log message units. If each log unit has a maximum length of `MSG_LEN`,
    // then the logger can take up memory up-to `(QUEUE_CAP x MSG_LEN)` in background.
    tp = std::make_shared<spdlog::details::thread_pool>(QUEUE_CAP, NUM_THREADS);

    // Instantiate an asynchronous `spdlog` logger that uses the thread pool `tp`, and writes to the provided
    // sink file at `output_file_path`.
    std::shared_ptr<spdlog::sinks::basic_file_sink_mt> sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(output_file_path);
    logger = std::make_shared<spdlog::async_logger>("async_output", sink, tp, spdlog::async_overflow_policy::block);

    // Set the log message pattern.
    logger->set_pattern("%v");
}


void Async_Logger_Wrapper::close_logger()
{
    // Note: For `spdlog`, `logger->flush()` posts a message to the queue requesting the flush operation,
    // so the function returns immediately. Hence a forceful eviction is necessary by dropping the `spdlog`
    // thread pools for the output to force-flush the pending logs. `spdlog::shutdown()` force-flushes
    // messages from the global pool only.

    logger->flush();
    tp.reset();
}
