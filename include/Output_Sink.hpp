
#ifndef OUTPUT_SINK_HPP
#define OUTPUT_SINK_HPP



#include "Async_Logger_Wrapper.hpp"
#include "spdlog/spdlog.h"

#include <fstream>
#include <string>


// A basic sink wrapper with minimal functionality — open, get reference to the wrapped sink, and close.
template <typename T_sink_>
class Output_Sink
{};


template <>
class Output_Sink<std::ofstream>
{
private:

    std::ofstream output_;


public:

    void init_sink(const std::string& output_file_path)
    {
        output_ = std::ofstream(output_file_path);
    }

    std::ofstream& sink()
    {
        return output_;
    }

    void close_sink()
    {
        output_.close();
    }
};


// TODO: redesign solution with a more informed scheme and skipping redundant `std::string` copying
// under the hood.
// Ref: https://github.com/gabime/spdlog/issues/643
template <>
class Output_Sink<Async_Logger_Wrapper>
{
private:

    Async_Logger_Wrapper output_;


public:

    void init_sink(const std::string& output_file_path)
    {
        output_.init_logger(output_file_path);
    }

    Async_Logger_Wrapper& sink()
    {
        return output_;
    }

    void close_sink()
    {
        output_.close_logger();
    }
};



#endif
