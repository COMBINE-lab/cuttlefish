
#include "utility.hpp"

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iterator>
#include <filesystem>
#include <fstream>
#include <cstdio>


std::string get_random_string(const size_t len, const char* const alphabet)
{
    std::string str;
    str.reserve(len);

    const unsigned int seed = static_cast<unsigned int>(std::time(NULL));
    std::srand(seed);
    for (size_t i = 0; i < len; ++i)
        str += alphabet[(std::rand() % (sizeof(alphabet) - 1))];

    return str;
}


bool is_prefix(const std::string& s, const std::string& pref)
{
    if(s.length() < pref.length())
        return false;

    size_t idx = 0;
    for(; idx < pref.length() && s[idx] == pref[idx]; ++idx);

    return idx == pref.length();
}


bool file_exists(const std::string& file_path)
{
    return std::filesystem::exists(file_path);
}


bool dir_exists(const std::string& dir_path)
{
    return std::filesystem::is_directory(dir_path);
}


std::size_t file_size(const std::string& file_path)
{
    // TODO: update policy for 0-sized files.
    std::error_code ec;
    const uintmax_t size = std::filesystem::file_size(file_path, ec);
    return ec ? 0 : static_cast<std::size_t>(size);
}


bool file_prefix_exists(const std::string& path, const std::string& prefix)
{
    for(const auto& entry: std::filesystem::directory_iterator(path))
        if(is_prefix(filename(entry.path()), prefix))
            return true;

    return false;
}


std::size_t load_file(const char* const file_path, char* const buf)
{
    std::ifstream input(file_path, std::ios::in | std::ios::binary);
    const auto sz = file_size(file_path);
    input.read(buf, sz);
    input.close();

    if(!input)
    {
        std::cerr << "Error reading of from file at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return sz;
}


std::size_t load_file(const std::string& file_path, char* const buf)
{
    return load_file(file_path.data(), buf);
}


void load_file(const std::string& file_path, const std::size_t sz, char* const buf)
{
    std::ifstream is(file_path, std::ios::binary);
    is.read(buf, sz);
    is.close();

    if(!is)
    {
        std::cerr << "Error reading " << sz << " bytes from file at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


std::string remove_whitespaces(const char* s)
{
    std::string str;
    str.reserve(strlen(s));

    for(const char* p = s; *p; ++p)
        if(!std::isspace(*p))
            str += *p;

    return str;
}


const std::string concat_strings(const std::vector<std::string>& s, const std::string& delimiter)
{
    std::ostringstream concat_stream;
    std::copy(s.begin(), s.end(), std::ostream_iterator<std::string>(concat_stream, delimiter.c_str()));

    std::string concat_str(concat_stream.str());
    concat_str.erase(concat_str.size() - delimiter.size(), delimiter.size());
    return concat_str;
}


bool remove_file(const std::string& file_path)
{
    std::error_code ec;
    return std::filesystem::remove(file_path, ec);
}


void clear_file(const std::string& file_path)
{
    std::ofstream file(file_path.c_str(), std::ofstream::out | std::ofstream::trunc);
    if(file.fail())
    {
        std::cerr << "Error opening file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    file.close();
}


const std::string filename(const std::string& file_path)
{
    return std::filesystem::path(file_path).filename().string();
}


const std::string dirname(const std::string& file_path)
{
    const std::string path = std::filesystem::path(file_path).remove_filename().string();
    return path.empty() ? "." : path;
}


void move_file(const std::string& from_path, const std::string& to_path)
{
    std::filesystem::copy(from_path, to_path);
    std::filesystem::remove(from_path);
}


std::size_t process_metric(const std::string& metric)
{
    const std::string process_file("/proc/self/status");
    const std::size_t field_len = metric.length();

    std::ifstream is(process_file);
    if(!is)
    {
        std::cerr << "Error opening the process information file.\n";
        return 0;
    }

    char line[1024];
    std::size_t val = 0;
    while(is.getline(line, sizeof(line) - 1), '\n')
        if(std::strncmp(line, metric.data(), field_len) == 0)
        {
            val = std::strtoul(line + field_len, NULL, 0);
            break;
        }

    is.close();
    if(!is)
    {
        std::cerr << "Error reading the process information file.\n";
        return 0;
    }


    return val;
}


std::size_t process_peak_memory()
{
    return process_metric("VmHWM:") * 1024;
}


std::size_t process_cur_memory()
{
    return process_metric("VmRSS:") * 1024;
}
