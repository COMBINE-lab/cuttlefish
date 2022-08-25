
#include "utility.hpp"

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
    return std::filesystem::remove(file_path);
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


std::size_t process_peak_memory()
{
    constexpr const char* process_file = "/proc/self/status";
    constexpr const char* peak_mem_field = "VmHWM:";
    const std::size_t field_len = std::strlen(peak_mem_field);

    std::FILE* fp = std::fopen(process_file, "r");
    if(fp == NULL)
    {
        std::cerr << "Error opening the process information file.\n";
        return 0;
    }

    char line[1024];
    std::size_t peak_mem = 0;
    while(std::fgets(line, sizeof(line) - 1, fp))
        if(std::strncmp(line, peak_mem_field, field_len) == 0)
        {
            peak_mem = std::strtoul(line + field_len, NULL, 0);
            break;
        }

    
    if(std::ferror(fp))
    {
        std::cerr << "Error reading the process information file.\n";
        return 0;
    }


    return peak_mem * 1024;
}
