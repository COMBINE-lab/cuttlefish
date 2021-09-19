
#include "utility.hpp"
#include "ghc/filesystem.hpp"

#include <cctype>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iostream>


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
    struct stat stat_buf;

    return stat(file_path.c_str(), &stat_buf) == 0;
}


bool file_prefix_exists(const std::string& path, const std::string& prefix)
{
    for(const auto& entry: ghc::filesystem::directory_iterator(path))
        if(is_prefix(entry.path(), prefix))
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
    return ghc::filesystem::remove(file_path);
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
    return ghc::filesystem::path(file_path).filename().string();
}


void move_file(const std::string& from_path, const std::string& to_path)
{
    ghc::filesystem::copy(from_path, to_path);
    ghc::filesystem::remove(from_path);
}
