
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


void remove_kmer_set(const std::string& kmc_file_pref)
{
    const std::string kmc_file1_path(kmc_file_pref + ".kmc_pre");
    const std::string kmc_file2_path(kmc_file_pref + ".kmc_suf");

    if(std::remove(kmc_file1_path.c_str()) || std::remove(kmc_file2_path.c_str()))
    {
        std::cerr << "Error removing the KMC database file from path prefix " << kmc_file_pref << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}
