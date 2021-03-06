
#include "utility.hpp"
#include "ghc/filesystem.hpp"

#include <cctype>
#include <cstring>
#include <cstdlib>
#include <ctime>


std::string get_random_string(const size_t len)
{
    static const char alphabet[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    char *s = new char[len + 1];

    const unsigned seed = time(NULL);
    srand(seed);
    for (size_t i = 0; i < len; ++i)
        s[i] = alphabet[(std::rand() % (sizeof(alphabet) - 1))];

    s[len] = '\0';


    const std::string str(s);
    delete[] s;

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
