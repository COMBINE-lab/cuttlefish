
#include <cstdlib>
#include <ctime>
#include <experimental/filesystem>


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


    return std::string(s);
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
    for(const auto& entry: std::experimental::filesystem::directory_iterator(path))
        if(is_prefix(entry.path(), prefix))
            return true;

    return false;
}
