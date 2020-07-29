
#ifndef UTILITY_HPP
#define UTILITY_HPP


#include <string>


// Returns a random string of length `len`.
std::string get_random_string(const size_t len);

// Returns `true` iff `pref` is a prefix of `s`.
bool is_prefix(const std::string& s, const std::string& pref);

// Returns `true` iff there exists some file in the file system path
// `path` with its name being prefixed by `prefix`.
bool file_prefix_exists(const std::string& path, const std::string& prefix);



#endif
