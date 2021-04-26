
#ifndef UTILITY_HPP
#define UTILITY_HPP


#include <string>


// Returns a random string of length `len`, using characters from `alphabet`.
std::string get_random_string(size_t len, const char* alphabet =    "0123456789"
                                                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                    "abcdefghijklmnopqrstuvwxyz");

// Returns `true` iff `pref` is a prefix of `s`.
bool is_prefix(const std::string& s, const std::string& pref);

// Returns `true` iff there exists a file in the file system with the path
// `file_path`.
bool file_exists(const std::string& file_path);

// Returns `true` iff there exists some file in the file system path
// `path` with its name being prefixed by `prefix`.
bool file_prefix_exists(const std::string& path, const std::string& prefix);

// Returns a string that is a copy of `s` but has all the whitespaces removed.
std::string remove_whitespaces(const char* s);

// Removes the k-mer set (KMC database) with the path prefix `kmc_file_pref`.
void remove_kmer_set(const std::string& kmc_file_pref);



#endif
