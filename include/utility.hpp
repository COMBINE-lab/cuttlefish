
#ifndef UTILITY_HPP
#define UTILITY_HPP



#include <cstddef>
#include <string>
#include <vector>

// TODO: wrap everything here in some namespaces.

// Returns a random string of length `len`, using characters from `alphabet`.
std::string get_random_string(size_t len, const char* alphabet =    "0123456789"
                                                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                    "abcdefghijklmnopqrstuvwxyz");

// Returns `true` iff `pref` is a prefix of `s`.
bool is_prefix(const std::string& s, const std::string& pref);

// Returns `true` iff there exists a file in the file system with the path
// `file_path`.
bool file_exists(const std::string& file_path);

// Returns `true` iff these exists a directory in the file system with the
// path `dir_path`.
bool dir_exists(const std::string& dir_path);

// Returns the file size is bytes of the file at path `file_path`. Returns
// `0` in case the file does not exist.
std::size_t file_size(const std::string& file_path);

// Returns `true` iff there exists some file in the file system path
// `path` with its name being prefixed by `prefix`.
bool file_prefix_exists(const std::string& path, const std::string& prefix);

// Returns a string that is a copy of `s` but has all the whitespaces removed.
std::string remove_whitespaces(const char* s);

// Given the collection of strings `s`, returns the concatenated string
// `s0 : s1 : ... : s_m`, where successive strings are separated by `delimiter`.
const std::string concat_strings(const std::vector<std::string>& s, const std::string& delimiter = ", ");

// Removes the file at path `file_path` from disk. Returns `true` iff the
// removal is successful.
bool remove_file(const std::string& file_path);

// Clears the content of the file at path `file_path`.
void clear_file(const std::string& file_path);

// Returns the name of the file present at the path `file_path`.
const std::string filename(const std::string& file_path);

// Returns the directory of the file present at the path `file_path`.
const std::string dirname(const std::string& file_path);

// Moves the file present at path `from_path` to the path `to_path`.
void move_file(const std::string& from_path, const std::string& to_path);

// Returns the maximum memory ("high-water-mark") used by the running
// process in bytes. Returns `0` in case of errors encountered.
std::size_t process_peak_memory();



#endif
