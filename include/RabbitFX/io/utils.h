/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#ifndef H_UTILS
#define H_UTILS

#include "Globals.h"

#include <string>

//-----------------
#include <stdlib.h>
#include <cstddef>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>
//-----------------

namespace rabbit {

namespace core {

template <uint32 _TBitNum>
class TBitMask {
 public:
  static const uint64 Value = ((uint64)1 << (_TBitNum - 1)) | TBitMask<_TBitNum - 1>::Value;
};

template <>
class TBitMask<0> {
 public:
  static const uint64 Value = 0;
};

template <uint32 _TNum>
struct TLog2 {
  static const uint32 Value = TLog2<(_TNum >> 1)>::Value + 1;
};

template <>
struct TLog2<1> {
  static const uint32 Value = 0;
};

template <>
struct TLog2<0> {};

template <typename _T>
inline void TSwap(_T &x_, _T &y_) {
  _T tmp = x_;
  x_ = y_;
  y_ = tmp;
}

template <typename _T>
inline void TFree(_T *p_) {
  if (p_ != (_T *)0) delete p_, p_ = (_T *)0;
}

inline uint32 to_string(uchar *str, uint32 value) {
  uint32 digits;
  uint32 power = 1;

  if (value == 0) {
    str[0] = '0';
    return 1;
  }

  for (digits = 0; digits < 10; ++digits) {
    if (value < power) break;
    power *= 10;
  }

  power /= 10;
  for (uint32 i = 0; power; ++i, power /= 10) {
    int32 d = value / power;
    str[i] = (uchar)('0' + d);
    value -= d * power;
  }

  return digits;
}

inline bool extend_string(uchar *&str, uint32 &size) {
  uint32 new_size = size * 2;
  uchar *p = new uchar[new_size + 1];

  if (!p) return false;

  std::copy(str, str + size, p);
  size = new_size;
  delete[] str;
  str = p;

  return true;
}

inline bool extend_string_to(uchar *&str, uint32 &size, uint32 new_size) {
  if (new_size <= size) return true;

  uchar *p = new uchar[new_size + 1];

  if (!p) return false;

  std::copy(str, str + size, p);

  size = new_size;
  delete[] str;
  str = p;

  return true;
}

inline bool ends_with(const std::string &str_, const std::string &suff_) {
  return str_.size() >= suff_.size() && str_.compare(str_.size() - suff_.size(), suff_.size(), suff_) == 0;
}

inline uint32 int_log(uint32 x, uint32 base) {
  uint32 r = 0;

  if (base == 0) return 1;
  if (base == 1) base++;

  for (uint32 tmp = base; tmp <= x; tmp *= base) ++r;

  return r;
}

inline uint32 to_num(const uchar *str, uint32 len) {
  uint32 r = 0;

  for (uint32 i = 0; i < len; ++i) r = r * 10 + (str[i] - '0');

  return r;
}

inline bool is_num(const uchar *str_, uint32 len_, uint32 &val_) {
  val_ = 0;
  uint32 i;
  for (i = 0; i < len_; ++i) {
    if (str_[i] < '0' || str_[i] > '9') break;
    val_ = val_ * 10 + (str_[i] - '0');
  }

  return i == len_ && (len_ == 1 || str_[0] != '0');
}

inline uint32 bit_length(uint64 x) {
  for (uint32 i = 0; i < 32; ++i) {
    if (x < (1ull << i)) return i;
  }
  return 64;
}

inline uint64_t seq2int(const char* data, int start, int keylen, bool& valid) {
	uint8_t mask = 0x06; //not general only works for DNA sequences, it's just a trick.
	uint64_t res = 0;
	const int end = start + keylen;
	for(int i = start; i < end; i++)
	{
		uint8_t meri = (uint8_t)data[i];
    if(data[i] == 'N'){
        valid = 0;
        return 0;
    }
		meri &= mask;
		meri >>= 1;
		res |= (uint64_t)meri;
		res <<= 2;
	}
	return res >> 2;
}

const int revmap[4] = {2,3,0,1};
//NOTICE: ensure there is no 'N' in key
inline uint64_t kmer_reverse_complete(uint64_t key, int keylen){
 	uint64_t res = 0;
 	for(int i = 0; i < keylen; ++i){
 		res |= revmap[key & 0x03];
 		res <<= 2;
 		key >>= 2;
 	}
 	return res;
}
//^0x02 -> complete
inline void reverse_complement(const char * src, char * dest, int length)
{
  char table[4] = {'T','G','A','C'};
  for ( int i = 0; i < length; i++ )
  {
    char base = src[i];

    base >>= 1;
    base &= 0x03;
    dest[length - i - 1] = table[static_cast<std::size_t>(base)];
  }
}

inline void seq_to_lower(char* seq, size_t len){
  for(std::size_t i = 0; i < len; ++i)
    seq[i] |= 0x20;
}
inline void seq_to_upper(char* seq, size_t len){
  for(std::size_t i = 0; i < len; ++i)
    seq[i] &= 0xdf;
}
inline void seq_to_lower(std::string &seq){
  for(char &c : seq)
    c |= 0x20;
}
inline void seq_to_upper(std::string &seq){
  for(char &c : seq)
    c &= 0xdf;
}

}  // namespace core
}  // namespace rabbit


//using namespace std;

namespace std{
inline char complement(char base) {
  switch (base) {
    case 'A':
    case 'a':
      return 'T';
    case 'T':
    case 't':
      return 'A';
    case 'C':
    case 'c':
      return 'G';
    case 'G':
    case 'g':
      return 'C';
    default:
      return 'N';
  }
}

inline bool starts_with(string const &value, string const &starting) {
  if (starting.size() > value.size()) return false;
  return equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with(string const &value, string const &ending) {
  if (ending.size() > value.size()) return false;
  return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string &str) {
  string::size_type pos = str.find_first_not_of(' ');
  if (pos == string::npos) {
    return string("");
  }
  string::size_type pos2 = str.find_last_not_of(' ');
  if (pos2 != string::npos) {
    return str.substr(pos, pos2 - pos + 1);
  }
  return str.substr(pos);
}

inline int split(const string &str, vector<string> &ret_, string sep = ",") {
  if (str.empty()) {
    return 0;
  }

  string tmp;
  string::size_type pos_begin = str.find_first_not_of(sep);
  string::size_type comma_pos = 0;

  while (pos_begin != string::npos) {
    comma_pos = str.find(sep, pos_begin);
    if (comma_pos != string::npos) {
      tmp = str.substr(pos_begin, comma_pos - pos_begin);
      pos_begin = comma_pos + sep.length();
    } else {
      tmp = str.substr(pos_begin);
      pos_begin = comma_pos;
    }

    ret_.push_back(tmp);
    tmp.clear();
  }
  return 0;
}

inline string replace(const string &str, const string &src, const string &dest) {
  string ret;

  string::size_type pos_begin = 0;
  string::size_type pos = str.find(src);
  while (pos != string::npos) {
    ret.append(str.data() + pos_begin, pos - pos_begin);
    ret += dest;
    pos_begin = pos + 1;
    pos = str.find(src, pos_begin);
  }
  if (pos_begin < str.length()) {
    ret.append(str.begin() + pos_begin, str.end());
  }
  return ret;
}

inline string reverse(const string &str) {
  string ret(str.length(), 0);
  for (std::size_t pos = 0; pos < str.length(); pos++) {
    ret[pos] = str[str.length() - pos - 1];
  }
  return ret;
}

inline string basename(const string &filename) {
  string::size_type pos = filename.find_last_of('/');
  if (pos == string::npos)
    return filename;
  else if (pos == filename.length() - 1)
    return "";  // a bad filename
  else
    return filename.substr(pos + 1, filename.length() - pos - 1);
}

inline string dirname(const string &filename) {
  string::size_type pos = filename.find_last_of('/');
  if (pos == string::npos) {
    return "./";
  } else
    return filename.substr(0, pos + 1);
}

inline string joinpath(const string &dirname, const string &basename) {
  if (dirname[dirname.length() - 1] == '/') {
    return dirname + basename;
  } else {
    return dirname + "/" + basename;
  }
}

// Check if a string is a file or directory
inline bool file_exists(const string &s) {
  bool exists = false;
  if (s.length() > 0) {
    struct stat status;
    int result = stat(s.c_str(), &status);
    if (result == 0) {
      exists = true;
    }
  }
  return exists;
}

// check if a string is a directory
inline bool is_directory(const string &path) {
  bool isdir = false;
  struct stat status;
  // visual studion use _S_IFDIR instead of S_IFDIR
  // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
  stat(path.c_str(), &status);
  if (status.st_mode & S_IFDIR) {
    isdir = true;
  }
  // #endif
  return isdir;
}

inline void check_file_valid(const string &s) {
  if (!file_exists(s)) {
    cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
    exit(-1);
  }
  if (is_directory(s)) {
    cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
    exit(-1);
  }
}

inline void check_file_writable(const string &s) {
  string dir = dirname(s);
  if (!file_exists(dir)) {
    cerr << "ERROR: '" << dir << " doesn't exist. Create this folder and run this command again." << endl;
    exit(-1);
  }
  if (is_directory(s)) {
    cerr << "ERROR: '" << s << "' is not a writable file, quit now" << endl;
    exit(-1);
  }
}

// Remove non alphabetic characters from a string
inline string str_keep_alpha(const string &s) {
  string new_str;
  for (size_t it = 0; it < s.size(); it++) {
    if (isalpha(s[it])) {
      new_str += s[it];
    }
  }
  return new_str;
}

// Remove invalid sequence characters from a string
inline string str_keep_valid_sequence(const string &s) {
  string new_str;
  for (size_t it = 0; it < s.size(); it++) {
    if (isalpha(s[it]) || s[it] == '-' || s[it] == '*') {
      new_str += s[it];
    }
  }
  return new_str;
}

inline int find_with_right_pos(const string &str, const string &pattern, int start = 0) {
  int pos = str.find(pattern, start);
  if (pos < 0)
    return -1;
  else
    return pos + pattern.length();
}

inline void str2upper(string &s) { transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper); }

inline void str2lower(string &s) { transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower); }

inline char num2qual(int num) {
  if (num > 127 - 33) num = 127 - 33;
  if (num < 0) num = 0;

  char c = num + 33;
  return c;
}

inline void error_exit(const string &msg) {
  cerr << "ERROR: " << msg << endl;
  exit(-1);
}

extern mutex logmtx;
inline void loginfo(const string s) {
  logmtx.lock();
  time_t tt = time(NULL);
  tm *t = localtime(&tt);
  cerr << "[" << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec << "] " << s << endl;
  logmtx.unlock();
}
} //namespace std

#endif
