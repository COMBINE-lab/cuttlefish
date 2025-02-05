/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#ifndef H_GLOBALS
#define H_GLOBALS

#ifndef NDEBUG
#define DEBUG 1
#endif

// Visual Studio warning supression
//
#if defined(_WIN32)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)  // D_SCL_SECURE
#pragma warning(disable : 4244)  // conversion uint64 to uint32
#pragma warning(disable : 4267)
#pragma warning(disable : 4800)  // conversion byte to bool
#endif

// assertions
//
#if defined(DEBUG) || defined(_DEBUG)
#include "assert.h"
#define ASSERT(x) assert(x)
#else
#define ASSERT(x) (void)(x)
#endif

#include <string>
#include <stdexcept>

namespace rabbit {

// basic types
//
typedef char int8;
typedef unsigned char uchar, byte, uint8;
typedef short int int16;
typedef unsigned short int uint16;
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

// exception class
//
class RioException : public std::exception {
  std::string message;

 public:
  RioException(const char *msg_) : message(msg_) {}

  RioException(const std::string &msg_) : message(msg_) {}

  ~RioException() throw() {}

  const char *what() const throw()  // for std::exception interface
  {
    return message.c_str();
  }
};

}  // namespace rabbit

#endif  // H_GLOBALS
