/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

/*
  This file is modified for an efficient I/O fasta and fastq framework.
  Zekun Yin
  zekun.yin@mail.sdu.edu.cn
*/

#ifndef H_COMMON
#define H_COMMON

#include "Globals.h"

#ifndef NDEBUG
#define DEBUG 1
#endif

#define BIT(x) (1 << (x))
#define BIT_ISSET(x, pos) ((x & BIT(pos)) != 0)
#define BIT_SET(x, pos) (x |= BIT(pos))
#define BIT_UNSET(x, pos) (x &= ~(BIT(pos)))
#ifndef MIN
    #define MIN(x, y) ((x) <= (y) ? (x) : (y))
#endif
#define MAX(x, y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(x) ((x) >= 0 ? 1 : -1)
#define REC_EXTENSION_FACTOR(size) (((size) / 4 > 1024) ? ((size) / 4) : 1024)
#define MEM_EXTENSION_FACTOR(size) REC_EXTENSION_FACTOR(size)

#if defined(_WIN32)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)  // D_SCL_SECURE
#pragma warning(disable : 4244)  // conversion uint64 to uint32
#pragma warning(disable : 4267)
#pragma warning(disable : 4800)  // conversion byte to bool
#endif

// TODO: refactor raw data structs to avoid using <string> as a member
#include <string>

#define COMPILE_TIME_ASSERT(COND, MSG) typedef char static_assertion_##MSG[(!!(COND)) * 2 - 1]
#define COMPILE_TIME_ASSERT1(X, L) COMPILE_TIME_ASSERT(X, static_assertion_at_line_##L)
#define COMPILE_TIME_ASSERT2(X, L) COMPILE_TIME_ASSERT1(X, L)
#define STATIC_ASSERT(X) COMPILE_TIME_ASSERT2(X, __LINE__)

#endif  // _COMMON_H
