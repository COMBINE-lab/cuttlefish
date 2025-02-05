/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.1
  Date   : 2022-01-04
*/

#ifndef _KMER_H
#define _KMER_H

// Important remark: there is no inheritance here to guarantee that all classes defined here are POD according to C++11

#ifndef FORCE_INLINE
#ifdef _WIN32
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif
#endif // !FORCE_INLINE

#ifdef _WIN32
#define _bswap64(x) _byteswap_uint64(x)
#else
#define _bswap64(x) __builtin_bswap64(x)
#endif

#ifdef _WIN32
#define EXPECT(x, y)	(x)
#else
#define EXPECT(x, y)	(__builtin_expect(x, y))
#endif

#define USE_META_PROG

#include "meta_oper.h"
#include <string>
#include <random>
#include <ostream>
#include <emmintrin.h>

inline uint32_t get_rev_compl_shift(uint32_t len)
{
	return (len + 31) / 32 * 32 - len;
}

//instead of x << 2*p there is (x<<p) << p, because
//for p=32 we would have x << 64 which is UB
FORCE_INLINE uint64_t shl_2p(uint64_t x, uint64_t p) {
	return (x << p) << p;
}
//instead of x >> 2*p there is (x>>p) >> p, because
//for p=32 we would have x >> 64 which is UB
FORCE_INLINE uint64_t shr_2p(uint64_t x, uint64_t p) {
	return (x >> p) >> p;
}

// *************************************************************************
// CKmer class for k > 32 with classic kmer counting
template<unsigned SIZE> struct CKmer {
	unsigned long long data[SIZE];

	typedef unsigned long long data_t;
	static uint32_t KMER_SIZE;

	FORCE_INLINE void clear(void);

	FORCE_INLINE unsigned long long& operator[](int idx);
	FORCE_INLINE const unsigned long long& operator[](int idx) const;
};

// *********************************************************************
template<unsigned SIZE> FORCE_INLINE void CKmer<SIZE>::clear(void)
{
#ifdef USE_META_PROG
	IterFwd([&](const int &i){
		data[i] = 0;
		}, uint_<SIZE-1>());
#else
	for(uint32_t i = 0; i < SIZE; ++i)
		data[i] = 0;
#endif
}

template<unsigned SIZE> FORCE_INLINE unsigned long long& CKmer<SIZE>::operator[](int idx)
{
	return data[idx];
}

template<unsigned SIZE> FORCE_INLINE const unsigned long long& CKmer<SIZE>::operator[](int idx) const
{
	return data[idx];
}

#endif

// ***** EOF
