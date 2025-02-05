/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#ifndef H_FORMATER
#define H_FORMATER

#include <vector>
#include <string>

#include "Globals.h"
#include "Common.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "FastxStream.h"
#include "FastxChunk.h"
//#include "Sequence.h"
#include "Reference.h"

namespace rabbit {

namespace fa {

std::string getSequence(FastaDataChunk *&chunk, uint64 &pos);  // addbyxxm
std::string getLine(FastaDataChunk *&chunk, uint64 &pos);
int chunkListFormat(FastaChunk &fachunk, std::vector<Reference> &refs);
int chunkFormat(FastaChunk &fachunk, std::vector<Reference> &refs);
int chunkFormat(FastaChunk &fachunk, std::vector<Reference> &refs, int kmerSize);
Reference getNextSeq(FastaChunk &fachunk, bool &done, uint64 &pos);

}  // namespace fa

namespace fq {

int chunkFormat(FastqChunk *chunk, std::vector<Reference> &, bool);
int chunkFormat(FastqChunk *chunk, std::vector<neoReference> &);
int chunkFormat(FastqDataChunk *fqChunk, std::vector<neoReference> &data);
int chunkFormat(FastqDataChunk *fqChunk, std::vector<Reference> &data, bool);

std::string getLine(FastqDataChunk *&chunk, int &pos);
int neoGetLine(FastqDataChunk *&chunk, uint64_t &pos, uint64_t &len);

}  // namespace fq

}  // namespace rabbit

#endif
