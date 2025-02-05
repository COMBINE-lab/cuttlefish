#ifndef __REFERENCE_H_
#define __REFERENCE_H_

#include <string>
#include <vector>
#include "Globals.h"
/**
 * @brief Reference struct that store the FASTA and FASTQ infomation
 */
struct Reference {
  std::string name;
  std::string comment;
  std::string seq;
  std::string quality;
  std::string strand;
  uint64_t length;
  uint64_t gid;
};

/**
 * @brief Reference struct that store the FASTA and FASTQ infomation, different from Reference
 *        neoReference only record the start position and length of name, sequence, strand and quality
 *        base on certain chunk data `base`
 */
struct neoReference {
  uint64_t pname;  /// name offset form base
  uint64_t pcom;   /// comment offset form base
  uint64_t pseq;   /// sequence offset form base
  uint64_t pqual;  /// quality offset form base
  uint64_t pstrand; ///strand offset form base

  uint64_t lname;   /// length of name
  uint64_t lcom;    /// length of comment
  uint64_t lseq;    /// length of sequence
	uint64_t lqual;   /// length of quality
  uint64_t lstrand; /// length of strand
  uint64_t gid;     /// global id
  rabbit::byte *base; /// base data pointer
};

/// Sequence infomations, only for FASTA data
typedef std::vector<Reference> SeqInfos;
/// One sequence sequence infomation, only for FASTA data
typedef Reference OneSeqInfo;

//TODO: it should support negtive number, but may cause some problem.
inline void cutseq_neoref(neoReference& ref, uint64_t from, uint64_t to){
	ref.pseq += from;
	ref.pqual += from;
	ref.lseq += (to - from + 1);
}

#endif
