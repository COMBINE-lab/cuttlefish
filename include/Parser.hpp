
#ifndef PARSER_HPP
#define PARSER_HPP



#include <string>
#include "zlib.h"


struct _KSEQ_DATA;


// Wrapper class to parse FASTA/FASTQ files using the kseq library.
class Parser
{
    typedef _KSEQ_DATA kseq_t;

private:

    gzFile const file_ptr;  // Pointer to the file to be parsed.
    kseq_t* const parser;   // The kseq parser.


public:

    // Constructs a parser to parse the file located `file_path`.
    Parser(const char* file_path);

    // Constructs a parser to parse the file located `file_path`.
    Parser(const std::string& file_path);

    // If sequences are remaining to be read, reads the next one
    // into memory and returns `true`. Returns `false` otherwise.
    bool read_next_seq() const;

    // Returns a pointer to the current sequence in the buffer.
    const char* seq() const;

    // Returns the length of the current sequence in the buffer.
    size_t seq_len() const;

    // Returns the current size of the buffer.
    size_t buff_sz() const;

    // Closes the parser.
    void close();
};



#endif
