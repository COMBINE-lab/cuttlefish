
#ifndef PARSER_HPP
#define PARSER_HPP



#include <string>
#include <queue>
#include "zlib.h"


struct _KSEQ_DATA;


// Wrapper class to parse FASTA/FASTQ files using the kseq library.
class Parser
{
    typedef _KSEQ_DATA kseq_t;

private:

    std::queue<std::string> ref_file_paths;    // Collection of the reference file paths.
    gzFile file_ptr;  // Pointer to the reference file being parsed.
    kseq_t* parser;   // The kseq parser for the reference file being parsed.


    // Opens the reference at path `reference_path`.
    void open_reference(const std::string& reference_path);

    // Opens the next reference to be parsed from the collection `ref_file_paths`.
    bool open_next_reference();


public:

    // Constructs a parser for the file at path `file_path`. Iff `is_list`
    // is `true`, then the file is treated as a collection of reference paths.
    Parser(const std::string& file_path, const bool is_list);

    // If sequences are remaining to be read, reads the next one
    // into memory and returns `true`. Returns `false` otherwise.
    bool read_next_seq();

    // Returns a pointer to the current sequence in the buffer.
    const char* seq() const;

    // Returns the length of the current sequence in the buffer.
    size_t seq_len() const;

    // Returns the current size of the buffer.
    size_t buff_sz() const;

    // Closes the internal kseq parser for the current reference.
    void close();
};



#endif
