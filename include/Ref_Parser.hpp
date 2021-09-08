
#ifndef REF_PARSER_HPP
#define REF_PARSER_HPP



#include "Sequence_Input.hpp"

#include <string>
#include <queue>
#include "zlib.h"


struct _KSEQ_DATA;  // Forward declaration for `kseq`'s sequence-data format.


// Wrapper class to parse FASTA/FASTQ files using the `kseq` library.
class Ref_Parser
{
    typedef _KSEQ_DATA kseq_t;

private:

    std::queue<std::string> ref_paths;  // Collection of the reference file paths.
    gzFile file_ptr = nullptr;  // Pointer to the reference file being parsed.
    kseq_t* parser = nullptr;   // The kseq parser for the reference file being parsed.

    std::string curr_ref_path;  // Path to the reference currently being parsed.
    uint64_t ref_count = 0; // Number of the reference currently being parsed.
    uint64_t seq_id_; // Number of the current sequence (in the current reference).


    // Opens the reference at path `reference_path`.
    void open_reference(const std::string& reference_path);

    // Opens the next reference to be parsed from the collection `ref_file_paths`.
    bool open_next_reference();


public:

    // Constructs a parser for the file at path `file_path`.
    Ref_Parser(const std::string& file_path);

    // Constructs a parser for the reference input collection present at `ref_input`.
    Ref_Parser(const Sequence_Input& ref_input);

    // Returns the path to the reference currently being parsed.
    const std::string& curr_ref() const;

    // If sequences are remaining to be read, reads the next one
    // into memory and returns `true`. Returns `false` otherwise.
    bool read_next_seq();

    // Returns a pointer to the current sequence in the buffer.
    const char* seq() const;

    // Returns the length of the current sequence in the buffer.
    size_t seq_len() const;

    // Returns the current size of the buffer.
    size_t buff_sz() const;

    // Returns the id (number) of the current reference being parsed.
    uint64_t ref_id() const;

    // Returns the id (number) of the current sequence in the buffer.
    uint64_t seq_id() const;

    // Returns the name (as parsed) of the current sequence in the buffer.
    const char* seq_name() const;

    // Closes the internal kseq parser for the current reference.
    void close();
};



#endif
