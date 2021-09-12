
#include "Ref_Parser.hpp"
#include "kseq/kseq.h"
#include "ghc/filesystem.hpp"

#include <fstream>
#include <iostream>


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(gzFile, gzread);


Ref_Parser::Ref_Parser(const std::string& file_path)
{
    ref_paths.push(file_path);

    // Open the first reference for subsequent parsing.
    open_next_reference();
}


Ref_Parser::Ref_Parser(const Seq_Input& ref_input): Ref_Parser(ref_input.seqs())
{
    // Open the first reference for subsequent parsing.
    open_next_reference();
}


Ref_Parser::Ref_Parser(const std::vector<std::string>& refs):
    ref_paths(std::deque<std::string>(refs.begin(), refs.end()))
{}


void Ref_Parser::open_reference(const std::string& reference_path)
{
    file_ptr = gzopen(reference_path.c_str(), "r");  // Open the file handler.
    if(file_ptr == nullptr)
    {
        std::cerr << "Error opening reference file " << reference_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    parser = kseq_init(file_ptr);   // Initialize the kseq parser.

    curr_ref_path = reference_path;
    ref_count++;
    seq_id_ = 0;

    std::cout << "\nOpened reference " << ref_count << " from " << curr_ref_path << "\n";
}


bool Ref_Parser::open_next_reference()
{
    if(ref_paths.empty())
        return false;


    open_reference(ref_paths.front());
    ref_paths.pop();

    return true;
}


const std::string& Ref_Parser::curr_ref() const
{
    return curr_ref_path;
}


bool Ref_Parser::read_next_seq()
{
    // Sequences still remain at the current reference being parsed.
    if(parser != nullptr && kseq_read(parser) >= 0)
    {
        seq_id_++;
        return true;
    }

    // The current reference has been parsed completely. Close its handles.
    close();

    // Start parsing the next reference, if exists.
    if(open_next_reference())
        return read_next_seq();

    return false;
}


const char* Ref_Parser::seq() const
{
    return parser->seq.s;
}


size_t Ref_Parser::seq_len() const
{
    return parser->seq.l;
}


size_t Ref_Parser::buff_sz() const
{
    return parser->seq.m;
}


uint64_t Ref_Parser::ref_id() const
{
    return ref_count;
}


uint64_t Ref_Parser::seq_id() const
{
    return seq_id_;
}


const char* Ref_Parser::seq_name() const
{
    return parser->name.s;
}


void Ref_Parser::close()
{
    if(file_ptr != nullptr)
    {
        kseq_destroy(parser);   // Close the kseq parser.
        gzclose(file_ptr);  // Close the file handler.

        parser = nullptr;
        file_ptr = nullptr;

        std::cerr << "\rClosed reference " << curr_ref() << ".\n";
    }
}
