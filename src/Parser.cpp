
#include "Parser.hpp"
#include "kseq/kseq.h"

#include <iostream>


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(gzFile, gzread);


Parser::Parser(const char* file_path):
    file_ptr(gzopen(file_path, "r")),   // Open the file handler.
    parser(kseq_init(file_ptr))   // Initialize the kseq parser.
{
    if(file_ptr == nullptr)
    {
        std::cerr << "Error opening input file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


Parser::Parser(const std::string& file_path):
    Parser(file_path.c_str())
{}


bool Parser::read_next_seq() const
{
    return (kseq_read(parser) >= 0);
}


const char* Parser::seq() const
{
    return parser->seq.s;
}


size_t Parser::seq_len() const
{
    return parser->seq.l;
}


size_t Parser::buff_sz() const
{
    return parser->seq.m;
}


void Parser::close()
{
    kseq_destroy(parser);   // Close the kseq parser.
    gzclose(file_ptr);  // Close the file handler.
}
