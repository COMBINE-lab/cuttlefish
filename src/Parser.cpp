
#include "Parser.hpp"
#include "kseq/kseq.h"

#include <fstream>
#include <iostream>


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(gzFile, gzread);


void Parser::open_reference(const std::string& reference_path)
{
    file_ptr = gzopen(reference_path.c_str(), "r");  // Open the file handler.
    if(file_ptr == nullptr)
    {
        std::cerr << "Error opening input file " << reference_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    parser = kseq_init(file_ptr);   // Initialize the kseq parser.
}


bool Parser::open_next_reference()
{
    if(ref_file_paths.empty())
        return false;


    open_reference(ref_file_paths.front());
    ref_file_paths.pop();

    return true;
}


Parser::Parser(const std::string& file_path, const bool is_list)
{
    if(!is_list)
        ref_file_paths.push(file_path);
    else
    {
        std::ifstream input(file_path.c_str(), std::ifstream::in);
        if(input.fail())
        {
            std::cerr << "Error opening input file " << file_path << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }


        std::string ref_path;
        while(input >> ref_path)
            ref_file_paths.push(ref_path);

        input.close();
    }


    open_next_reference();
}


bool Parser::read_next_seq()
{
    // Sequences still remain at the current reference being parsed.
    if(kseq_read(parser) >= 0)
        return true;

    // The current reference has been parsed completely. Close its handles.
    close();

    // Start parsing the next reference, if exists.
    if(open_next_reference())
        return kseq_read(parser) >= 0;

    return false;
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
    if(file_ptr != nullptr)
    {
        kseq_destroy(parser);   // Close the kseq parser.
        gzclose(file_ptr);  // Close the file handler.

        file_ptr = nullptr;
    }
}
