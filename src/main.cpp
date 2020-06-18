
#include "CdBG.hpp"
#include "CdBG_Builder.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>

int main(int argc, const char** argv)
{
    if(argc != 6)
    {
        std::cerr << "Command format: cuttlefish <ref_file> <k> <kmer_database> <thread_count> <out_file>\n";
        std::exit(EXIT_FAILURE);
    }


    std::string refs(argv[1]);
    uint16_t k(std::atoi(argv[2]));
    std::string kmer_database(argv[3]);
    uint16_t thread_count(std::atoi(argv[4]));
    std::string output_file(argv[5]);


    std::cout << "Constructing compacted de Bruijn graph for references at " << argv[1] << ", with k = " << k << "\n";

    CdBG_Builder cdbg(refs, k);

    cdbg.construct(kmer_database, thread_count, output_file);

    std::cout << "Constructed the compacted de Bruijn graph.\n";


    return EXIT_SUCCESS;
}