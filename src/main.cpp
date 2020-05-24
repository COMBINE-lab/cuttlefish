
#include "CdBG.hpp"
#include "CdBG_Builder.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>

int main(int argc, const char** argv)
{
    std::string refs(argv[1]);
    uint16_t k(std::atoi(argv[2]));
    std::string output_file(argv[3]);

    if(argc != 4)
    {
        std::cerr << "Command format: cuttlefish <ref_file> <k> <out_file>\n";
        std::exit(EXIT_FAILURE);
    }


    std::cout << "Constructing compacted de Bruijn graph for references at " << argv[1] << ", with k = " << k << "\n";

    CdBG_Builder cdbg(refs, k);

    // Classify the vertices.
    cdbg.construct(output_file);

    std::cout << "Constucted the compacted de Bruijn graph.\n";

    // For debugging.
    // cdbg.print_vertices();


    return EXIT_SUCCESS;
}