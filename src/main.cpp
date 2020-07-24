
#include "CdBG.hpp"
#include "cxxopts/cxxopts.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish", "efficiently constructed the compacted de Bruijn graph.");
    options.add_options()
    ("r,ref", "the refrence sequence (FASTA)", cxxopts::value<std::string>())
    ("k,klen", "the k-mer length", cxxopts::value<uint16_t>())
    ("d,database", "the KMC database", cxxopts::value<std::string>())
    ("t,threads", "the number of threads to use", cxxopts::value<uint16_t>())
    ("o,output", "the output file", cxxopts::value<std::string>())
    ("b,bbhash", "the BBHash file", cxxopts::value<std::string>())
    ("h,help", "print usage");

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        /*if(argc != 6)
    {
        std::cerr << "Command format: cuttlefish <ref_file> <k> <kmer_database> <thread_count> <out_file>\n";
        std::exit(EXIT_FAILURE);
    }*/

        auto refs = result["ref"].as<std::string>();
        auto k = result["klen"].as<uint16_t>();
        auto kmer_database = result["database"].as<std::string>();
        auto thread_count = result["threads"].as<uint16_t>();
        auto output_file = result["output"].as<std::string>();
        auto bbhash_file = result["bbhash"].as<std::string>();
        /*
    std::string refs(argv[1]);
    uint16_t k(std::atoi(argv[2]));
    std::string kmer_database(argv[3]);
    uint16_t thread_count(std::atoi(argv[4]));
    std::string output_file(argv[5]);
    */

        std::cout << "Constructing compacted de Bruijn graph for references at " << refs << ", with k = " << k << "\n";

        CdBG cdbg(refs, k, kmer_database);

        cdbg.construct(bbhash_file, thread_count, output_file);

        std::cout << "Constructed the compacted de Bruijn graph.\n";
    } catch (std::exception &e) {
        std::cerr << "[exception] : " << e.what() << "\n";
        std::cerr << "\nUsage : \n\n";
        std::cerr << options.help() << std::endl;
    }

    return EXIT_SUCCESS;
}