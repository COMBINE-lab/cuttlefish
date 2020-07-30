
#include "CdBG.hpp"
#include "cxxopts/cxxopts.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>


// Driver function for the CdBG build.
void build(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish", "Efficiently construct the compacted de Brijn graph from references");
    options.add_options()
        ("r,ref", "reference sequence (in FASTA)", cxxopts::value<std::string>())
        ("k,kmer_len", "k-mer length", cxxopts::value<uint16_t>())
        ("d,kmc_db", "KMC database prefix", cxxopts::value<std::string>())
        ("t,threads", "number of threads to use", cxxopts::value<uint16_t>()->default_value("1"))
        ("o,output", "output file", cxxopts::value<std::string>())
        ("f,format", "output format (0: txt, 1: GFAv1, 2: GFAv2)", cxxopts::value<uint16_t>()->default_value("0"))
        ("w,work_dir", "working directory", cxxopts::value<std::string>()->default_value("."))
        ("b,bbhash", "BBHash file (optional)", cxxopts::value<std::string>()->default_value(""))
        ("h,help", "print usage");

    try
    {
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return;
        }

        auto ref = result["ref"].as<std::string>();
        auto k = result["kmer_len"].as<uint16_t>();
        auto kmer_database = result["kmc_db"].as<std::string>();
        auto thread_count = result["threads"].as<uint16_t>();
        auto output_file = result["output"].as<std::string>();
        auto format = result["format"].as<uint16_t>();
        auto working_dir = result["work_dir"].as<std::string>();
        auto bbhash_file = result["bbhash"].as<std::string>();

        if((k & 1) == 0)  // Discard even `k`.
        {
            std::cout << "The k-mer length (k) needs to be odd." << std::endl;
            return;
        }

        if(format > 2)  // Discard invalid output formats.
        {
            std::cout << "Invalid output file format." << std::endl;
            return;
        }

        
        std::cout << "Constructing compacted de Bruijn graph for the reference at " << ref << ", with k = " << k << "\n";

        CdBG cdbg(ref, k, kmer_database);
        cdbg.construct(bbhash_file, thread_count, output_file, format, working_dir);

        std::cout << "Constructed the compacted de Bruijn graph at " << output_file << "\n";
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << std::endl << "Usage :" << std::endl;
        std::cerr << options.help() << std::endl;
    }
}


// Driver function for the CdBG validation.
int validate(int argc, char** argv);


int main(int argc, char** argv)
{
    build(argc, argv);

    return EXIT_SUCCESS;
}