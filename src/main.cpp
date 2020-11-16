
#include "Input_Defaults.hpp"
#include "CdBG.hpp"
#include "Validator.hpp"
#include "Build_Params.hpp"
#include "Validation_Params.hpp"
#include "Application.hpp"
#include "cxxopts/cxxopts.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>


// Driver function for the CdBG build.
void build(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish build", "Efficiently construct the compacted de Bruijn graph from references");
    options.add_options()
        ("r,refs", "reference files", cxxopts::value<std::vector<std::string>>()->default_value(cuttlefish::_default::EMPTY))
        ("l,lists", "reference file lists", cxxopts::value<std::vector<std::string>>()->default_value(cuttlefish::_default::EMPTY))
        ("d,dirs", "reference file directories", cxxopts::value<std::vector<std::string>>()->default_value(cuttlefish::_default::EMPTY))
        ("k,kmer_len", "k-mer length", cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::K)))
        ("s,kmc_db", "set of k-mers (KMC database) prefix", cxxopts::value<std::string>())
        ("t,threads", "number of threads to use", cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::THREAD_COUNT)))
        ("o,output", "output file", cxxopts::value<std::string>())
        ("f,format", "output format (0: txt, 1: GFA 1.0, 2: GFA 2.0)", cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::OP_FORMAT)))
        ("w,work_dir", "working directory", cxxopts::value<std::string>()->default_value(cuttlefish::_default::WORK_DIR))
        ("mph", "minimal perfect hash (BBHash) file (optional)", cxxopts::value<std::string>()->default_value(cuttlefish::_default::EMPTY))
        ("buckets", "hash table buckets (cuttlefish) file (optional)", cxxopts::value<std::string>()->default_value(cuttlefish::_default::EMPTY))
        ("h,help", "print usage");

    try
    {
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return;
        }

        const auto refs = result["refs"].as<std::vector<std::string>>();
        const auto lists = result["lists"].as<std::vector<std::string>>();
        const auto dirs = result["dirs"].as<std::vector<std::string>>();
        const auto k = result["kmer_len"].as<uint16_t>();
        const auto kmer_database = result["kmc_db"].as<std::string>();
        const auto thread_count = result["threads"].as<uint16_t>();
        const auto output_file = result["output"].as<std::string>();
        const auto format = result["format"].as<uint16_t>();
        const auto working_dir = result["work_dir"].as<std::string>();
        const auto mph_file = result["mph"].as<std::string>();
        const auto buckets_file = result["buckets"].as<std::string>();

        const Build_Params params(refs, lists, dirs, k, kmer_database, thread_count, output_file, format, working_dir, mph_file, buckets_file);
        if(!params.is_valid())
        {
            std::cerr << "Invalid input configuration. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
        

        std::cout << "\nConstructing the compacted de Bruijn graph for k = " << k << ".\n";

        const Application<cuttlefish::MAX_K> app(params);
        app.execute();

        std::cout << "\nConstructed the compacted de Bruijn graph at " << output_file << "\n";
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << std::endl << "Usage :" << std::endl;
        std::cerr << options.help() << std::endl;
    }
}


// Driver function for the CdBG validation.
void validate(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish validate", "Validate a compacted de Bruijn graph constructed by cuttlefish");
    options.add_options()
        ("r,refs", "reference files", cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("l,lists", "reference file lists", cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("d,dirs", "reference file directories", cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("k,kmer_len", "k-mer length", cxxopts::value<uint16_t>())
        ("s,kmc_db", "set of k-mers (KMC database) prefix", cxxopts::value<std::string>())
        ("g,cdbg", "compacted de Bruijn graph file", cxxopts::value<std::string>())
        ("t,threads", "number of threads to use", cxxopts::value<uint16_t>()->default_value("1"))
        ("w,work_dir", "working directory", cxxopts::value<std::string>()->default_value("."))
        ("mph", "minimal perfect hash (BBHash) file (optional)", cxxopts::value<std::string>()->default_value(""))
        ("h,help", "print usage");

    try
    {
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return;
        }

        const auto refs = result["refs"].as<std::vector<std::string>>();
        const auto lists = result["lists"].as<std::vector<std::string>>();
        const auto dirs = result["dirs"].as<std::vector<std::string>>();
        const auto k = result["kmer_len"].as<uint16_t>();
        const auto kmer_database = result["kmc_db"].as<std::string>();
        const auto cdbg = result["cdbg"].as<std::string>();
        const auto thread_count = result["threads"].as<uint16_t>();
        const auto working_dir = result["work_dir"].as<std::string>();
        const auto mph_file = result["mph"].as<std::string>();


        const Validation_Params params(refs, lists, dirs, k, kmer_database, cdbg, thread_count, working_dir, mph_file);
        if(!params.is_valid())
        {
            std::cerr << "Invalid input configuration. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        std::cout << "\nValidating the compacted de Bruijn graph for k = " << k << "\n";

        const Application<cuttlefish::MAX_K> app(params);
        std::cout << (app.validate() ? "\nValidation successful" : "\nValidation failed") << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << std::endl << "Usage :" << std::endl;
        std::cerr << options.help() << std::endl;
    }
}


int main(int argc, char** argv)
{
    if(argc < 2)
        std::cout << "Usage:\ncuttlefish <command> [OPTIONS]" << std::endl;
    else
    {
        const std::string command = argv[1];
        if(command == "build")
            build(argc - 1, argv + 1);
        else if(command == "validate")
            validate(argc - 1, argv + 1);
        else
            std::cout << "Invalid command. Supported commands: `build` and `validate`" << std::endl;
    }

    return EXIT_SUCCESS;
}