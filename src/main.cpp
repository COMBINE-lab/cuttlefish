
#include "Input_Defaults.hpp"
#include "CdBG.hpp"
#include "Read_CdBG.hpp"
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
#include <iomanip>


// Driver function for the CdBG build.
void build(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish build", "Efficiently construct the compacted de Bruijn graph from sequencing reads or reference sequences");

    options.add_options("common")
        ("r,refs", "input files", cxxopts::value<std::vector<std::string>>()->default_value(cuttlefish::_default::EMPTY))
        ("l,lists", "input file lists", cxxopts::value<std::vector<std::string>>()->default_value(cuttlefish::_default::EMPTY))
        ("d,dirs", "input file directories", cxxopts::value<std::vector<std::string>>()->default_value(cuttlefish::_default::EMPTY))
        ("k,kmer-len", "k-mer length", cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::K)))
        ("t,threads", "number of threads to use", cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::THREAD_COUNT)))
        ("o,output", "output file", cxxopts::value<std::string>())
        ("w,work-dir", "working directory", cxxopts::value<std::string>()->default_value(cuttlefish::_default::WORK_DIR))
        ("h,help", "print usage")
        ;

    options.add_options("cuttlefish 2.0")
        ("read", "construct a compacted read de Bruijn graph")
        ("ref", "construct a compacted reference de Bruijn graph")
        ("c,cutoff", "frequency cutoff for (k + 1)-mers", cxxopts::value<uint32_t>()->default_value(std::to_string(cuttlefish::_default::CUTOFF_FREQ)))
        ("m,max-memory", "soft maximum memory limit (in GB)", cxxopts::value<std::size_t>()->default_value(std::to_string(cuttlefish::_default::MAX_MEMORY)))
        ("unrestrict-memory", "do not impose memory usage restriction")
        ("path-cover", "extract a maximal path cover of the de Bruijn graph")
        ;
    
    options.add_options("cuttlefish 1.0")
        ("s,kmc-db", "set of vertices, i.e. k-mers (KMC database) prefix", cxxopts::value<std::string>()->default_value(cuttlefish::_default::WORK_DIR))
        ("f,format", "output format (0: txt, 1: GFA 1.0, 2: GFA 2.0, 3: GFA-reduced)", cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::OP_FORMAT)))
        ;

    options.add_options("specialized")
        // TODO: repurpose the following two options
        ("mph", "minimal perfect hash (BBHash) file (optional)", cxxopts::value<std::string>()->default_value(cuttlefish::_default::EMPTY))
        ("buckets", "hash table buckets (cuttlefish) file (optional)", cxxopts::value<std::string>()->default_value(cuttlefish::_default::EMPTY))
        ("save-vertices", "save the vertex set of the graph")
        ;

    options.add_options("debug")
        ("e,edge-db", "set of edges, i.e. (k + 1)-mers (KMC database) prefix", cxxopts::value<std::string>()->default_value(cuttlefish::_default::EMPTY))
#ifdef CF_DEVELOP_MODE
        ("gamma", "gamma for the BBHash MPHF", cxxopts::value<double>()->default_value(std::to_string(cuttlefish::_default::GAMMA)))
#endif
        ;

    try
    {
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return;
        }

        const auto is_read_graph = result["read"].as<bool>();
        const auto is_ref_graph = result["ref"].as<bool>();
        const auto refs = result["refs"].as<std::vector<std::string>>();
        const auto lists = result["lists"].as<std::vector<std::string>>();
        const auto dirs = result["dirs"].as<std::vector<std::string>>();
        const auto k = result["kmer-len"].as<uint16_t>();
        const auto cutoff = result["cutoff"].as<uint32_t>();
        const auto kmer_database = result["kmc-db"].as<std::string>();
        const auto edge_database = result["edge-db"].as<std::string>();
        const auto thread_count = result["threads"].as<uint16_t>();
        const auto max_memory = result["max-memory"].as<std::size_t>();
        const auto strict_memory = !result["unrestrict-memory"].as<bool>();
        const auto output_file = result["output"].as<std::string>();
        const auto format = result["format"].as<uint16_t>();
        const auto working_dir = result["work-dir"].as<std::string>();
        const auto path_cover = result["path-cover"].as<bool>();
        const auto mph_file = result["mph"].as<std::string>();
        const auto buckets_file = result["buckets"].as<std::string>();
        const auto save_vertices = result["save-vertices"].as<bool>();
#ifdef CF_DEVELOP_MODE
        const double gamma = result["gamma"].as<double>();
#endif

        const Build_Params params(  is_read_graph, is_ref_graph,
                                    refs, lists, dirs,
                                    k, cutoff, kmer_database, edge_database, thread_count, max_memory, strict_memory,
                                    output_file, format, working_dir,
                                    path_cover,
                                    mph_file, buckets_file, save_vertices
#ifdef CF_DEVELOP_MODE
                                    , gamma
#endif
                                );
        if(!params.is_valid())
        {
            std::cerr << "Invalid input configuration. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        // std::cout.precision(3);


        const std::string dBg_type(params.is_read_graph() ? "read" : "reference");

        std::cout << "\nConstructing the compacted " << dBg_type << " de Bruijn graph for k = " << k << ".\n";

        (params.is_read_graph() || params.is_ref_graph()) ?
            Application<cuttlefish::MAX_K, Read_CdBG>(params).execute() :
            Application<cuttlefish::MAX_K, CdBG>(params).execute();

        std::cout << "\nConstructed the " << dBg_type << " compacted de Bruijn graph at " << output_file << ".\n";
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

        std::cout << (Application<cuttlefish::MAX_K, CdBG>(params).validate() ?
                        "\nValidation successful" : "\nValidation failed") << std::endl;
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
#ifdef CF_DEVELOP_MODE
    std::cout << "Warning: Executing in Develop Mode.\n";
#endif

    if(argc < 2)
    {
        std::cout << "Usage:\ncuttlefish <command> [OPTIONS]" << std::endl;
        std::cout << "Supported commands: `build` and `validate`." << std::endl;
    }
    else
    {
        std::string command(argv[1]);
        std::transform(command.begin(), command.end(), command.begin(), [](const char ch) { return std::tolower(ch); });

        if(command == "build")
            build(argc - 1, argv + 1);
        else if(command == "validate")
            validate(argc - 1, argv + 1);
        else
            std::cout << "Invalid command. Supported commands: `build` and `validate`" << std::endl;
    }

    return EXIT_SUCCESS;
}