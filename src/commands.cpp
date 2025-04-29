#include "Input_Defaults.hpp"
#include "dBG_Contractor.hpp"
// #include "Kmer_Index.hpp"
#include "Build_Params.hpp"
#include "Application.hpp"
#include "version.hpp"
#include "cxxopts/cxxopts.hpp"
#include "parlay/parallel.h"

#include <string>
#include <vector>
#include <iostream>
#include <optional>

#ifdef __cplusplus
extern "C" {
#endif
  int cf_build(int argc, char** argv);
  int cf_validate(int argc, char** argv);
  int print_cf_version();
#ifdef __cplusplus
}
#endif


// Driver function for the CdBG build.
int cf_build(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish build", "Efficiently construct the compacted de Bruijn graph from sequencing reads or reference sequences");

    std::optional<std::vector<std::string>> seqs;
    std::optional<std::vector<std::string>> lists;
    std::optional<std::vector<std::string>> dirs;
    std::optional<uint32_t> cutoff;
    options.add_options("common")
        ("s,seq", "input files",
            cxxopts::value<std::optional<std::vector<std::string>>>(seqs))
        ("l,list", "input file lists",
            cxxopts::value<std::optional<std::vector<std::string>>>(lists))
        ("d,dir", "input file directories",
            cxxopts::value<std::optional<std::vector<std::string>>>(dirs))
        ("k,kmer-len", "k-mer length",
            cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::K)))
        // ("t,threads", "number of threads to use",
        //     cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::THREAD_COUNT)))
        // ("idx", "construct a k-mer index of the de Bruijn graph")
        ("min-len", "minimizer length",
            cxxopts::value<uint16_t>()->default_value(std::to_string(cuttlefish::_default::MIN_LEN)))
        ("o,output", "output file",
            cxxopts::value<std::string>())
        ("w,work-dir", "working directory",
            cxxopts::value<std::string>()->default_value(cuttlefish::_default::WORK_DIR))
/*
    options.add_options("cuttlefish_2")
        ("path-cover", "extract a maximal path cover of the de Bruijn graph")
        ;
*/

    // options.add_options("cuttlefish_3")
        ("read", "construct a compacted read de Bruijn graph (for FASTQ input)")
        ("ref", "construct a compacted reference de Bruijn graph (for FASTA input)")
        ("c,cutoff", "frequency cutoff for (k + 1)-mers (default: refs: " + std::to_string(cuttlefish::_default::CUTOFF_FREQ_REFS) + ", reads: " + std::to_string(cuttlefish::_default::CUTOFF_FREQ_READS) + ")",
            cxxopts::value<std::optional<uint32_t>>(cutoff))
        ("color", "whether to color the compacted graph or not")
/*
        ("vertex-part-count", "number of vertex-partitions in the discontinuity graph; needs to be a power of 2",
            cxxopts::value<std::size_t>()->default_value(std::to_string(cuttlefish::_default::VERTEX_PART_COUNT)))
        ("lmtig-bucket-count", "number of buckets storing literal locally-maximal unitigs",
            cxxopts::value<std::size_t>()->default_value(std::to_string(cuttlefish::_default::LMTIG_BUCKET_COUNT)))
        ("gmtig-bucket-count", "number of buckets for global maximal unitigs",
            cxxopts::value<std::size_t>()->default_value(std::to_string(cuttlefish::_default::GMTIG_BUCKET_COUNT)))
*/
        ("h,help", "print usage")
        ;

    std::optional<uint16_t> format_code;
/*
    options.add_options("cuttlefish_1")
        ("f,format", "output format (0: FASTA, 1: GFA 1.0, 2: GFA 2.0, 3: GFA-reduced)",
            cxxopts::value<std::optional<uint16_t>>(format_code))
        ("track-short-seqs", "track existence of sequences shorter than k bases")
        ("poly-N-stretch", "includes information of polyN stretches in the tiling output")
        ;
*/
        ;

    try
    {
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0;
        }

        const auto is_read_graph = result["read"].as<bool>();
        const auto is_ref_graph = result["ref"].as<bool>();
        const auto k = result["kmer-len"].as<uint16_t>();
        const auto color = result["color"].as<bool>();
        const auto vertex_part_count = cuttlefish::_default::VERTEX_PART_COUNT;     // result["vertex-part-count"].as<std::size_t>();
        const auto lmtig_bucket_count = cuttlefish::_default::LMTIG_BUCKET_COUNT;   // result["lmtig-bucket-count"].as<std::size_t>();
        const auto gmtig_bucket_count = cuttlefish::_default::GMTIG_BUCKET_COUNT;   // result["gmtig-bucket-count"].as<std::size_t>();
        // const auto thread_count = result["threads"].as<uint16_t>();
        const auto idx = false; // result["idx"].as<bool>();
        const auto min_len = result["min-len"].as<uint16_t>();
        const auto output_file = result["output"].as<std::string>();
        const auto format = format_code ?   std::optional<cuttlefish::Output_Format>(cuttlefish::Output_Format(format_code.value())) :
                                            std::optional<cuttlefish::Output_Format>();
        const auto track_short_seqs = false;    // result["track-short-seqs"].as<bool>();
        const auto poly_n_stretch = false;  // result["poly-N-stretch"].as<bool>();
        const auto working_dir = result["work-dir"].as<std::string>();
        const auto path_cover = false;  // result["path-cover"].as<bool>();

        const Build_Params params(  is_read_graph, is_ref_graph,
                                    seqs, lists, dirs,
                                    k, cutoff,
                                    color,
                                    vertex_part_count, lmtig_bucket_count, gmtig_bucket_count,
                                    parlay::num_workers(),
                                    idx, min_len,
                                    output_file, format, track_short_seqs, poly_n_stretch, working_dir,
                                    path_cover
                                );
        if(!params.is_valid())
        {
            std::cerr << "Invalid input configuration. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        // std::cout.precision(3);


        if(params.idx())
        {
            std::cout << "\nConstructing a k-mer index of the de Bruijn graph for k = " << k << ".\n";

            // Application<cuttlefish::MAX_K, Kmer_Index>(params).execute();

            std::cout << "\nConstructed a k-mer index of the de Bruijn graph at " << params.output_prefix() << ".\n";
        }
        else
        {
            const std::string dBg_type(params.is_read_graph() ? "read" : "reference");

            std::cout << "\nConstructing the compacted " << dBg_type << " de Bruijn graph for k = " << k << ".\n";

            // (params.is_read_graph() || params.is_ref_graph()) ?
            //     Application<cuttlefish::MAX_K, Read_CdBG>(params).execute() :
            //     Application<cuttlefish::MAX_K, CdBG>(params).execute();
            Application<cuttlefish::MAX_K, cuttlefish::dBG_Contractor>(params).execute();

            std::cout << "\nConstructed the " << dBg_type << " compacted de Bruijn graph at " << output_file << ".\n";
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << std::endl << "Usage :" << std::endl;
        std::cerr << options.help() << std::endl;
    }

    return EXIT_SUCCESS;
}


// Driver function for the CdBG validation.
/*
int cf_validate(int argc, char** argv)
{
    cxxopts::Options options("cuttlefish validate", "Validate a compacted de Bruijn graph constructed by cuttlefish");
    options.add_options()
        ("r,refs", "reference files",
            cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("l,lists", "reference file lists",
            cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("d,dirs", "reference file directories",
            cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("k,kmer_len", "k-mer length",
            cxxopts::value<uint16_t>())
        ("s,kmc_db", "set of k-mers (KMC database) prefix",
            cxxopts::value<std::string>())
        ("g,cdbg", "compacted de Bruijn graph file",
            cxxopts::value<std::string>())
        ("t,threads", "number of threads to use",
            cxxopts::value<uint16_t>()->default_value("1"))
        ("w,work_dir", "working directory",
            cxxopts::value<std::string>()->default_value("."))
        ("mph", "minimal perfect hash (BBHash) file (optional)",
            cxxopts::value<std::string>()->default_value(""))
        ("h,help", "print usage");

    try
    {
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0;
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

    return 0;
}
*/

int print_cf_version()
{
    std::cout << "cuttlefish " VERSION << std::endl;
    return 0;
}