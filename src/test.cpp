
#include "Kmer.hpp"
#include "Directed_Kmer.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"
#include "BBHash/BooPHF.h"
#include "Kmer_Hasher.hpp"
#include "Validator.hpp"
#include "kseq/kseq.h"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <cstring>
#include <set>
#include <map>

//#include "Validator.cpp"


// Debug
// #include "Validator.cpp"


// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(int, read);


void test_kseq(const char* fileName)
{
    // STEP 2: open the file handler
    FILE* fp = fopen(fileName, "r");

    // STEP 3: initialize seq
    kseq_t* parser = kseq_init(fileno(fp));


    // STEP 4: read sequence
    int c = 0;
    size_t max_len = 0, max_size = 0;
    size_t len = 0;
    while(true)
    {
        if(kseq_read(parser) < 0)
            break;

        c++;

        std::cout << "Name: " << parser->name.s << "\n";
        if(parser->comment.l)
            std::cout << "Comment: " << parser->comment.s << "\n";

        // std::cout << "Seq: " << seq->seq.s << "\n";
        std::cout << "seq.m: " << parser->seq.m << "\n";
        std::cout << "seq.l: " << parser->seq.l << "\n";

        len += parser->seq.l;
        max_len = std::max(max_len, parser->seq.l);
        max_size = std::max(max_size, parser->seq.m);
        
        
        if(parser->qual.l)
            std::cout << "Quality: " << parser->qual.s << "\n";
    }


    std::cout << "Line count: " << c << "\n";
    std::cout << "Max seq length: " << max_len << "\n";
    std::cout << "Max seq buffer size: " << max_size << "\n";
    std::cout << "Total reference length: " << len << "\n";

    kseq_destroy(parser);
    fclose(fp);
}


void check_repeated_kmers(const char* file_name)
{
    std::ifstream input(file_name, std::ifstream::in);
    if(!input)
    {
        std::cerr << "Cannot open file " << file_name <<". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::set<uint64_t> kmers;
    uint64_t kmer;
    uint64_t kmer_count = 0;
    while(input >> kmer)
    {
        if(kmers.find(kmer) != kmers.end())
        {
            std::cerr << "Repeated k-mers found. Aborting\n";
            std::exit(EXIT_FAILURE);
        }

        kmer_count++;
        if(kmer_count % 10000000 == 0)
            std::cout << "\rKmer-count: " << kmer_count;


        kmers.insert(kmer);
    }


    std::cout << "\nNo repeated k-mers found.\n";

    input.close();
}


void check_N_base(const char* file_name, uint16_t k)
{
    Kmer::set_k(k);

    std::ifstream input(file_name, std::ifstream::in);
    if(!input)
    {
        std::cerr << "Cannot open file " << file_name <<". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::string label;
    uint64_t count;
    uint64_t kmer_count = 0;
    while(input >> label >> count)
    {
        kmer_count++;

        if(label.find_first_of('N') != std::string::npos)
            std::cout << "N nucleotide found in " << label << "\n";

        if(kmer_count % 10000000 == 0)
            std::cout << "Kmer-count: " << kmer_count << "\n";
        // if(label.find_first_of("N") != std::string::npos)
        //     std::cout << "N is present\n";
    }

    std::cout << "\nDone k-mers coding.\n"; 


    input.close();
}


void check_absent_nucleotide(const char* file_name)
{
    // STEP 2: open the file handler
    FILE* fp = fopen(file_name, "r");

    // STEP 3: initialize seq
    kseq_t* parser = kseq_init(fileno(fp));


    // STEP 4: read sequence
    int c = 0;
    while(kseq_read(parser) >= 0)
    {
        c++;

        std::cout << "Name: " << parser->name.s << "\n";
        if(parser->comment.l)
            std::cout << "Comment: " << parser->comment.s << "\n";

        if(parser->qual.l)
            std::cout << "Quality: " << parser->qual.s << "\n";

        std::cout << "Length: " << parser->seq.l << "\n";

        const char* seq = parser->seq.s;
        for(size_t i = 0; seq[i];)
            if(seq[i] == 'N')
            {
                size_t start_idx = i;
                
                while(seq[i] == 'N')
                    i++;

                size_t end_idx = i - 1;

                std::cout << "[" << start_idx << ", " << end_idx << "]: " << end_idx - start_idx + 1 << "\n";
            }
            else
                i++;

        
        // break;
    }

    kseq_destroy(parser);
    fclose(fp);
}


void test_unipaths(const char* file_name)
{
    std::ifstream input(file_name);
    if(!input)
    {
        std::cerr << "Error opening file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::string unitig;
    std::map<size_t, uint64_t> L;
    std::set<std::string> U;
    while(input >> unitig)
    {
        if(U.find(unitig) != U.end())
            std::cerr << "Repitition of the same unipath\n";
        else
            U.insert(unitig);

        L[unitig.length()]++;
    }


    input.close();

    
    std::cout << "Total unitigs: " << U.size() << "\n";
    for(auto p: L)
        std::cout << p.first << " " << p.second << "\n";
}


void test_kmer_iterator(const char* file_name)
{
    const std::string kmc_file(file_name);

    // Open the k-mers container.
    Kmer_Container kmers(kmc_file);

    Kmer::set_k(kmers.kmer_length());

    std::cout << "k-mers length: " << kmers.kmer_length() << "\n";
    std::cout << "k-mers count: " << kmers.size() << "\n";


    // Iterating over the k-mers on the database from disk.

    std::cout << "\nPerforming some non-trivial task with dereferenced iterators.\n";
    auto it_beg = kmers.begin();
    auto it_end = kmers.end();
    uint64_t count = 0;
    cuttlefish::kmer_t max_kmer = cuttlefish::kmer_t();
    for(auto it = it_beg; it != it_end; ++it)
    {
        // Use the iterator from here
        // std::cout << *it << "\n";
        max_kmer = std::max(max_kmer, *it);
        count++;
    }

    std::cout << "Max k-mer: " << max_kmer.string_label() << "\n";
    std::cout << "k-mers count found using iterators: " << count << "\n";
}


void check_uint64_BBHash(const char* file_name, uint16_t thread_count)
{
    typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
    typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

    std::vector<uint64_t> input_keys;
    std::ifstream input;
    input.open(file_name);
    uint64_t key;
    std::cout << "Loading k-mers.\n";
    while(input >> key)
        input_keys.emplace_back(key);
    std::cout << "Loading k-mers done.\n";
    input.close();
    //build the mphf  
    boophf_t * bphf = new boomphf::mphf<uint64_t, hasher_t>(input_keys.size(), input_keys, thread_count);
    delete bphf;
}


void test_async_writer(const char* log_file_name)
{
    // Clear the log file first, as `spdlog` logger appends messages.
    std::ofstream temp(log_file_name);
    temp.close();


    auto f = [](uint64_t thread_ID, std::shared_ptr<spdlog::logger> file_writer)
    {
        for(int i = 0; i < 100; ++i)
            file_writer->info("Writing {} from thread {}", i, thread_ID);
    };

    try
    {
        auto async_file = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", log_file_name);

        // Set log message pattern for the writer.
        async_file->set_pattern("%v");


        std::vector<std::thread> writer;
        for(int i = 0; i < 5; ++i)
            writer.emplace_back(f, i, async_file);

        for(int i = 0; i < 5; ++i)
            writer[i].join();


        // Close the loggers?
        spdlog::drop_all();
    }
    catch(const spdlog::spdlog_ex& ex)
    {
        std::cerr << "Logger initialization failed with: " << ex.what() << "\n";
    }
}


void count_kmers_in_unitigs(const char* file_name, uint16_t k)
{
    std::ifstream input(file_name);
    if(!input)
    {
        std::cerr << "Error opening file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    
    uint64_t kmer_count = 0;
    std::string unitig;
    while(input >> unitig)
    {
        if(unitig.length() < k)
        {
            std::cerr << "Unitig length " << unitig.length() << " found, < k-value.\n";
            std::exit(EXIT_FAILURE);
        }

        kmer_count += unitig.length() - k + 1;
    }


    std::cout << "Total k-mers found: " << kmer_count << "\n";
}


void validate_result(char **argv)
{
    const std::string ref_file_name(argv[1]);
    const uint16_t k(atoi(argv[2]));
    const std::string kmc_db_name(argv[3]);
    const std::string cdbg_file_name(argv[4]);
    const std::string bbhash_file_name(argv[5]);
    const uint16_t thread_count(atoi(argv[6]));
    cuttlefish::logger_t console = spdlog::stdout_color_mt("Validator");
    
    Validator validator(ref_file_name, k, kmc_db_name, cdbg_file_name, console);
    std::cout << "Validation " << (validator.validate(bbhash_file_name, thread_count) ? "successful" : "failed") << "\n";
}


int main(int argc, char** argv)
{
    (void)argc;
    // const char* fileName = argv[1];

    // test_kseq(argv[1]);

    // test_rolling_kmer(argv[1], atoi(argv[2]));

    // test_kmer_sampling(std::string(argv[1]), atoi(argv[2]), std::string(argv[3]), std::atoi(argv[4]));

    // convert_kmers_to_int(argv[1], atoi(argv[2]), argv[3]);

    // check_repeated_kmers(argv[1]);

    // check_N_base(argv[1], atoi(argv[2]));

    // check_absent_nucleotide(argv[1]);

    // test_unipaths(argv[1]);

    // check_uint64_BBHash(argv[1], atoi(argv[2]));

    // test_kmer_iterator(argv[1]);

    // test_async_writer(argv[1]);

    // count_kmers_in_unitigs(argv[1], atoi(argv[2]));

    validate_result(argv);


    return 0;
}
