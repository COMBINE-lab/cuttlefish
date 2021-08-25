
#include "Directed_Kmer.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"
#include "Kmer_Buffered_Iterator.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "BBHash/BooPHF.h"
#include "Kmer_Hasher.hpp"
#include "Validator.hpp"
#include "Character_Buffer.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "FASTA_Record.hpp"
#include "kseq/kseq.h"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <chrono> 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <zlib.h>
#include <cstring>
#include <set>
#include <map>


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


// void check_N_base(const char* file_name, uint16_t k)
// {
//     Kmer_u64::set_k(k);

//     std::ifstream input(file_name, std::ifstream::in);
//     if(!input)
//     {
//         std::cerr << "Cannot open file " << file_name <<". Aborting.\n";
//         std::exit(EXIT_FAILURE);
//     }


//     std::string label;
//     uint64_t count;
//     uint64_t kmer_count = 0;
//     while(input >> label >> count)
//     {
//         kmer_count++;

//         if(label.find_first_of('N') != std::string::npos)
//             std::cout << "N nucleotide found in " << label << "\n";

//         if(kmer_count % 10000000 == 0)
//             std::cout << "Kmer-count: " << kmer_count << "\n";
//         // if(label.find_first_of("N") != std::string::npos)
//         //     std::cout << "N is present\n";
//     }

//     std::cout << "\nDone k-mers coding.\n"; 


//     input.close();
// }


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


template <uint16_t k>
void test_kmer_iterator(const char* file_name)
{
    const std::string kmc_file(file_name);

    // Open the k-mers container.
    Kmer_Container<k> kmers(kmc_file);

    // Kmer_u64::set_k(kmers.kmer_length());
    Kmer<k>::set_k(kmers.kmer_length());

    std::cout << "k-mers length: " << kmers.kmer_length() << "\n";
    std::cout << "k-mers count: " << kmers.size() << "\n";


    // Iterating over the k-mers on the database from disk.

    std::cout << "\nPerforming some non-trivial task with dereferenced iterators.\n";
    auto it_beg = kmers.begin();
    auto it_end = kmers.end();
    uint64_t count = 0;
    Kmer<k> max_kmer = Kmer<k>();
    for(auto it = it_beg; it != it_end; ++it)
    {
        // Use the iterator from here
        // std::cout << *it << "\n";
        max_kmer = std::max(max_kmer, *it);
        count++;
        if(count % 100000000 == 0)
            std::cout << "Processed " << count << " k-mers\n";
    }

    std::cout << "Max k-mer: " << max_kmer.string_label() << "\n";
    std::cout << "k-mers count found using iterators: " << count << "\n";
}


/*
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
    boophf_t * bphf = new boomphf::mphf<uint64_t, hasher_t>(input_keys.size(), input_keys, ".", thread_count);
    delete bphf;
}
*/


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


template <uint16_t k>
void test_buffered_iterator_performance(const char* const file_name)
{
    const std::string kmc_file(file_name);

    // Open the k-mers container.
    Kmer_Container<k> kmers(kmc_file);

    std::cout << "k-mers length: " << kmers.kmer_length() << "\n";
    std::cout << "k-mers count: " << kmers.size() << "\n";


    // Iterating over the k-mers on the database from disk.

    std::cout << "\nPerforming some non-trivial task with dereferenced buffered iterators.\n";
    auto it_beg(kmers.buf_begin());
    auto it_end(kmers.buf_end());
    uint64_t count = 0;
    Kmer<k> max_kmer = Kmer<k>();

    for(auto it = it_beg; it != it_end; ++it)
    {
        // Use the iterator from here
        max_kmer = std::max(max_kmer, *it);
        count++;
        if(count % 100000000 == 0)
            std::cout << "Processed " << count << " k-mers\n";
    }

    std::cout << "Max k-mer: " << max_kmer.string_label() << "\n";
    std::cout << "k-mers count found using iterators: " << count << "\n";
}


template <uint16_t k>
void test_SPMC_iterator_performance(const char* const db_path, const size_t consumer_count)
{
    Kmer_Container<k> kmer_container(db_path);

    // Kmer_SPMC_Iterator<k> it(kmer_container.spmc_begin(consumer_count));
    Kmer_SPMC_Iterator<k> it(&kmer_container, consumer_count);
    it.launch_production();

    std::cout << "\nProduction ongoing\n";

    std::vector<std::unique_ptr<std::thread>> T(consumer_count);
    std::vector<Kmer<k>> max_kmer(consumer_count);

    std::atomic<uint64_t> ctr{0};
    for(size_t i = 0; i < consumer_count; ++i)
    {
        const size_t consumer_id = i;
        auto& mk = max_kmer[consumer_id];
        T[consumer_id].reset(
            new std::thread([&kmer_container, &it, &mk, &ctr, consumer_id]()
            // new std::thread([&kmer_container, &it, &max_kmer, i]()
                {
                    std::cout << "Launched consumer " << consumer_id << ".\n";
                    Kmer<k> kmer;
                    Kmer<k> max_kmer;
                    // uint64_t local_count{0};
                    while(it.tasks_expected(consumer_id))
                        if(it.value_at(consumer_id, kmer))
                        {
                            max_kmer = std::max(max_kmer, kmer);
                            // local_count++;
                            // if (local_count % 5000000 == 0) {
                            //     ctr += local_count;
                            //     local_count = 0;
                            //     std::cerr << "parsed " << ctr << " k-mers\n";
                            // }
                        }

                    // ctr += local_count;
                    mk = max_kmer;
                }
            )
        );
    }


    it.seize_production();
    for(size_t i = 0; i < consumer_count; ++i)
        T[i]->join();

    //Kmer<k> global_max;
    //for (size_t i = 0; i < consumer_count; ++i) {
    //    global_max = std::max(global_max, max_kmer[i]);
    //}
    std::cout << "Parsed " << ctr << " k-mers\n";
    std::cout << "Max k-mer: " << std::max_element(max_kmer.begin(), max_kmer.end())->string_label() << "\n";
}


template <uint16_t k>
void test_iterator_correctness(const char* const db_path, const size_t consumer_count)
{
    const std::string kmc_path(db_path);

    // Open the k-mers container.
    Kmer_Container<k> kmer_container(kmc_path);

    std::cout << "k-mers length: " << kmer_container.kmer_length() << "\n";
    std::cout << "k-mers count: " << kmer_container.size() << "\n";


    // Iterating over the k-mers on the database from disk.

    std::cout << "\nCollecting the k-mers with a buffered iterator.\n";

    std::vector<Kmer<k>> buf_kmers;
    buf_kmers.reserve(kmer_container.size());

    auto it_end(kmer_container.buf_end());
    uint64_t count = 0;
    for(auto it(kmer_container.buf_begin()); it != it_end; ++it)
    {
        // Use the iterator from here
        buf_kmers.emplace_back(*it);
        count++;
        if(count % 100000000 == 0)
            std::cout << "Processed " << count << " k-mers\n";
    }


    std::cout << "\nCollecting the k-mers with an SPMC iterator.\n";
    std::vector<std::vector<Kmer<k>>> kmers(consumer_count);
    Kmer_SPMC_Iterator<k> it(kmer_container.spmc_begin(consumer_count));

    auto start = std::chrono::high_resolution_clock::now(); 
    it.launch_production();

    std::vector<std::unique_ptr<std::thread>> T(consumer_count);

    for(size_t i = 0; i < consumer_count; ++i)
        T[i].reset(
            new std::thread([&it, &kmers, i]()
                {
                    std::cout << "Launched consumer " << i << ".\n";
                    Kmer<k> kmer;
                    
                    while(it.tasks_expected(i))
                        if(it.task_available(i) && it.value_at(i, kmer))
                            kmers[i].emplace_back(kmer);
                }
            )
        );


    it.seize_production();
    for(size_t i = 0; i < consumer_count; ++i)
        T[i]->join();

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    // To get the value of duration use the count() 
    // member function on the duration object 
    std::cout << "actual parsing took: " << duration.count() << std::endl; 

    std::vector<Kmer<k>> spmc_kmers;
    spmc_kmers.reserve(kmer_container.size());
    for(size_t i = 0; i < consumer_count; ++i)
        spmc_kmers.insert(spmc_kmers.end(), kmers[i].begin(), kmers[i].end());

    
    if(buf_kmers.size() != spmc_kmers.size())
        std::cout << "Incorrect set of k-mers\n";
    else
    {
        std::cout << "k-mer counts match\n";

        std::cout << "Sorting the k-mer sets (of both collection).\n";
        std::sort(buf_kmers.begin(), buf_kmers.end());
        std::sort(spmc_kmers.begin(), spmc_kmers.end());
        std::cout << "Done sorting\n";

        size_t mis = 0;
        for(size_t i = 0; i < buf_kmers.size(); ++i)
            if(!(buf_kmers[i] == spmc_kmers[i]))
            {
                // std::cout << "Mismatching k-mers found\n";
                mis++;
            }

        std::cout << "#mismatching_k-mers = " << mis << "\n";
        std::cout << (mis > 0 ? "Incorrect" : "Correct") << " k-mer set parsed.\n";
    }
}


template <uint16_t k>
void write_kmers(const std::string& kmc_db_path, const uint16_t thread_count, const std::string& output_file_path)
{
    const Kmer_Container<k> kmer_container(kmc_db_path);
    Kmer_SPMC_Iterator<k> parser(&kmer_container, thread_count);

    parser.launch_production();
    
    std::ofstream output(output_file_path);

    std::vector<std::unique_ptr<std::thread>> T(thread_count);

    for(size_t i = 0; i < thread_count; ++i)
    {
        const size_t consumer_id = i;
        
        T[consumer_id].reset(
            new std::thread([&parser, consumer_id, &output]()
                {
                    Kmer<k> kmer;
                    std::vector<char> str;
                    str.reserve(k + 2);

                    uint64_t local_count{0};
                    Character_Buffer<10485760, std::ofstream> buffer(output);

                    while(parser.tasks_expected(consumer_id))
                        if(parser.value_at(consumer_id, kmer))
                        {
                            kmer.get_label(str);
                            str.emplace_back('\n');
                            // buffer += str;
                            buffer += FASTA_Record<std::vector<char>>(0, str);

                            local_count++;
                            if(local_count % 10000000 == 0)
                                std::cout << "Thread " << consumer_id << " parsed " << local_count << " k-mers\n";
                        }
                }
            )
        );
    }


    parser.seize_production();
    for(std::size_t id = 0; id < thread_count; ++id)
        T[id]->join();

    output.close();
}


int main(int argc, char** argv)
{
    (void)argc;
    (void)argv;
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

    // test_kmer_iterator<61>(argv[1]);

    // test_async_writer(argv[1]);

    // count_kmers_in_unitigs(argv[1], atoi(argv[2]));

    static constexpr uint16_t k = 28;
    static const size_t consumer_count = std::atoi(argv[2]);

    // test_buffered_iterator_performance<k>(argv[1]);
    test_SPMC_iterator_performance<k>(argv[1], consumer_count);

    // write_kmers<32>(argv[1], std::atoi(argv[2]), argv[3]);

    return 0;
}
