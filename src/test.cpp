
/*
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(int, read)


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
                    uint64_t local_count{0};
                    while(it.tasks_expected(consumer_id))
                        if(it.value_at(consumer_id, kmer))
                        {
                            max_kmer = std::max(max_kmer, kmer);
                            local_count++;
                            if (local_count % 10000000 == 0) {
                                ctr += local_count;
                                local_count = 0;
                                std::cerr << "\rparsed " << ctr << " k-mers";
                            }
                        }

                    ctr += local_count;
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
    std::cout << "\nParsed " << ctr << " k-mers\n";
    std::cout << "Max k-mer: " << std::max_element(max_kmer.begin(), max_kmer.end())->string_label() << "\n";
}


template <uint16_t k>
void test_SPSC_iterator_performance(const char* const db_path)
{
    Kmer_Container<k> kmer_container(db_path);
    std::cout << "kmer_count: " << kmer_container.size() << "\n";

    Kmer_SPSC_Iterator<k> it(&kmer_container);
    it.launch();

    std::cout << "\nIteration ongoing\n";

    Kmer<k> max_kmer;
    // Kmer<k> kmer;
    std::vector<Kmer<k>> kmers;
    uint64_t kmer_count = 0;

    // while(it.parse_kmer(kmer))
    while(it.parse_kmers_atomic(kmers))
    {
        // if(max_kmer < kmer)
        //     max_kmer = kmer;
        max_kmer = std::max(max_kmer, *std::max_element(kmers.begin(), kmers.end()));

        // kmer_count++;
        kmer_count += kmers.size();
        // if(kmer_count % 100000 == 0)
        //     std::cout << "\rParsed " << kmer_count << " k-mers.";
    }

    it.close();

    std::cout << "\rParsed " << kmer_count << " k-mers.\n";
    std::cout << "\nMax k-mer: " << max_kmer.string_label() << "\n";
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


template <uint16_t k>
void test_multiway_merge(const std::string& db_list)
{
    std::vector<std::string> dbs;

    std::ifstream input(db_list);
    std::string path;
    while(input >> path)
        dbs.emplace_back(path);

    input.close();

    std::cout << "DB-count: " << dbs.size() << "\n";

    Multiway_Merger<k> merger(dbs);
    merger.launch();

    Kmer<k> kmer;
    std::vector<uint16_t> color;
    std::vector<uint16_t> empty_sources;
    uint32_t kmer_count = 0;
    Kmer<k> last_kmer;

    while(merger.template next<std::vector<uint16_t>>(kmer, color, empty_sources))
    {
        if(!empty_sources.empty())
        {
            merger.top_up(empty_sources);
            empty_sources.clear();
        }

        kmer_count++;

        if(kmer_count % 1000000 == 0)
            std::cerr << "\r" << kmer_count;

        // if(kmer < last_kmer)
        // {
        //     std::cout << "k-mers found not in sorted order. Aborting.\n";
        //     std::exit(EXIT_FAILURE);
        // }

        last_kmer = kmer;
    }

    std::cout << "k-mer count: " << kmer_count << "\n";
}


template <uint16_t k, uint16_t l>
void minimizer_iterator_test(const char* const file_path)
{
    Ref_Parser input(file_path);

    while(input.read_next_seq())
    {
        std::cout << "\nAt seq. " << input.seq_id() << ", of length " << input.seq_len() << ".\n";

        if(input.seq_len() < k)
            continue;

        const char* const seq = input.seq();
        Minimizer_Iterator min_iterator(seq, input.seq_len(), k, l);

        std::size_t kmer_idx = 0;
        cuttlefish::minimizer_t min;
        std::size_t idx;

        // std::cout << seq << "\n";

        do
        {
            std::cerr << "\r" << kmer_idx;

            // Get minimizer from iterator.
            min_iterator.value_at(min, idx);
            const uint64_t hash = Minimizer_Iterator::hash(min);

            // std::cout << "k-mer idx: " << kmer_idx << ", minimizer idx: " << idx << "\n";

            // Get minimizer with complete search.
            Kmer<l> min_lmer(seq, kmer_idx);
            std::size_t min_idx = kmer_idx;
            uint64_t min_hash = Minimizer_Iterator::hash(min_lmer.as_int());

            for(std::size_t i = kmer_idx + 1; i + l - 1 < kmer_idx + k; ++i)
            {
                const Kmer<l> lmer(seq, i);
                const uint64_t lmer_hash = Minimizer_Iterator::hash(lmer.as_int());

                if(min_hash < lmer_hash)
                    continue;

                if(min_hash > lmer_hash || min_lmer > lmer)
                    min_lmer = lmer, min_idx = i, min_hash = lmer_hash;

                // Basic lexicographic minimizer.
                // if(min_lmer > lmer)
                //     min_lmer = lmer, min_idx = i;
            }

            if(min_lmer.as_int() != min || min_idx != idx || min_hash != hash)
            {
                std::cout << "\nBug in computing minimizers. Aborting.\n";
                std::cout << "k-mer index: " << kmer_idx << ".\n";
                std::cout << "From iterator: " << min << ", at index " << idx << ".\n";
                std::cout << "Manually:      " << min_lmer.as_int() << ", at index " << min_idx << ".\n";

                std::exit(EXIT_FAILURE);
            }


            kmer_idx++;
        }
        while(++min_iterator);
    }


    input.close();
}


template <uint16_t k>
void cross_check_terminal_kmers(const std::string& path_pref)
{
    Kmer_Index<k> kmer_idx(path_pref);

    std::vector<Kmer<k>> from_idx, from_fa;

    Kmer<k> kmer;
    for(std::size_t id = 0; id < kmer_idx.path_count(); ++id)
    {
        kmer_idx.get_kmer(id, 0, kmer);
        from_idx.emplace_back(kmer);

        if(kmer_idx.path_size(id) > k)
        {
            kmer_idx.get_kmer(id, kmer_idx.path_size(id) - k, kmer);
            from_idx.emplace_back(kmer);
        }
    }

    std::cout << "Collected terminal k-mers from the k-mer index.\n";


    std::ifstream input(path_pref + ".fa");
    std::string path;
    while(input >> path)
    {
        if(path.front() == '>')
            continue;

        from_fa.emplace_back(path, 0);

        if(path.length() > k)
            from_fa.emplace_back(path, path.length() - k);
    }

    std::cout << "Collected terminal k-mers from the CdBG output.\n";

    std::cout << "# from k-mer idx: " << from_idx.size() << "\n";
    std::cout << "# from CdBG o/p:  " << from_fa.size() << "\n";


    std::sort(from_idx.begin(), from_idx.end());
    std::sort(from_fa.begin(), from_fa.end());
    std::cout << "Terminal k-mers " << (from_idx == from_fa ? "matched" : " doesn't match") << "\n";
}


template <uint16_t k>
void cross_check_index_kmers(const std::string& op_pref)
{
    Kmer_Index<k> kmer_idx(op_pref);

    std::vector<Kmer<k>> from_idx, from_fa;

    Kmer<k> kmer;
    for(std::size_t id = 0; id < kmer_idx.path_count(); ++id)
    {
        const std::size_t path_size = kmer_idx.path_size(id);
        for(std::size_t idx = 0; idx + k <= path_size; ++idx)
        {
            kmer_idx.get_kmer(id, idx, kmer);
            from_idx.emplace_back(kmer);
        }
    }

    std::cout << "Collected all k-mers from the k-mer index.\n";


    std::ifstream input(op_pref + ".fa");
    std::string path;
    while(input >> path)
    {
        if(path.front() == '>')
            continue;

        for(std::size_t idx = 0; idx + k <= path.length(); ++idx)
            from_fa.emplace_back(path, idx);
    }

    std::cout << "Collected all k-mers from the CdBG output.\n";

    std::cout << "# from k-mer idx: " << from_idx.size() << "\n";
    std::cout << "# from CdBG o/p:  " << from_fa.size() << "\n";


    std::sort(from_idx.begin(), from_idx.end());
    std::sort(from_fa.begin(), from_fa.end());
    std::cout << "All k-mers " << (from_idx == from_fa ? "matched" : " doesn't match") << "\n";
}


template <uint16_t k, uint16_t l>
uint64_t compute_breakpoints(const std::string& file_path)
{
    std::ifstream input(file_path);

    uint64_t br = 0;
    uint64_t min_inst = 0;
    uint64_t unitig_c = 0;
    std::string unitig;
    while(input >> unitig)
    {
        if(unitig.front() == '>')
            continue;

        // Minimizer_Iterator min_it(unitig.c_str(), unitig.length(), k - 1, l);
        Minimizer_Iterator<const char*, true> min_it(unitig.c_str(), unitig.length(), k - 1, l);
        cuttlefish::minimizer_t last_min, min, last_idx, idx;

        min_it.value_at(last_min, last_idx);
        min_inst++;
        while(++min_it)
        {
            min_it.value_at(min, idx);
            if(min != last_min)
                br++, min_inst++, last_min = min, last_idx = idx;
            else if(idx != last_idx)
                min_inst++, last_idx = idx;
        }

        unitig_c++;
        if(unitig_c % 1000000 == 0)
            std::cerr <<    "\rProcessed " << unitig_c << " unitigs. "
                            "#minimizer-instances: " << min_inst << ", #breakpoints: " << br << ".";
    }


    std::cout << "\n#minimizer-instances: " << min_inst << ", #breakpoints: " << br << ".\n";

    return br;
}


template <uint16_t k, uint16_t l>
bool validate_canonical_minimizer(const std::string& file_path)
{
    Ref_Parser parser(file_path);
    uint64_t matched = 0;
    std::size_t parsed_len = 0;

    while(parser.read_next_seq())
    {
        const char* p = parser.seq();

        Minimizer_Iterator<decltype(p), true> min_it(p, parser.seq_len(), k, l);
        cuttlefish::minimizer_t min;
        std::size_t min_idx;
        cuttlefish::minimizer_t min_manual;

        do
        {
            min_it.value_at(min, min_idx);
            min_manual = Minimizer_Utility::canonical_minimizer(Kmer<k>(p), l);

            if(min != min_manual)
            {
                std::cerr << "mismatched at k-mer " << matched << "\n";
                std::cerr << "Iterator-hash: " << Minimizer_Utility::hash(min) << ", Manual-hash: " << Minimizer_Utility::hash(min_manual) << "\n";
                return false;
            }

            matched++;
            p++;
        }
        while(++min_it);

        parsed_len += parser.seq_len();
        std::cerr << "\rParsed length: " << parsed_len;
    }

    return true;
}


template <uint16_t k, uint16_t l>
void get_discontinuity_edges(const std::string& file_path)
{
    Ref_Parser parser(file_path);
    std::vector<cuttlefish::Discontinuity_Edge<k>> E;
    uint64_t disc_count = 0;
    uint64_t trivial_unitig_count = 0;
    std::size_t max_trivial_len = 0;

    while(parser.read_next_seq())
    {
        const auto seq = parser.seq();
        Minimizer_Iterator<decltype(seq), true> min_it(seq, parser.seq_len(), k - 1, l);
        cuttlefish::minimizer_t last_min, last_idx, min, idx;

        std::size_t kmer_idx = 0;
        Kmer<k> last_disc_kmer, curr_disc_kmer;
        bool on_chain = false;

        min_it.value_at(last_min, last_idx);
        while(++min_it)
        {
            min_it.value_at(min, idx);
            if(idx != last_idx)
            {
                disc_count++;

                curr_disc_kmer = Kmer<k>(seq, kmer_idx);
                if(on_chain)
                {
                    const Kmer<k>& x = last_disc_kmer;
                    const Kmer<k>& y = curr_disc_kmer;
                    cuttlefish::side_t s_x = (x == x.canonical() ? cuttlefish::side_t::back : cuttlefish::side_t::front);
                    cuttlefish::side_t s_y = (y == y.canonical() ? cuttlefish::side_t::front : cuttlefish::side_t::back);

                    E.emplace_back(x, s_x, y, s_y, 1, 0);
                }

                last_disc_kmer = curr_disc_kmer;
                on_chain = true;
            }

            last_min = min, last_idx = idx;
            kmer_idx++;
        }


        if(!on_chain)
            trivial_unitig_count++,
            max_trivial_len = std::max(max_trivial_len, parser.seq_len());
    }

    std::cerr << "#discontinuity-k-mers: " << disc_count << "\n";
    std::cerr << "#discontinuity-edge:   " << E.size() << "\n";
    std::cerr << "#trivial-unitigs: " << trivial_unitig_count << "\n";
    std::cerr << "#max-trivial-op-len: " << max_trivial_len << "\n";
}


void test_unitig_file(const std::string cdbg_path, const std::string& op_path)
{
    Ref_Parser parser(cdbg_path);

    std::vector<std::string> U;
    std::size_t deposited_chars = 0;

    cuttlefish::Unitig_File_Writer uni_writer(op_path + ".written");

    while(parser.read_next_seq())
    {
        const auto seq = parser.seq();
        const auto seq_len = parser.seq_len();
        U.emplace_back(seq);
        deposited_chars += seq_len;
        uni_writer.add(seq, seq + seq_len);
    }

    uni_writer.close();


    cuttlefish::Unitig_File_Reader uni_reader(op_path + ".written");
    std::string unitig;
    std::ofstream output(op_path + ".read");

    std::size_t idx = 0;
    while(uni_reader.read_next_unitig(unitig))
    {
        if(U[idx] != unitig)
        {
            std::cerr << "Mismatch at unitig-index " << idx << "\n";
            std::cerr << "expected length: " << U[idx].length() << "; parsed length: " << unitig.length() << "\n";
            // std::cerr << "Expected " << U[idx] << ";\nParsed " << unitig << "\n";
            std::exit(EXIT_FAILURE);
        }

        output << unitig << "\n";
        idx++;
    }

    output.close();

    std::cerr << "All unitigs matched\n";

    // std::cerr << "Deposited " << deposited_chars << "; written " << uni_writer.size() << " chars.\n";
}


template <uint16_t k>
void bootstrap_discontinuity_graph(const uint16_t l, const std::string& cdbg_path, const std::string& dg_path, const std::size_t part_count, const std::size_t unitig_bucket_count)
{
    cuttlefish::Edge_Matrix<k> E(part_count, dg_path + "E_");
    cuttlefish::Discontinuity_Graph_Bootstrap<k> dgb(cdbg_path, l, E, dg_path, unitig_bucket_count);
    dgb.generate();
    E.serialize();
}


template <uint16_t k>
void benchmark_hash_table(std::size_t elem_count, double load_factor = 0.8)
{
    cuttlefish::Concurrent_Hash_Table<Kmer<k>, std::size_t, Kmer_Hasher<k>> ht(elem_count, load_factor);
    std::vector<Kmer<k>> kmers;
    std::vector<std::size_t> vals;

    kmers.reserve(elem_count);
    vals.reserve(elem_count);

    std::random_device random_device;
    std::mt19937 random_engine(random_device());
    std::uniform_int_distribution<std::size_t> distribution(0lu, std::numeric_limits<std::size_t>::max());

    for(std::size_t i = 0; i < elem_count; ++i)
        kmers.emplace_back(Kmer<k>(distribution(random_engine))),
        vals.emplace_back(distribution(random_engine));

    constexpr auto now = std::chrono::high_resolution_clock::now;
    constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };

    auto t_st = now();
    for(std::size_t i = 0; i < elem_count; ++i)
        ht.template insert<false>(kmers[i], vals[i]);
    auto t_en = now();

    std::cerr << "Inserted " << elem_count << " elements into custom hash-table in time " << duration(t_en - t_st) << " seconds.\n";


    std::unordered_map<Kmer<k>, std::size_t, Kmer_Hasher<k>> M;
    M.reserve(elem_count);
    t_st = now();
    for(std::size_t i = 0; i < elem_count; ++i)
        M.emplace(kmers[i], vals[i]);
    t_en = now();

    std::cerr << "Inserted " << elem_count << " elements into std::unordered_map in time " << duration(t_en - t_st) << " seconds.\n";


    std::size_t val;

    t_st = now();
    for(std::size_t i = 0; i < elem_count; ++i)
        if(!ht.find(kmers[i], val) || val != vals[i])
        {
            std::cerr << "True positive not found or value incorrect in custom hash-table\n";
            return;
        }
    t_en = now();
    std::cerr << "Searched " << elem_count << " true-positive elements into custom hash-table in time " << duration(t_en - t_st) << " seconds.\n";

    t_st = now();
    for(std::size_t i = 0; i < elem_count; ++i)
        if(M.find(kmers[i]) == M.end())
        {
            std::cerr << "True positive not found in std::unordered_map\n";
            return;
        }
    t_en = now();
    std::cerr << "Searched " << elem_count << " true-positive elements into std::unordered_map in time " << duration(t_en - t_st) << " seconds.\n";

    std::vector<Kmer<k>> kmers_false;
    kmers_false.reserve(elem_count);
    for(std::size_t i = 0; i < elem_count; ++i)
    {
        std::size_t x;
        do
            x = distribution(random_engine);
        while(M.find(x) != M.end());

        kmers_false.emplace_back(Kmer<k>(x));
    }

    t_st = now();
    for(std::size_t i = 0; i < elem_count; ++i)
        if(ht.find(kmers_false[i], val))
        {
            std::cerr << "True negative found in custom hash-table\n";
            return;
        }
    t_en = now();
    std::cerr << "Searched " << elem_count << " true-negative elements into custom hash-table in time " << duration(t_en - t_st) << " seconds.\n";

    t_st = now();
    for(std::size_t i = 0; i < elem_count; ++i)
        if(M.find(kmers_false[i]) != M.end())
        {
            std::cerr << "True negative found in std::unordered_map\n";
            return;
        }
    t_en = now();
    std::cerr << "Searched " << elem_count << " true-negative elements into std::unordered_map in time " << duration(t_en - t_st) << " seconds.\n";

}


template <uint16_t k>
void iterate_subgraphs(const std::string& bin_dir, const std::size_t bin_c)
{
    cuttlefish::Unitig_Write_Distributor lmtigs(bin_dir + "lmtigs", 1024, parlay::num_workers());

    typedef Async_Logger_Wrapper sink_t;
    Output_Sink<sink_t> output_sink;
    clear_file("output.cf3");
    output_sink.init_sink("output.cf3");

    // 100 KB (soft limit) worth of maximal unitig records (FASTA) can be retained in memory per worker, at most, before flushes.
    typedef Character_Buffer<sink_t> op_buf_t;
    std::vector<Padded_Data<op_buf_t>> output_buf(parlay::num_workers(), op_buf_t(output_sink.sink())); // Worker-specific output buffers.

    std::cerr << bin_dir << "; " << bin_c << "\n";
    cuttlefish::Edge_Matrix<k> E(64, "E_");
    std::atomic_uint64_t solved = 0;
    std::atomic_uint64_t v_c = 0;
    std::atomic_uint64_t isolated = 0;
    std::atomic_uint64_t e_c = 0;
    std::atomic_uint64_t max_graph_sz = 0;
    std::atomic_uint64_t label_sz = 0;
    std::atomic_uint64_t disc_edge_c = 0;

    // cuttlefish::Subgraphs_Processor<k> subgraphs(bin_dir, bin_c, E, output_buf);
    parlay::parallel_for(0, bin_c,
        [&](const std::size_t bin_id)
        {
            cuttlefish::Subgraph<k, false> G(bin_dir, bin_id, E, lmtigs, output_buf[parlay::worker_id()].data());
            G.construct();
            v_c += G.size();
            e_c += G.edge_count();

            while(true)
            {
                uint64_t max_sz = max_graph_sz;
                if(max_sz >= G.size())
                    break;

                if(max_graph_sz.compare_exchange_strong(max_sz, G.size()))
                    break;
            }


            G.contract();
            isolated += G.isolated_vertex_count();
            label_sz += G.label_size();
            disc_edge_c += G.discontinuity_edge_count();

            if(++solved % 8 == 0)
                std::cerr << "\rProcessed " << solved << " subgraphs.";
        }
    , 1);
    std::cerr << "\n";

    E.serialize();

    parlay::parallel_for(0, parlay::num_workers(),
                        [&](const std::size_t idx){ output_buf[idx].data().close(); }, 1);

    output_sink.close_sink();

    lmtigs.close();


    std::cerr << "Total vertex count: " << v_c << "\n";
    std::cerr << "Total isolated count: " << isolated << "\n";
    std::cerr << "Total edge count:   " << e_c << "\n";
    std::cerr << "Maximum subgraph-size: " << max_graph_sz << ".\n";
    std::cerr << "Total label size: " << label_sz << "\n";
    std::cerr << "Discontinuity edge count: " << disc_edge_c << ".\n";
}

*/

#include <memory>
#include <functional>
#include "rapidgzip/ParallelGzipReader.hpp"


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

    // static constexpr uint16_t k = 27;
    // static constexpr uint16_t l = 11;

    // const uint8_t f_th = 1;
    // const std::size_t subgraph_count = 2000;
    // const std::size_t parts = 64;
    // const std::size_t unitig_buckets = 1024;
    // const std::string cdbg_path(argv[1]);
    // const std::string output_path(argv[2]);
    // const std::string temp_path(argv[3]);

    // cuttlefish::State_Config::set_edge_threshold(f_th);

    // bootstrap_discontinuity_graph<k>(l, cdbg_path, temp_path, parts, unitig_buckets);

    // cuttlefish::dBG_Contractor<k> dbg_contractor(subgraph_count, parts, unitig_buckets, output_path, temp_path);
    // dbg_contractor.construct();
    // std::cerr << "sizeof(disc-edge<31>): " << sizeof(cuttlefish::Discontinuity_Edge<k>) << "\n";
    // const std::size_t elem_count = std::atoi(argv[1]);
    // const double lf = 0.75;
    // benchmark_hash_table<k>(elem_count, lf);

    // test_unitig_file(argv[1], argv[2]);

    // const uint64_t br = compute_breakpoints<k, l>(argv[1]);
    // std::cout << "#breakpoints: " << br << "\n";
    // std::cout << (validate_canonical_minimizer<k, l>(argv[1]) ? "Canonical minimizers correctly computed.\n" : "Canonical minimizers found wrong.\n");
    // minimizer_iterator_test<k, l>(argv[1]);
    // std::cout << (Index_Validator<k, l>::validate(argv[1], argv[2]) ? "Index cross-checking successful.\n" : "Index is incorrect.\n");


    // const std::string bin_dir(argv[1]);
    // const std::size_t bin_c(std::atoi(argv[2]));
    // iterate_subgraphs<k>(bin_dir, bin_c);
    // cuttlefish::Parser(argv[1], std::atoi(argv[2])).parse();
/*
	auto file_reader = std::make_unique<StandardFileReader>(argv[1]);
	const size_t thread_count = std::atoi(argv[2]);

	const size_t chunk_size = 8 * 1024 * 1024;
	auto chunk = new char[chunk_size];
    std::size_t total_bytes_read = 0;

	rapidgzip::ParallelGzipReader<> reader(std::move(file_reader), thread_count);
	reader.setCRC32Enabled(false);

	// std::function<void( const std::shared_ptr<rapidgzip::ChunkData>, size_t, size_t )> writeFunctor;
	// auto totalBytesRead = reader.read( writeFunctor );


	while (true)
    {
		const auto read = reader.read(chunk, chunk_size);
		if (!read)
			break;

        total_bytes_read += read;
	}

    delete[] chunk;

	std::cerr << "Total bytes decompressed: " << total_bytes_read << ".\n";
*/
    return 0;
}


/*
/usr/bin/time rapidgzip -dk -P 1 /ssd2/jamshed/seq-reads/SRR3440495__NIST_sequencing_of_the_AJ_Trio-_2x250bp_overlapping_libraries_with_nominal_350bp_insert_size_designed_for_DISCOVAR_assembly_NIST_HG004_NA24143-D3_S3_L002_004__2.fastq.gz -o /dev/null
40.67user 2.65system 0:43.34elapsed 99%CPU (0avgtext+0avgdata 41472maxresident)k
0inputs+0outputs (0major+119600minor)pagefaults 0swaps

/usr/bin/time ./cuttlefish-private/build/src/test /ssd2/jamshed/seq-reads/SRR3440495__NIST_sequencing_of_the_AJ_Trio-_2x250bp_overlapping_libraries_with_nominal_350bp_insert_size_designed_for_DISCOVAR_assembly_NIST_HG004_NA24143-D3_S3_L002_004__2.fastq.gz 1 1 > /dev/null
eof?: 1
193.61user 4.74system 3:18.40elapsed 99%CPU (0avgtext+0avgdata 678684maxresident)k
0inputs+0outputs (0major+814642minor)pagefaults 0swaps
*/
