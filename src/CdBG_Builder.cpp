
#include "CdBG_Builder.hpp"
#include "kseq/kseq.h"

#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <ctime>
#include <chrono>
#include "zlib.h"


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(int, read);


CdBG_Builder::CdBG_Builder(const std::string& ref_file, const uint16_t k):
        ref_file(ref_file), k(k)
{
    Kmer::set_k(k);
}


void CdBG_Builder::construct(const std::string& output_file)
{
    std::cout << "Classifying the vertices.\n";
    classify_vertices();

    std::cout << "Classification of the vertices complete. Now outputting the maximal unitigs.\n";
    output_maximal_unitigs(output_file);
}


void CdBG_Builder::classify_vertices()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA file containing the reference.
    FILE* input = fopen(ref_file.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Initialize the parser.
    kseq_t* parser = kseq_init(fileno(input));


    // Mark all the k-mers as unvisited.
    Vertices.clear();


    // Parse sequences one-by-one.
    uint32_t seqCount = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        process_first_kmer(seq);

        for(uint32_t kmer_idx = 1; kmer_idx < seq_len - k; ++kmer_idx)
            process_internal_kmer(seq, kmer_idx);

        process_last_kmer(seq, seq_len);
    }


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done classifying the vertices. Time taken = " << elapsed_seconds << " seconds.\n";
}


bool CdBG_Builder::is_self_loop(const char* seq, const cuttlefish::kmer_t& kmer, const uint32_t kmer_idx) const
{
    return kmer.is_same_kmer(cuttlefish::kmer_t(seq, kmer_idx + 1));
}


void CdBG_Builder::process_first_kmer(const char* seq)
{
    const cuttlefish::kmer_t kmer(seq, 0);
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t next_nucl = seq[k];


    // The k-mer is already classified as a complex node.
    if(Vertices.is_present(kmer_hat) && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(seq, kmer, 0))
    {
        Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_MULTI_OUT));
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!Vertices.is_present(kmer_hat))
            Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, next_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = Vertices[kmer_hat].decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.exit == next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != next_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!Vertices.is_present(kmer_hat))
            Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, complement(next_nucl)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = Vertices[kmer_hat].decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == complement(next_nucl))
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != complement(next_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
        }
        
    }
}


void CdBG_Builder::process_last_kmer(const char* seq, const uint32_t seq_len)
{
    const cuttlefish::kmer_t kmer(seq, seq_len - k);
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t prev_nucl = seq[seq_len - k - 1];


    // The k-mer is already classified as a complex node.
    if(Vertices.is_present(kmer_hat) && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!Vertices.is_present(kmer_hat))
            Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, prev_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = Vertices[kmer_hat].decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == prev_nucl)
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
                
                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != prev_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!Vertices.is_present(kmer_hat))
            Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, complement(prev_nucl)));
        else    // The sentinel k-mer had been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = Vertices[kmer_hat].decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.exit == complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
        }
    }
}


void CdBG_Builder::process_internal_kmer(const char* seq, const uint32_t kmer_idx)
{
    const cuttlefish::kmer_t kmer(seq, kmer_idx);
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t prev_nucl = seq[kmer_idx - 1];
    const cuttlefish::nucleotide_t next_nucl = seq[kmer_idx + k];


    // The k-mer is already classified as a complex node.
    if(Vertices.is_present(kmer_hat) && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(seq, kmer, kmer_idx))
    { 
        Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_MULTI_OUT));
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        // The k-mer is encountered for the first time, and in the forward direction.
        if(!Vertices.is_present(kmer_hat))
            Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, prev_nucl, next_nucl));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = Vertices[kmer_hat].decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == prev_nucl && vertex.exit == next_nucl)
                    return;

                if(vertex.enter != prev_nucl && vertex.exit != next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
                else if(vertex.enter != prev_nucl)
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else    // vertex.exit != next_nucl
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != next_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != prev_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
        }
    }
    else
    {
        // The k-mer is encountered for the first time, and in the backward direction.
        if(!Vertices.is_present(kmer_hat))
            Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, complement(next_nucl), complement(prev_nucl)));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = Vertices[kmer_hat].decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == complement(next_nucl) && vertex.exit == complement(prev_nucl))
                    return;

                if(vertex.enter != complement(next_nucl) && vertex.exit != complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
                else if(vertex.enter != complement(next_nucl))
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else if(vertex.exit != complement(prev_nucl))
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;

                Vertices[kmer_hat] = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    Vertices[kmer_hat] = Vertex_Encoding(vertex);
                }
            }
        }
    }
}


void CdBG_Builder::output_maximal_unitigs(const std::string& output_file)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA file containing the reference.
    FILE* input = fopen(ref_file.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Initialize the parser.
    kseq_t* parser = kseq_init(fileno(input));


    // Open the output file.
    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if(!output)
    {
        std::cerr << "Error opening output file " << output_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Parse sequences one-by-one.
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;



        uint32_t unipath_start_idx, unipath_end_idx;

        for(uint32_t kmer_idx = 0; kmer_idx <= seq_len - k; ++kmer_idx)
        {
            if(is_unipath_start(seq, kmer_idx))
                unipath_start_idx = kmer_idx;

            if(is_unipath_end(seq, kmer_idx))
            {
                unipath_end_idx = kmer_idx;
                output_unitig(seq, unipath_start_idx, unipath_end_idx, output);
            }
        }
    }

    output.close();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done classifying the vertices. Time taken = " << elapsed_seconds << " seconds.\n";
}


bool CdBG_Builder::is_unipath_start(const char* seq, const uint32_t kmer_idx) const
{
    const cuttlefish::kmer_t kmer(seq, kmer_idx);
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::state_t state = Vertices[kmer_hat].state();

    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::FWD && state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    if(dir == cuttlefish::BWD && state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;


    // assert(kmer_idx > 0);

    
    cuttlefish::kmer_t prev_kmer(seq, kmer_idx - 1);
    cuttlefish::kmer_t prev_kmer_hat = prev_kmer.canonical();
    cuttlefish::kmer_dir_t prev_kmer_dir = prev_kmer.direction(prev_kmer_hat);
    cuttlefish::state_t prev_state = Vertices[prev_kmer_hat].state();

    if(prev_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD && prev_state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::BWD && prev_state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    
    return false;
}


bool CdBG_Builder::is_unipath_end(const char* seq, const uint32_t kmer_idx) const
{
    cuttlefish::kmer_t kmer(seq, kmer_idx);
    cuttlefish::kmer_t kmer_hat = kmer.canonical();
    cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    cuttlefish::state_t state = Vertices[kmer_hat].state();

    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::FWD && state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::BWD && state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;


    // assert(kmer_idx < ref.length() - k);

    cuttlefish::kmer_t next_kmer(seq, kmer_idx + 1);
    cuttlefish::kmer_t next_kmer_hat = next_kmer.canonical();
    cuttlefish::kmer_dir_t next_kmer_dir = next_kmer.direction(next_kmer_hat);
    cuttlefish::state_t next_state = Vertices[next_kmer_hat].state();

    if(next_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::FWD && next_state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::BWD && next_state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;


    return false;
}


void CdBG_Builder::output_unitig(const char* seq, const uint32_t start_idx, const uint32_t end_idx, std::ofstream& output)
{
    cuttlefish::kmer_t u(seq, start_idx);
    cuttlefish::kmer_t v(seq, end_idx);
    cuttlefish::kmer_t u_hat = u.canonical();
    cuttlefish::kmer_t v_hat = v.canonical();


    if(Vertices[u_hat].is_outputted())   // or Vertices[v_hat].is_outputted()
        return;

    
    Vertices[u_hat] = Vertices[u_hat].outputted();
    Vertices[v_hat] = Vertices[v_hat].outputted();

    const uint32_t unipath_len = end_idx - start_idx + 1 + k - 1;
    std::vector<char> unipath;
    unipath.reserve(unipath_len);


    if(u < v.reverse_complement())
        for(int32_t idx = start_idx; idx <= (int32_t)end_idx + k - 1; ++idx)
            unipath.push_back(seq[idx]);
    else
        for(int32_t idx = end_idx + k - 1; idx >= (int32_t)start_idx; --idx)
            unipath.push_back(complement(seq[idx]));

    
    output << std::string(unipath.begin(), unipath.end()) << "\n";
}


void CdBG_Builder::print_vertices() const
{
    Vertices.print_hash_table();
    // for(auto vertex: Vertices)
    //     std::cout << vertex.first << " : " << vertex.second.decode() << "\n";
}
