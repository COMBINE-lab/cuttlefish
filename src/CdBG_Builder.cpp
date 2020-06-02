
#include "CdBG_Builder.hpp"
#include "Directed_Kmer.hpp"
#include "Annotated_Kmer.hpp"
#include "kseq/kseq.h"

#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>
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


void CdBG_Builder::construct(const std::string& kmers_file, const uint64_t kmer_count, const std::string& output_file)
{
    std::cout << "Constructing the minimal perfect hash function.\n";
    Vertices.construct(kmers_file, kmer_count);

    std::cout << "Classifying the vertices.\n";
    classify_vertices();

    std::cout << "Outputting the maximal unitigs.\n";
    output_maximal_unitigs(output_file);

    Vertices.clear();
}


void CdBG_Builder::classify_vertices()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
    FILE* input = fopen(ref_file.c_str(), "r");
    if(input == NULL)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    
    // Initialize the parser.
    kseq_t* parser = kseq_init(fileno(input));

    // Parse sequences one-by-one, and continue partial classification of the k-mers through them.
    uint32_t seqCount = 0;
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;

        std::cout << "Processing sequence " << ++seqCount << ", with length " << seq_len << ".\n";

        Directed_Kmer curr_kmer(cuttlefish::kmer_t(seq, 0));
        Directed_Kmer next_kmer(cuttlefish::kmer_t(seq, 1));

        process_first_kmer(curr_kmer.canonical, curr_kmer.dir, next_kmer.canonical, seq[k]);

        for(uint32_t kmer_idx = 1; kmer_idx < seq_len - k; ++kmer_idx)
        {
            curr_kmer = next_kmer;
            next_kmer.roll_to_next_kmer(seq[kmer_idx + k]);

            process_internal_kmer(curr_kmer.canonical, curr_kmer.dir, next_kmer.canonical, seq[kmer_idx - 1], seq[kmer_idx + k]);
        }

        process_last_kmer(next_kmer.canonical, next_kmer.dir, seq[seq_len - k - 1]);
    }


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done classifying the vertices. Time taken = " << elapsed_seconds << " seconds.\n";
}


bool CdBG_Builder::is_self_loop(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_t& next_kmer_hat) const
{
    return kmer_hat == next_kmer_hat;
}


void CdBG_Builder::process_first_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t next_nucl)
{
    Vertex_Encoding& vertex_encoding = Vertices[kmer_hat];


    // The k-mer is already classified as a complex node.
    // if(Vertices.is_visited(kmer_hat) && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
    if(vertex_encoding.is_visited() && vertex_encoding.state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
    {
        vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_MULTI_OUT));
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, next_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.exit == next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != next_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding = Vertex_Encoding(vertex);
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, complement(next_nucl)));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == complement(next_nucl))
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != complement(next_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
        
    }
}


void CdBG_Builder::process_last_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t dir, const cuttlefish::nucleotide_t prev_nucl)
{
    Vertex_Encoding& vertex_encoding = Vertices[kmer_hat];


    // The k-mer is already classified as a complex node.
    if(vertex_encoding.is_visited() && vertex_encoding.state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, prev_nucl));
        else    // The sentinel k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == prev_nucl)
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
                
                vertex_encoding = Vertex_Encoding(vertex);
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != prev_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }
    else
    {
        // The sentinel k-mer is encountered for the first time, and in the backward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, complement(prev_nucl)));
        else    // The sentinel k-mer had been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.exit == complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding= Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                vertex_encoding = Vertex_Encoding(vertex);
            }
        }
    }
}


void CdBG_Builder::process_internal_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_dir_t& dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t prev_nucl, const cuttlefish::nucleotide_t next_nucl)
{
    Vertex_Encoding& vertex_encoding = Vertices[kmer_hat];


    // The k-mer is already classified as a complex node.
    if(vertex_encoding.is_visited() && vertex_encoding.state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(kmer_hat, next_kmer_hat))
    { 
        vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_MULTI_OUT));
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        // The k-mer is encountered for the first time, and in the forward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, prev_nucl, next_nucl));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

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

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != next_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != prev_nucl)
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }
    else
    {
        // The k-mer is encountered for the first time, and in the backward direction.
        if(!vertex_encoding.is_visited())
            vertex_encoding = Vertex_Encoding(Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, complement(next_nucl), complement(prev_nucl)));
        else    // The k-mer has been visited earlier and has some state; modify it accordingly.
        {
            Vertex vertex = vertex_encoding.decode();

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

                vertex_encoding = Vertex_Encoding(vertex);
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.exit != complement(prev_nucl))
                {
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;

                    vertex_encoding = Vertex_Encoding(vertex);
                }
            }
        }
    }
}


void CdBG_Builder::output_maximal_unitigs(const std::string& output_file)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the file handler for the FASTA / FASTQ file containing the reference.
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


    // Parse sequences one-by-one, and output each unique maximal unitig encountered through them.
    while(kseq_read(parser) >= 0)
    {
        const char* seq = parser->seq.s;
        const size_t seq_len = parser->seq.l;


        Annotated_Kmer prev_annot_kmer;
        Annotated_Kmer curr_annot_kmer(cuttlefish::kmer_t(seq, 0), 0, Vertices);
        Annotated_Kmer next_annot_kmer(cuttlefish::kmer_t(seq, 1), 1, Vertices);

        Annotated_Kmer unipath_start_kmer;


        // Process the first k-mer.

        // A maximal unitig starts at the beginning of a sequence.
        unipath_start_kmer = curr_annot_kmer;
        if(is_unipath_end(curr_annot_kmer.state, curr_annot_kmer.dir, next_annot_kmer.state, next_annot_kmer.dir))
            output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);


        // Process the internal k-mers.

        for(uint32_t kmer_idx = 1; kmer_idx < seq_len - k; ++kmer_idx)
        {
            prev_annot_kmer = curr_annot_kmer;
            curr_annot_kmer = next_annot_kmer;
            next_annot_kmer.roll_to_next_kmer(seq[kmer_idx + k], Vertices);


            if(is_unipath_start(curr_annot_kmer.state, curr_annot_kmer.dir, prev_annot_kmer.state, prev_annot_kmer.dir))
                unipath_start_kmer = curr_annot_kmer;

            if(is_unipath_end(curr_annot_kmer.state, curr_annot_kmer.dir, next_annot_kmer.state, next_annot_kmer.dir))
                output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);
        }

        prev_annot_kmer = curr_annot_kmer;
        curr_annot_kmer = next_annot_kmer;


        // Process the last k-mer.
        
        if(is_unipath_start(curr_annot_kmer.state, curr_annot_kmer.dir, prev_annot_kmer.state, prev_annot_kmer.dir))
            unipath_start_kmer = curr_annot_kmer;

        // A maximal unitig ends at the ending of a sequence.
        output_unitig(seq, unipath_start_kmer, curr_annot_kmer, output);
    }

    
    // Close the output file.
    if(output.fail())
        std::cerr << "Errors had been encountered for the output stream.\n";

    output.close();


    // Close the parser and the input file.
    kseq_destroy(parser);
    fclose(input);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done outputting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


bool CdBG_Builder::is_unipath_start(const cuttlefish::state_t state, const cuttlefish::kmer_dir_t dir, const cuttlefish::state_t prev_kmer_state, const cuttlefish::kmer_dir_t prev_kmer_dir) const
{
    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(state == cuttlefish::MULTI_IN_SINGLE_OUT)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(state == cuttlefish::SINGLE_IN_MULTI_OUT)
            return true;


    // assert(kmer_idx > 0);


    if(prev_kmer_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD)
    {
        if(prev_kmer_state == cuttlefish::SINGLE_IN_MULTI_OUT)
            return true;
    }
    else    // prev_kmer_dir == cuttlefish::BWD
        if(prev_kmer_state == cuttlefish::MULTI_IN_SINGLE_OUT)
            return true;

    
    return false;
}


bool CdBG_Builder::is_unipath_end(const cuttlefish::state_t state, const cuttlefish::kmer_dir_t dir, const cuttlefish::state_t next_kmer_state, const cuttlefish::kmer_dir_t next_kmer_dir) const
{
    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::FWD)
    {
        if(state == cuttlefish::SINGLE_IN_MULTI_OUT)
            return true;
    }
    else    // dir == cuttlefish::BWD
        if(state == cuttlefish::MULTI_IN_SINGLE_OUT)
            return true;


    // assert(kmer_idx < ref.length() - k);


    if(next_kmer_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::FWD)
    {
        if(next_kmer_state == cuttlefish::MULTI_IN_SINGLE_OUT)
            return true;
    }
    else    // next_kmer_dir == cuttlefish::BWD
        if(next_kmer_state == cuttlefish::SINGLE_IN_MULTI_OUT)
            return true;


    return false;
}


void CdBG_Builder::output_unitig(const char* seq, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, std::ofstream& output)
{
    Vertex_Encoding& vertex_encoding_s = Vertices[start_kmer.canonical];
    Vertex_Encoding& vertex_encoding_e = Vertices[end_kmer.canonical];


    // This is to avoid race conditions that may arise while multi-threading.
    // If two threads try to output the same unitig at the same time but
    // encounter it in the opposite orientations, then data races may arise.
    Vertex_Encoding& vertex_encoding = (start_kmer.canonical < end_kmer.canonical ?
                                        vertex_encoding_s : vertex_encoding_e);
    if(vertex_encoding.is_outputted())
        return;
    

    vertex_encoding_s = vertex_encoding_s.outputted();
    vertex_encoding_e = vertex_encoding_e.outputted();

    write_path(seq, start_kmer.idx, end_kmer.idx, start_kmer.kmer < end_kmer.rev_compl, output);
}


void CdBG_Builder::write_path(const char* seq, const uint32_t start_kmer_idx, const uint32_t end_kmer_idx, const bool in_forward, std::ofstream& output) const
{
    if(in_forward)
        for(uint32_t idx = start_kmer_idx; idx <= end_kmer_idx + k - 1; ++idx)
            output << seq[idx];
    else
    {
        // To avoid underflow of unsigned integers, the flanking indices are incremented by 1.
        uint32_t idx = end_kmer_idx + k;
        while(idx > start_kmer_idx)
        {
            idx--;
            output << complement(seq[idx]);
        }
    }

    output << "\n";
}
