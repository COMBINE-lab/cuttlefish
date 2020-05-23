
#include "CdBG_Builder.hpp"

#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>


CdBG_Builder::CdBG_Builder(const std::string& ref_file, const uint16_t k):
        ref_file(ref_file), k(k)
{
    Kmer::set_k(k);
}


void CdBG_Builder::construct(const std::string& output_file)
{
    classify_vertices();

    output_maximal_unitigs(output_file);
}


void CdBG_Builder::classify_vertices()
{
    // Open the file containing newline-separated referenes.
    std::ifstream refs(ref_file.c_str(), std::ifstream::in);
    if(!refs)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Mark all the k-mers as unvisited.
    Vertices.clear();


    std::string ref;

    while(refs >> ref)
    {
        process_first_kmer(ref);

        for(uint32_t kmer_idx = 1; kmer_idx < ref.length() - k; ++kmer_idx)
            process_internal_kmer(ref, kmer_idx);

        process_last_kmer(ref);
    }


    refs.close();
}


bool CdBG_Builder::is_self_loop(const std::string& ref, const cuttlefish::kmer_t& kmer, const uint32_t kmer_idx) const
{
    return kmer.is_same_kmer(cuttlefish::kmer_t(ref.substr(kmer_idx + 1, k)));
}


void CdBG_Builder::process_first_kmer(const std::string& ref)
{
    const cuttlefish::kmer_t kmer(ref.substr(0, k));
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t next_nucl = ref[k];


    // The k-mer is already classified as a complex node.
    if(Vertices.find(kmer_hat) != Vertices.end() && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(ref, kmer, 0))
    {
        Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_MULTI_OUT));
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(Vertices.find(kmer_hat) == Vertices.end())
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
        if(Vertices.find(kmer_hat) == Vertices.end())
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


void CdBG_Builder::process_last_kmer(const std::string& ref)
{
    const cuttlefish::kmer_t kmer(ref.substr(ref.length() - k, k));
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t prev_nucl = ref[ref.length() - k - 1];


    // The k-mer is already classified as a complex node.
    if(Vertices.find(kmer_hat) != Vertices.end() && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;


    if(dir == cuttlefish::FWD)
    {
        // The sentinel k-mer is encountered for the first time, and in the forward direction.
        if(Vertices.find(kmer_hat) == Vertices.end())
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
        if(Vertices.find(kmer_hat) == Vertices.end())
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


void CdBG_Builder::process_internal_kmer(const std::string& ref, const uint32_t kmer_idx)
{
    const cuttlefish::kmer_t kmer(ref.substr(kmer_idx, k));
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t prev_nucl = ref[kmer_idx - 1];
    const cuttlefish::nucleotide_t next_nucl = ref[kmer_idx + k];


    // The k-mer is already classified as a complex node.
    if(Vertices.find(kmer_hat) != Vertices.end() && Vertices[kmer_hat].state() == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    // The k-mer forms a self-loop with the next k-mer.
    if(is_self_loop(ref, kmer, kmer_idx))
    { 
        Vertices[kmer_hat] = Vertex_Encoding(Vertex(cuttlefish::MULTI_IN_MULTI_OUT));
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        // The k-mer is encountered for the first time, and in the forward direction.
        if(Vertices.find(kmer_hat) == Vertices.end())
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
        if(Vertices.find(kmer_hat) == Vertices.end())
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
    // Open the input file containing newline-separated references.
    std::ifstream refs(ref_file.c_str(), std::ifstream::in);
    if(!refs)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Open the output file.
    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if(!output)
    {
        std::cerr << "Error opening output file " << output_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::string ref;
    while(refs >> ref)
    {
        uint32_t unipath_start_idx, unipath_end_idx;

        for(uint32_t kmer_idx = 0; kmer_idx <= ref.length() - k; ++kmer_idx)
        {
            if(is_unipath_start(ref, kmer_idx))
                unipath_start_idx = kmer_idx;

            if(is_unipath_end(ref, kmer_idx))
            {
                unipath_end_idx = kmer_idx;
                output_unitig(ref, unipath_start_idx, unipath_end_idx, output);
            }
        }
    }


    output.close();
}


bool CdBG_Builder::is_unipath_start(const std::string& ref, const uint32_t kmer_idx) const
{
    const cuttlefish::kmer_t kmer(ref.substr(kmer_idx, k));
    const cuttlefish::kmer_t kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::state_t state = (Vertices.find(kmer_hat) -> second).state();

    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::FWD && state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    if(dir == cuttlefish::BWD && state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;


    assert(kmer_idx > 0);

    
    cuttlefish::kmer_t prev_kmer(ref.substr(kmer_idx - 1, k));
    cuttlefish::kmer_t prev_kmer_hat = prev_kmer.canonical();
    cuttlefish::kmer_dir_t prev_kmer_dir = prev_kmer.direction(prev_kmer_hat);
    cuttlefish::state_t prev_state = (Vertices.find(prev_kmer_hat) -> second).state();

    if(prev_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD && prev_state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::BWD && prev_state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    
    return false;
}


bool CdBG_Builder::is_unipath_end(const std::string& ref, const uint32_t kmer_idx) const
{
    cuttlefish::kmer_t kmer(ref.substr(kmer_idx, k));
    cuttlefish::kmer_t kmer_hat = kmer.canonical();
    cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    cuttlefish::state_t state = (Vertices.find(kmer_hat) -> second).state();

    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::FWD && state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    if(dir == cuttlefish::BWD && state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;


    assert(kmer_idx < ref.length() - k);

    cuttlefish::kmer_t next_kmer(ref.substr(kmer_idx + 1, k));
    cuttlefish::kmer_t next_kmer_hat = next_kmer.canonical();
    cuttlefish::kmer_dir_t next_kmer_dir = next_kmer.direction(next_kmer_hat);
    cuttlefish::state_t next_state = (Vertices.find(next_kmer_hat) -> second).state();

    if(next_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::FWD && next_state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::BWD && next_state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;


    return false;
}


void CdBG_Builder::output_unitig(const std::string& ref, const uint32_t start_idx, const uint32_t end_idx, std::ofstream& output)
{
    cuttlefish::kmer_t u(ref.substr(start_idx, k));
    cuttlefish::kmer_t v(ref.substr(end_idx, k));
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
            unipath.push_back(ref[idx]);
    else
        for(int32_t idx = end_idx + k - 1; idx >= (int32_t)start_idx; --idx)
            unipath.push_back(complement(ref[idx]));

    
    output << std::string(unipath.begin(), unipath.end()) << "\n";
}


void CdBG_Builder::print_vertices() const
{
    for(auto vertex: Vertices)
        std::cout << vertex.first << " : " << vertex.second.decode() << "\n";
}
