
#include "globals.hpp"
#include "CdBG.hpp"

#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>


bool CdBG::is_self_loop(const std::string &ref, const Kmer_Str& kmer, const uint32_t kmer_idx) const
{
    return kmer.is_same_kmer(Kmer_Str(ref.substr(kmer_idx + 1, k)));
}


void CdBG::process_first_kmer(const std::string& ref)
{
    const Kmer_Str kmer(ref.substr(0, k));
    const Kmer_Str kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t next_nucl = ref[k];


    if(Vertices.find(kmer_hat) != Vertices.end() && Vertices[kmer_hat].state == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    if(is_self_loop(ref, kmer, 0))
    {
        Vertices[kmer_hat] = Vertex(cuttlefish::MULTI_IN_MULTI_OUT);
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        if(Vertices.find(kmer_hat) == Vertices.end())
            Vertices[kmer_hat] = Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, next_nucl);
        else
        {
            Vertex& vertex = Vertices[kmer_hat];

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.exit == next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
        }
    }
    else
    {
        if(Vertices.find(kmer_hat) == Vertices.end())
            Vertices[kmer_hat] = Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, complement(next_nucl));
        else
        {
            Vertex& vertex = Vertices[kmer_hat];

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == complement(next_nucl))
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != complement(next_nucl))
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
        }
    }
}


void CdBG::process_last_kmer(const std::string& ref)
{
    const Kmer_Str kmer(ref.substr(ref.length() - k, k));
    const Kmer_Str kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    const cuttlefish::nucleotide_t prev_nucl = ref[ref.length() - k - 1];


    if(Vertices.find(kmer_hat) != Vertices.end() && Vertices[kmer_hat].state == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    
    if(dir == cuttlefish::FWD)
    {
        if(Vertices.find(kmer_hat) == Vertices.end())
            Vertices[kmer_hat] = Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, prev_nucl);
        else
        {
            Vertex& vertex = Vertices[kmer_hat];

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter == prev_nucl)
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != prev_nucl)
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
        }
    }
    else
    {
        if(Vertices.find(kmer_hat) == Vertices.end())
            Vertices[kmer_hat] = Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, complement(prev_nucl));
        else
        {
            Vertex& vertex = Vertices[kmer_hat];

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.exit == complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
                vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
        }
    }
}


void CdBG::process_internal_kmer(const std::string& ref, const uint32_t kmer_idx)
{
    const Kmer_Str kmer(ref.substr(kmer_idx, k));
    const Kmer_Str kmer_hat = kmer.canonical();
    cuttlefish::kmer_dir_t dir = kmer.direction(kmer_hat);
    cuttlefish::nucleotide_t prev_nucl = ref[kmer_idx - 1];
    cuttlefish::nucleotide_t next_nucl = ref[kmer_idx + k];


    if(Vertices.find(kmer_hat) != Vertices.end() && Vertices[kmer_hat].state == cuttlefish::MULTI_IN_MULTI_OUT)
        return;

    if(is_self_loop(ref, kmer, kmer_idx))
    {
        Vertices[kmer_hat] = Vertex(cuttlefish::MULTI_IN_MULTI_OUT);
        return;
    }


    if(dir == cuttlefish::FWD)
    {
        if(Vertices.find(kmer_hat) == Vertices.end())
            Vertices[kmer_hat] = Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, prev_nucl, next_nucl);
        else
        {
            Vertex& vertex = Vertices[kmer_hat];

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter != prev_nucl && vertex.exit != next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
                else if(vertex.enter != prev_nucl)
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else if(vertex.exit != next_nucl)
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != next_nucl)
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != prev_nucl)
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
        }
    }
    else
    {
        if(Vertices.find(kmer_hat) == Vertices.end())
            Vertices[kmer_hat] = Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, complement(next_nucl), complement(prev_nucl));
        else
        {
            Vertex& vertex = Vertices[kmer_hat];

            if(vertex.state == cuttlefish::SINGLE_IN_SINGLE_OUT)
            {
                if(vertex.enter != complement(next_nucl) && vertex.exit != complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
                else if(vertex.enter != complement(next_nucl))
                    vertex.state = cuttlefish::MULTI_IN_SINGLE_OUT;
                else if(vertex.exit != complement(prev_nucl))
                    vertex.state = cuttlefish::SINGLE_IN_MULTI_OUT;
            }
            else if(vertex.state == cuttlefish::MULTI_IN_SINGLE_OUT)
            {
                if(vertex.exit != complement(prev_nucl))
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
            else    // vertex.state == cuttlefish::SINGLE_IN_MULTI_OUT
            {
                if(vertex.enter != complement(next_nucl))
                    vertex.state = cuttlefish::MULTI_IN_MULTI_OUT;
            }
        }
    }
}


void CdBG::classify_vertices()
{
    std::ifstream references(ref_file.c_str(), std::ifstream::in);
    if(!references)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Mark all the k-mers as unvisited.
    Vertices.clear();


    std::string ref;

    while(references >> ref)
    {
        process_first_kmer(ref);

        for(uint32_t kmer_idx = 1; kmer_idx < ref.length() - k; ++kmer_idx)
            process_internal_kmer(ref, kmer_idx);

        process_last_kmer(ref);
    }

    
    references.close();   
}


void CdBG::construct(const std::string& output_file)
{
    classify_vertices();

    output_maximal_unitigs(output_file);
}


bool CdBG::is_unipath_start(const std::string& ref, const uint32_t kmer_idx) const
{
    const Kmer_Str kmer(ref.substr(kmer_idx, k));
    const Kmer_Str kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t kmer_dir = kmer.direction(kmer_hat);
    const cuttlefish::state_t state = (Vertices.find(kmer_hat) -> second).state;

    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(kmer_dir == cuttlefish::FWD && state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    if(kmer_dir == cuttlefish::BWD && state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;


    assert(kmer_idx > 0);

    const Kmer_Str prev_kmer(ref.substr(kmer_idx - 1, k));
    const Kmer_Str prev_kmer_hat = prev_kmer.canonical();
    const cuttlefish::kmer_dir_t prev_kmer_dir = prev_kmer.direction(prev_kmer_hat);
    const cuttlefish::state_t prev_state = (Vertices.find(prev_kmer_hat) -> second).state;

    if(prev_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::FWD && prev_state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    if(prev_kmer_dir == cuttlefish::BWD && prev_state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    
    return false;
}


bool CdBG::is_unipath_end(const std::string& ref, const uint32_t kmer_idx) const
{
    const Kmer_Str kmer(ref.substr(kmer_idx, k));
    const Kmer_Str kmer_hat = kmer.canonical();
    const cuttlefish::kmer_dir_t kmer_dir = kmer.direction(kmer_hat);
    const cuttlefish::state_t state = (Vertices.find(kmer_hat) -> second).state;

    if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(kmer_dir == cuttlefish::FWD && state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    if(kmer_dir == cuttlefish::BWD && state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;


    assert(kmer_idx < ref.length() - k);

    const Kmer_Str next_kmer(ref.substr(kmer_idx + 1, k));
    const Kmer_Str next_kmer_hat = next_kmer.canonical();
    const cuttlefish::kmer_dir_t next_kmer_dir = next_kmer.direction(next_kmer_hat);
    const cuttlefish::state_t next_state = (Vertices.find(next_kmer_hat) -> second).state;

    if(next_state == cuttlefish::MULTI_IN_MULTI_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::FWD && next_state == cuttlefish::MULTI_IN_SINGLE_OUT)
        return true;

    if(next_kmer_dir == cuttlefish::BWD && next_state == cuttlefish::SINGLE_IN_MULTI_OUT)
        return true;

    
    return false;
}


void CdBG::output_maximal_unitigs(const std::string& output_file)
{
    std::ifstream references(ref_file.c_str(), std::ifstream::in);
    if(!references)
    {
        std::cerr << "Error opening input file " << ref_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if(!output)
    {
        std::cerr << "Error opening output file " << output_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::string ref;
    while(references >> ref)
    {
        uint32_t unipath_start_idx, unipath_end_idx;

        for(uint32_t kmer_idx = 0; kmer_idx <= ref.length() - k; ++kmer_idx)
        {
            if(is_unipath_start(ref, kmer_idx))
                unipath_start_idx = kmer_idx;

            if(is_unipath_end(ref, kmer_idx))
            {
                unipath_end_idx = kmer_idx;
                output_unipath(ref, output, unipath_start_idx, unipath_end_idx);
            }
        }
    }


    references.close();
    output.close();
}


void CdBG::output_unipath(const std::string& ref, std::ofstream &output, const uint32_t start_idx, const uint32_t end_idx)
{
    const Kmer_Str u(ref.substr(start_idx, k));
    const Kmer_Str v(ref.substr(end_idx, k));
    const Kmer_Str u_hat = u.canonical();
    const Kmer_Str v_hat = v.canonical();

    if(Vertices[u_hat].outputted)   // or Vertices[v_hat].outputted
        return;

    
    Vertices[u_hat].outputted = true;
    Vertices[v_hat].outputted = true;

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


void CdBG::print_vertices() const
{
    for(auto p = Vertices.begin(); p != Vertices.end(); ++p)
        std::cout << p -> first << " : " << p -> second << "\n";
}
