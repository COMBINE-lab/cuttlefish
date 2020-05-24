
#ifndef GLOBALS_HPP
#define GLOBALS_HPP


#include <string>
#include <iostream>


class Kmer_Str;
class Kmer;

namespace cuttlefish
{
    // typedef Kmer_Str kmer_t;
    typedef Kmer kmer_t;
    typedef bool kmer_dir_t;
    typedef char nucleotide_t;
    typedef uint8_t state_t;


    constexpr nucleotide_t PLACEHOLDER_NUCLEOTIDE = 'N';

    constexpr kmer_dir_t FWD = true;
    constexpr kmer_dir_t BWD = false;

    constexpr uint8_t SINGLE_IN_SINGLE_OUT = 0;
    constexpr uint8_t MULTI_IN_SINGLE_OUT = 1;
    constexpr uint8_t SINGLE_IN_MULTI_OUT = 2;
    constexpr uint8_t MULTI_IN_MULTI_OUT = 3;
}


inline cuttlefish::nucleotide_t complement(const cuttlefish::nucleotide_t nucleotide)
{
    switch (nucleotide)
    {
    case 'A':
        return 'T';

    case 'C':
        return 'G';

    case 'G':
        return 'C';

    case 'T':
        return 'A';
    
    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        return 'A';
        
        // std::cerr << "Invalid nucleotide " << nucleotide << " encountered. Aborting.";
        // std::exit(EXIT_FAILURE);
    }
}


#endif
