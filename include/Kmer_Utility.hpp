
#ifndef KMER_UTILITY_HPP
#define KMER_UTILITY_HPP



#include <cstdint>
#include <cstdlib>
#include <iostream>


class Kmer_Utility
{
protected:

    // A = 0, C = 1, G = 2, T = 3.
    // Note that, this is not possible to change this mapping w/o modifications to the
    // interfacing of our code with the KMC api. This mapping is essential for some
    // crucial performance hacks in the interfacing.
    enum DNA_Base: uint8_t
    {
        A = 0b00,   // 0b00
        C = 0b01,   // 0b01
        G = 0b10,   // 0b11
        T = 0b11,   // 0b11
        N = 0b100   // 0b100
    };

    // Returns the mapping integer value of the given character `nucleotide`.
    static DNA_Base map_nucleotide(const char nucleotide);

    // Returns the mapping integer value of the complement of `nucleotide`.
    static DNA_Base complement_nucleotide(const DNA_Base nucleotide);


public:

    // Returns the DNA-complement character of the character `nucl`.
    static char complement(const char nucl);
};


inline Kmer_Utility::DNA_Base Kmer_Utility::map_nucleotide(const char nucleotide)
{
    switch(nucleotide)
    {
    case 'A':
        return DNA_Base::A;
    
    case 'C':
        return DNA_Base::C;

    case 'G':
        return DNA_Base::G;

    case 'T':
        return DNA_Base::T;

    default:
        // Placeholder rule to handle `N` nucleotides. Currently, as per the rule used by the KMC tool.
        
        std::cerr << "Encountered invalid nucleotide " << nucleotide << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline Kmer_Utility::DNA_Base Kmer_Utility::complement_nucleotide(const DNA_Base nucleotide)
{
    switch(nucleotide)
    {
    case DNA_Base::A:
        return DNA_Base::T;

    case DNA_Base::C:
        return DNA_Base::G;

    case DNA_Base::G:
        return DNA_Base::C;

    case DNA_Base::T:
        return DNA_Base::A;

    default:
        // Placeholder rule to handle `N` nucleotides. Currently, as per the rule used by the KMC tool.
        
        std::cerr << "Encountered invalid DNA_Base. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline char Kmer_Utility::complement(const char nucl)
{
    switch (nucl)
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
        std::cerr << "Invalid nucleotide " << nucl << " encountered. Aborting.";
        std::exit(EXIT_FAILURE);
    }
}




#endif
