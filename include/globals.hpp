
#ifndef GLOBALS_HPP
#define GLOBALS_HPP



#include <memory>
#include <iostream>


// Forward declaration of k-mer type.
class Kmer;

// Forward declarations of the type of the bitvector used and the type to access its entries (mutable).
namespace compact
{
    template<typename IDX, unsigned BITS, typename W, typename Allocator> class cas_vector;

    namespace iterator_imp
    {
        template<typename IDX, unsigned BITS, typename W, bool TS, unsigned UB> class lhs_setter;
    }
}


namespace cuttlefish
{
    typedef Kmer kmer_t;
    typedef bool kmer_dir_t;
    typedef char nucleotide_t;
    typedef uint8_t state_t;
    typedef uint8_t vertex_code_t;

    constexpr nucleotide_t PLACEHOLDER_NUCLEOTIDE = 'N';

    constexpr kmer_dir_t FWD = true;
    constexpr kmer_dir_t BWD = false;

    constexpr uint8_t SINGLE_IN_SINGLE_OUT = 0;
    constexpr uint8_t MULTI_IN_SINGLE_OUT = 1;
    constexpr uint8_t SINGLE_IN_MULTI_OUT = 2;
    constexpr uint8_t MULTI_IN_MULTI_OUT = 3;

    constexpr uint8_t BITS_PER_KMER = 5;
    typedef compact::cas_vector<uint8_t, BITS_PER_KMER, uint64_t, std::allocator<uint64_t>> bitvector_t;
    typedef compact::iterator_imp::lhs_setter<uint8_t, BITS_PER_KMER, uint64_t, true, 63U> bitvector_entry_t;
}


// Returns the plain DNA-complement character of the provided `nucleotide` character.
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
        
        std::cerr << "Invalid nucleotide " << nucleotide << " encountered. Aborting.";
        std::exit(EXIT_FAILURE);
    }
}


#endif
