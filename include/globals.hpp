
#ifndef GLOBALS_HPP
#define GLOBALS_HPP



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
    typedef uint8_t vertex_code_t;

    constexpr nucleotide_t PLACEHOLDER_NUCLEOTIDE = 'N';

    constexpr kmer_dir_t FWD = true;
    constexpr kmer_dir_t BWD = false;

    enum class Vertex_Class: uint8_t
    {
        single_in_single_out = 0,
        multi_in_single_out = 1,
        single_in_multi_out = 2,
        multi_in_multi_out = 3
    };

    constexpr uint8_t BITS_PER_KMER = 5;
    typedef compact::cas_vector<uint8_t, BITS_PER_KMER, uint64_t, std::allocator<uint64_t>> bitvector_t;
    typedef compact::iterator_imp::lhs_setter<uint8_t, BITS_PER_KMER, uint64_t, true, 63U> bitvector_entry_t;
}


// Returns the plain DNA-complement character of the provided `nucleotide` character.
// TODO: Move it to the consumer classes (only `CdBG_Builder`).
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
