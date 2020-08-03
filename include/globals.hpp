
#ifndef GLOBALS_HPP
#define GLOBALS_HPP



#include <memory>


// Forward declaration of the k-mer type.
class Kmer_u64;

// Forward declaration of the k-mer hasher type.
class Kmer_Hasher;

// Forward declaration of the minimal perfect hash function type.
namespace boomphf
{
    template<typename elem_t, typename Hasher_t> class mphf;
}

// Forward declarations of the type of the bitvector used and the type to access its entries (mutable).
namespace compact
{
    template<typename IDX, unsigned BITS, typename W, typename Allocator> class ts_vector;

    namespace iterator_imp
    {
        template<typename IDX, unsigned BITS, typename W, bool TS, unsigned UB> class lhs_setter;
    }
}

// Forward declaration of the output writer type.
namespace spdlog
{
    class logger;
}


namespace cuttlefish
{
    typedef Kmer_u64 kmer_t;
    typedef bool dir_t;
    typedef char nucleotide_t;
    typedef uint8_t state_code_t;

    constexpr nucleotide_t PLACEHOLDER_NUCLEOTIDE = 'N';

    constexpr dir_t FWD = true;
    constexpr dir_t BWD = false;

    enum class Vertex_Class: uint8_t
    {
        single_in_single_out = 0,
        multi_in_single_out = 1,
        single_in_multi_out = 2,
        multi_in_multi_out = 3
    };

    typedef boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> mphf_t;    // The MPH function type.

    constexpr uint8_t BITS_PER_KMER = 5;
    typedef compact::ts_vector<state_code_t, BITS_PER_KMER, uint64_t, std::allocator<uint64_t>> bitvector_t;
    typedef compact::iterator_imp::lhs_setter<state_code_t, BITS_PER_KMER, uint64_t, true, 64U> bitvector_entry_t;

    typedef std::shared_ptr<spdlog::logger> logger_t;
}



#endif
