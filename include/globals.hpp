
#ifndef GLOBALS_HPP
#define GLOBALS_HPP



#include "Kmer.hpp"

#include <memory>
#include "boost/preprocessor/repetition/repeat.hpp"


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
    constexpr uint16_t MAX_K = (2 * INSTANCE_COUNT - 1);


    typedef bool dir_t;
    typedef DNA::Base base_t;
    typedef uint8_t state_code_t;


    constexpr dir_t FWD = true;
    constexpr dir_t BWD = false;


    enum class State_Class: uint8_t
    {
        single_in_single_out = 0,
        multi_in_single_out = 1,
        single_in_multi_out = 2,
        multi_in_multi_out = 3,
    };


    typedef enum class Side: bool
    {
        front = false,
        back = true
    } side_t;


    constexpr uint8_t BITS_PER_REF_KMER = 5;
    typedef compact::ts_vector<state_code_t, BITS_PER_REF_KMER, uint64_t, std::allocator<uint64_t>> ref_bitvector_t;
    typedef compact::iterator_imp::lhs_setter<state_code_t, BITS_PER_REF_KMER, uint64_t, true, 64U> ref_bitvector_entry_t;

    constexpr uint8_t BITS_PER_READ_KMER = 6;


    typedef std::shared_ptr<spdlog::logger> logger_t;
}


// The macro `INSTANCE_COUNT` must be set exactly to `(MAX_K + 1) / 2` for a required maximum k-value.
// Also, the `MAX_K` value must be odd (as the k-values used in the algorithm) for correct results.
#ifndef INSTANCE_COUNT
    #define INSTANCE_COUNT 32
#endif

#define INSTANTIATE(z, k, class_name) template class class_name<2 * k + 1>;
#define ENUMERATE(count, instantiator, class_name) BOOST_PP_REPEAT(count, instantiator, class_name)

#define INSTANTIATE_PER_BIT(z, k, class_name)   template class class_name<2 * k + 1, cuttlefish::BITS_PER_REF_KMER>;\
                                                template class class_name<2 * k + 1, cuttlefish::BITS_PER_READ_KMER>;
#define ENUMERATE_PER_BIT(count, instantiator, class_name) BOOST_PP_REPEAT(count, instantiator, class_name)

// BOOST_PP_REPEAT reference: https://www.boost.org/doc/libs/1_55_0/libs/preprocessor/doc/ref/repeat.html



#endif
