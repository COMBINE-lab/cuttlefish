
#ifndef APPLICATION_HPP
#define APPLICATION_HPP



#include "Validator.hpp"

#include <cstdint>


class Build_Params;
class Validation_Params;


// The top-level application class for the compaction algorithm.
template <uint16_t k, template<uint16_t> typename T_App>
class Application
{
private:

    // Pointer to an application instance of the next `Application` class in the top-down hierarchy (on `k`).
#ifndef FIXED_K
    Application<k - 2, T_App>* const app_next_level;
#else
    Application<1, T_App>* const app_next_level;
#endif


    // Pointer to a driver object that operates with the k-value `k`.
    T_App<k>* const app;

    // TODO: Make the validator member generic, like `T_App`.
    // Pointer to a `Validator` object that operates with the k-value `k`.
    Validator<k>* const validator;


public:

    // Constructs an `Application` instance with the provided build-parameters,
    // if the provided `k` parameter matches to the specialized template argument `k`.
    Application(const Build_Params& params);

    // Constructs an `Application` instance with the provided validation-parameters,
    // if the provided `k` parameter matches to the specialized template argument `k`.
    Application(const Validation_Params& params);

    ~Application();

    // Executes the compaction algorithm.
    void execute() const;

    // Validates the result of the compaction algorithm.
    bool validate() const;
};


template <template<uint16_t> typename T_App>
class Application<1, T_App>
{
public:

    Application(const Build_Params&) {}
    Application(const Validation_Params&) {}

    void execute() const {}
    bool validate() const { return false; }
};



#endif
