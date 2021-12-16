
#ifndef APPLICATION_HPP
#define APPLICATION_HPP



#include "Validator.hpp"
#include "Build_Params.hpp"
#include "Validation_Params.hpp"


// The top-level application class for the compaction algorithm.
template <uint16_t k, template<uint16_t> typename T_App>
class Application
{
private:

    // Pointer to an application instance of the next `Application` class in the top-down hierarchy (on `k`).
    Application<k - 2, T_App>* const app_next_level;

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
private:

    T_App<1>* const app;

    Validator<1>* const validator;


public:

    Application(const Build_Params& params);

    Application(const Validation_Params& params);

    ~Application();

    void execute() const;

    bool validate() const;
};



#endif
