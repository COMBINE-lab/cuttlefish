
#ifndef APPLICATION_HPP
#define APPLICATION_HPP



#include "CdBG.hpp"
#include "Validator.hpp"
#include "Build_Params.hpp"
#include "Validation_Params.hpp"


// The top-level application class for the compaction algorithm.
template <uint16_t k>
class Application
{
private:

    // Pointer to an application instance of the next `Application` class in the top-down hierarchy (on `k`).
    Application<k - 2>* const app_next_level;

    // Pointer to a `CdBG` object that operates with the k-value `k`.
    CdBG<k>* const cdbg;

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


template <>
class Application<1>
{
private:

    CdBG<1>* const cdbg;

    Validator<1>* const validator;


public:

    Application(const Build_Params& params):
        cdbg(params.k() == 1 ? new CdBG<1>(params) : nullptr),
        validator(nullptr)
    {}


    Application(const Validation_Params& params):
        cdbg(nullptr),
        validator(params.k() == 1 ? new Validator<1>(params) : nullptr)
    {}


    ~Application()
    {
        if(cdbg != nullptr)
            delete cdbg;

        if(validator != nullptr)
            delete validator;
    }


    void execute() const
    {
        if(cdbg != nullptr)
            cdbg->construct();
        else
        {
            std::cerr << "The provided k is not valid. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }


    bool validate() const
    {
        if(validator != nullptr)
            return validator->validate();

        std::cerr << "The provided k is not valid. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
};


template <uint16_t k>
inline Application<k>::Application(const Build_Params& params):
    app_next_level(new Application<k - 2>(params)),
    cdbg(params.k() == k ? new CdBG<k>(params) : nullptr),
    validator(nullptr)
{}


template <uint16_t k>
inline Application<k>::Application(const Validation_Params& params):
    app_next_level(new Application<k - 2>(params)),
    cdbg(nullptr),
    validator(params.k() == k ? new Validator<k>(params): nullptr)
{}


template <uint16_t k>
inline Application<k>::~Application()
{
    delete app_next_level;

    if(cdbg != nullptr)
        delete cdbg;

    if(validator != nullptr)
        delete validator;
}


template <uint16_t k>
inline void Application<k>::execute() const
{
    if(cdbg != nullptr)
        cdbg->construct();
    else
        app_next_level->execute();
}


template <uint16_t k>
inline bool Application<k>::validate() const
{
    if(validator != nullptr)
        return validator->validate();

    return app_next_level->validate();
}



#endif
