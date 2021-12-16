
#include "Application.hpp"
#include "globals.hpp"
#include "CdBG.hpp"
#include "Read_CdBG.hpp"


template <uint16_t k, template <uint16_t> typename T_App>
Application<k, T_App>::Application(const Build_Params& params):
    app_next_level(new Application<k - 2, T_App>(params)),
    app(params.k() == k ? new T_App<k>(params) : nullptr),
    validator(nullptr)
{}


template <uint16_t k, template <uint16_t> typename T_App>
Application<k, T_App>::Application(const Validation_Params& params):
    app_next_level(new Application<k - 2, T_App>(params)),
    app(nullptr),
    validator(params.k() == k ? new Validator<k>(params): nullptr)
{}


template <uint16_t k, template <uint16_t> typename T_App>
Application<k, T_App>::~Application()
{
    delete app_next_level;

    if(app != nullptr)
        delete app;

    if(validator != nullptr)
        delete validator;
}


template <uint16_t k, template <uint16_t> typename T_App>
void Application<k, T_App>::execute() const
{
    if(app != nullptr)
        app->construct();
    else
        app_next_level->execute();
}


template <uint16_t k, template <uint16_t> typename T_App>
bool Application<k, T_App>::validate() const
{
    if(validator != nullptr)
        return validator->validate();

    return app_next_level->validate();
}


template <template<uint16_t> typename T_App>
Application<1, T_App>::Application(const Build_Params& params):
    app(params.k() == 1 ? new T_App<1>(params) : nullptr),
    validator(nullptr)
{}


template <template<uint16_t> typename T_App>
Application<1, T_App>::Application(const Validation_Params& params):
    app(nullptr),
    validator(params.k() == 1 ? new Validator<1>(params) : nullptr)
{}


template <template<uint16_t> typename T_App>
Application<1, T_App>::~Application()
{
    if(app != nullptr)
        delete app;

    if(validator != nullptr)
        delete validator;
}


template <template<uint16_t> typename T_App>
void Application<1, T_App>::execute() const
{
    if(app != nullptr)
        app->construct();
    else
    {
        std::cerr << "The provided k is not valid. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <template<uint16_t> typename T_App>
bool Application<1, T_App>::validate() const
{
    if(validator != nullptr)
        return validator->validate();

    std::cerr << "The provided k is not valid. Aborting.\n";
    std::exit(EXIT_FAILURE);
}



// Template instantiations for the required instances.
template class Application<cuttlefish::MAX_K, CdBG>;
template class Application<cuttlefish::MAX_K, Read_CdBG>;
