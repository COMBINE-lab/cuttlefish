
#include "Application.hpp"
#include "globals.hpp"
#include "CdBG.hpp"
#include "Read_CdBG.hpp"
#include "dBG_Contractor.hpp"
#include "Kmer_Index.hpp"


template <uint16_t k, template <uint16_t> typename T_App>
Application<k, T_App>::Application(const Build_Params& params):
#ifndef FIXED_K
    app_next_level(new Application<k - 2, T_App>(params))
#else
    app_next_level(new Application<1, T_App>(params))
#endif
    , app(params.k() == k ? new T_App<k>(params) : nullptr)
    , validator(nullptr)
{}


/*
template <uint16_t k, template <uint16_t> typename T_App>
Application<k, T_App>::Application(const Validation_Params& params):
#ifndef FIXED_K
    app_next_level(new Application<k - 2, T_App>(params))
#else
    app_next_level(new Application<1, T_App>(params))
#endif
    , app(nullptr)
    , validator(params.k() == k ? new Validator<k>(params): nullptr)
{}
*/


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


/*
template <uint16_t k, template <uint16_t> typename T_App>
bool Application<k, T_App>::validate() const
{
    if(validator != nullptr)
        return validator->validate();

    return app_next_level->validate();
}
*/



// Template instantiations for the required instances.
// template class Application<cuttlefish::MAX_K, CdBG>;
// template class Application<cuttlefish::MAX_K, Read_CdBG>;
// template class Application<cuttlefish::MAX_K, Kmer_Index>;
template class Application<cuttlefish::MAX_K, cuttlefish::dBG_Contractor>;
