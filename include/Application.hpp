
#ifndef APPLICATION_HPP
#define APPLICATION_HPP



#include "CdBG.hpp"
#include "Build_Params.hpp"


// The top-level application class for the compaction algorithm.
template <uint16_t k>
class Application
{
private:

    // Pointer to an application instance of the next `Application` class in the top-down hierarchy (on `k`).
    Application<k - 2>* const app_next_level;

    // Pointer to a `CdBG` object that operates with the k-value `k`.
    CdBG<k>* const cdbg;


public:

    // Constructs an `Application` instance with the provided build-parameters,
    // if the provided `k` parameter matches to the specialized template argument `k`.
    Application(const Build_Params& params);

    ~Application();

    // Executes the compaction algorithm.
    void execute() const;
};


template <>
class Application<1>
{
private:

    CdBG<1>* const cdbg;


public:

    Application(const Build_Params& params):
        cdbg(params.k() == 1 ? new CdBG<1>(params) : nullptr)
    {}

    ~Application()
    {
        if(cdbg != nullptr)
            delete cdbg;
    }

    void execute() const
    {
        if(cdbg != nullptr)
            cdbg->construct();
    }
};


template <uint16_t k>
inline Application<k>::Application(const Build_Params& params):
    app_next_level(new Application<k - 2>(params)),
    cdbg(params.k() == k ? new CdBG<k>(params) : nullptr)
{}


template <uint16_t k>
inline Application<k>::~Application()
{
    delete app_next_level;

    if(cdbg != nullptr)
        delete cdbg;
}


template <uint16_t k>
inline void Application<k>::execute() const
{
    if(cdbg != nullptr)
        cdbg->construct();
    else
        app_next_level->execute();
}



#endif
