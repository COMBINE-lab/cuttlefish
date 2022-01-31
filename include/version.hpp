

#include <string>


// https://stackoverflow.com/a/20632065/2007834
#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)


inline std::string version()
{
    return STRINGIFY(PROJECT_VERSION);
}
