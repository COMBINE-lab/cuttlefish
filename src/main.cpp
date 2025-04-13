
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>


#ifdef __cplusplus
extern "C" {
#endif
    int cf_build(int argc, char** argv);
    int cf_validate(int argc, char** argv);
    int print_cf_version();
#ifdef __cplusplus
}
#endif


void display_help_message()
{
    print_cf_version();
    std::cout << "Supported commands: `build`, `help`, `version`.\n";

    std::cout << "Usage:\n";
    std::cout << "\tcuttlefish build [options]\n";
}


int main(int argc, char** argv)
{
#ifdef CF_DEVELOP_MODE
    std::cout << "Warning: Executing in Develop Mode.\n";
#endif

#ifndef NDEBUG
    std::cout << "Warning: Executing in Debug Mode.\n";
#endif

    if(argc < 2)
        display_help_message();
    else
    {
        std::string command(argv[1]);
        std::transform(command.begin(), command.end(), command.begin(), [](const char ch) { return std::tolower(ch); });

        if(command == "build")
            return cf_build(argc - 1, argv + 1);
        // else if(command == "validate")
        //     return cf_validate(argc - 1, argv + 1);
        else if(command == "help")
            display_help_message();
        else if(command == "version")
            return print_cf_version();
        else
            display_help_message();
    }

    return EXIT_SUCCESS;
}
