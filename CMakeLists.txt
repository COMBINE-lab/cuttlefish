
# Specify the minimum version of CMake to use.
# CMake can have different behaviors (policies) based on its version used.
cmake_minimum_required(VERSION 3.14)


# Name the project, its version, and languages used in it.
set(PROJECT_NAME cuttlefish)
project(${PROJECT_NAME}
        VERSION 0.1.0
        LANGUAGES CXX
        )


# Fix language standards, and set hard requirements for such.
# All targets defined from this point onward will pick up these requirements.
set(CMAKE_CXX_STANDARD 11)


# Bundle the warning flags that we want to pass on to the compiler.
# Disable unknown pragmas, b/c bbhash uses them extensively
# Same with maybe-uninitialized (might want this one back later though)
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
set(WARNING_FLAGS -Wall -Wextra -Wno-unknown-pragmas -Wno-maybe-uninitialized)
add_definitions(-D__STDC_FORMAT_MACROS)

# Search the file system for the appropriate threads package for this platform, and then set
# the `CMAKE_THREAD_LIBS_INIT` variable (and some other variables as well).
find_package(Threads REQUIRED)  # The threads package is required for the BBHash library used in the project.
set(THREADS_PREFER_PTHREAD_FLAG TRUE)   # The BBHash library uses `pthread`.


# The `Debug` configuration optimizes the program for debugging and enables full debug information.
# The `Release` configuration enables most compiler optimizations for speed and defines `NDEBUG`
# (No Debug) which will remove all traces of the standard library assert calls.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


# Add subdirectory `src` to the build; CMake will open `src/CMakeLists.txt` for such.
add_subdirectory(src)
