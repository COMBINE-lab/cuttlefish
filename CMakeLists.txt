
# Specify the minimum version of CMake to use.
# CMake can have different behaviors (policies) based on its version used.
cmake_minimum_required(VERSION 3.14)


# Name the project, its version, and languages used in it.
set(PROJECT_NAME cuttlefish)
project(${PROJECT_NAME}
        VERSION 1.0.0
        LANGUAGES CXX C
        )


# Fix language standards, and set hard requirements for such.
# All targets defined from this point onward will pick up these requirements.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Turn off using platform-specific compiler standards.
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_EXTENSIONS OFF)


# Bundle the warning flags that we want to pass on to the compiler.
# Disable unknown pragmas, b/c bbhash uses them extensively.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
set(WARNING_FLAGS -Wall -Wextra)
set(SUPPRESS_WARNING_FLAGS -Wno-unknown-pragmas)

# Bundle the extra optimization flags (not associated with the `-O` levels)
# that we want to pass on to the compiler.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
set(OPTIMIZE_FLAGS -funroll-loops)


# Add the required preprocessor definitions (`#define`s) to pass on.
add_definitions(-D__STDC_FORMAT_MACROS -DSPDLOG_FMT_EXTERNAL_HO -DFMT_HEADER_ONLY -DXXH_INLINE_ALL)
if(INSTANCE_COUNT)
    add_definitions(-DINSTANCE_COUNT=${INSTANCE_COUNT})
endif()


# Search the file system for the appropriate threads package for this platform, and then set
# the `CMAKE_THREAD_LIBS_INIT` variable (and some other variables as well).
find_package(Threads REQUIRED)  # The threads package is required for the BBHash library used in the project.
set(THREADS_PREFER_PTHREAD_FLAG TRUE)   # The BBHash library uses `pthread`.

# Search and load setting for the `zlib` library. The library is required to seamlessly adapt
# the `kseq` library to gzipped compressed files.
find_package(ZLIB REQUIRED)


# The `Debug` configuration optimizes the program for debugging and enables full debug information.
# The `Release` configuration enables most compiler optimizations for speed and defines `NDEBUG`
# (No Debug) which will remove all traces of the standard library assert calls.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3")


# Add subdirectory `src` to the build; CMake will open `src/CMakeLists.txt` for such.
add_subdirectory(src)
