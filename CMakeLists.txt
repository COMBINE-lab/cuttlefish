
# Specify the minimum version of CMake to use.
# CMake can have different behaviors (policies) based on its version used.
cmake_minimum_required(VERSION 3.14)


# Name the project, its version, and languages used in it.
set(PROJECT_NAME cuttlefish)
project(${PROJECT_NAME}
        VERSION 0.1.0
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
# the `kseq` and the `kmc` libraries to gzip-compressed files.
include(FindZLIB)
if(NOT ZLIB_FOUND)
    message(FATAL_ERROR "zlib (https://zlib.net/) is required. Aborting.")
endif()

# Search and load setting for the `bzip2` library. It is required to seamlessly adapt the
# `kmc` library to bzip-compressed files.
include(FindBZip2)
if(NOT BZIP2_FOUND)
    message(FATAL_ERROR "bzip2 (https://sourceware.org/bzip2/) is required. Aborting.")
endif()


# Module required to download and install external projects.
include(ExternalProject)
set(EXT_LIB ${CMAKE_SOURCE_DIR}/external/lib/)
set(EXT_INCLUDE ${CMAKE_SOURCE_DIR}/external/include/)


# Prepare the `jemalloc` library. It provides scalable concurrency support and better avoidance
# of fragmentation.
message("Build system will fetch and install jemalloc")
ExternalProject_Add(prj_jemalloc
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/jemalloc/jemalloc/archive/5.2.1.tar.gz -o jemalloc-5.2.1.tar.gz &&
                        tar -xzf jemalloc-5.2.1.tar.gz &&
                        rm jemalloc-5.2.1.tar.gz

    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/jemalloc-5.2.1
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   sh -c "CC=${CMAKE_C_COMPILER} ./autogen.sh --disable-debug --enable-static --prefix=<INSTALL_DIR> --silent"
    INSTALL_COMMAND     cp -r lib ${EXT_LIB} &&
                        cp -r include ${EXT_INCLUDE}
)

add_library(jemalloc STATIC IMPORTED)
set_target_properties(jemalloc PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libjemalloc.a)
add_dependencies(jemalloc prj_jemalloc)


# Prepare the `kmc` library — required by the Cuttlefish algorithm implementation.
message("Build system will fetch and install KMC3")
ExternalProject_Add(prj_kmc
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    # DOWNLOAD_COMMAND    curl -k -L https://github.com/refresh-bio/KMC/archive/refs/heads/master.zip -o KMC-master.zip &&
    #                     unzip -qq KMC-master.zip &&
    #                     rm KMC-master.zip
    DOWNLOAD_COMMAND    git clone https://github.com/refresh-bio/kmc &&
                        cd kmc &&
                        git checkout f36cf7ca452c08e5fbb694bde3ee7d430f7bdb7a

    # SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/KMC-master
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/kmc
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external/
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       make -j16
    INSTALL_COMMAND     cp bin/libkmc_core.a ${EXT_LIB} &&
                        cp include/kmc_runner.h ${EXT_INCLUDE}
)

add_library(kmc STATIC IMPORTED)
set_target_properties(kmc PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libkmc_core.a)
add_dependencies(kmc prj_kmc)


# The `Debug` configuration optimizes the program for debugging and enables full debug information.
# The `Release` configuration enables most compiler optimizations for speed and defines `NDEBUG`
# (No Debug) which will remove all traces of the standard library assert calls.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


# Add subdirectory `src` to the build; CMake will open `src/CMakeLists.txt` for such.
add_subdirectory(src)
