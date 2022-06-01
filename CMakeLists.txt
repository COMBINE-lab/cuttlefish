

# Specify the minimum version of CMake to use.
# CMake can have different behaviors (policies) based on its version used.
cmake_minimum_required(VERSION 3.14)


# Name the project, its version, and languages used in it.
set(PROJECT_NAME cuttlefish)
project(${PROJECT_NAME}
        VERSION 2.0.0
        LANGUAGES CXX C
        )


# Fix language standards, and set hard requirements for such.
# All targets defined from this point onward will pick up these requirements.
set(CMAKE_CXX_STANDARD 17)
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
# TODO: find out what are `__STDC_FORMAT_MACROS` and `SPDLOG_FMT_EXTERNAL_HO` for.
add_compile_definitions(__STDC_FORMAT_MACROS SPDLOG_FMT_EXTERNAL_HO FMT_HEADER_ONLY XXH_INLINE_ALL)

add_compile_definitions(PROJECT_VERSION=${CMAKE_PROJECT_VERSION})

if(INSTANCE_COUNT)
    add_compile_definitions(INSTANCE_COUNT=${INSTANCE_COUNT})
endif()

if(CF_VALIDATION_MODE)
    add_compile_definitions(CF_VALIDATION_MODE)
endif()

if(CF_DEVELOP_MODE)
    add_compile_definitions(CF_DEVELOP_MODE)
endif()


# Here, we have some platform-specific considerations
# of which we must take care.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    execute_process(
        COMMAND getconf LEVEL1_DCACHE_LINESIZE
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    )
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # OSX has `getconf`, but doesn't report `LEVEL1_DCACHE_LINESIZE`;
    # so we instead use `sysctl` for the corresponding variable.
    execute_process(
        COMMAND sysctl machdep.cpu.cache.linesize
        COMMAND awk "{print $2}"
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    )

    # Later on, jemalloc will complain if the C compiler 
    # hasn't been properly set.
    if(NOT CONDA_BUILD)
        set(CMAKE_C_COMPILER clang)
        set(CMAKE_CXX_COMPILER clang++)
    endif()
endif()

add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})


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


# Set path for modules required to search for existing packages in the system.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# Module required to download and install external projects.
include(ExternalProject)

set(EXT_LIB ${CMAKE_SOURCE_DIR}/external/lib/)
set(EXT_INCLUDE ${CMAKE_SOURCE_DIR}/external/include/)

file(MAKE_DIRECTORY ${EXT_LIB} ${EXT_INCLUDE})


# Prepare the `jemalloc` library. It provides scalable concurrency support and better avoidance
# of fragmentation.
set(MALLOC_LIB "")
set(JEMALLOC_MIN_VERSION "5.2.1")
find_package(Jemalloc ${JEMALLOC_MIN_VERSION})
if(JEMALLOC_FOUND)
    message("Found jemalloc (v${JEMALLOC_VERSION}) in the system")
    set(MALLOC_LIB ${JEMALLOC_LIBRARIES})
else()
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
        INSTALL_COMMAND     cp lib/libjemalloc.a ${EXT_LIB}
    )

    set(MALLOC_LIB ${EXT_LIB}/libjemalloc.a)
endif()

add_library(jemalloc STATIC IMPORTED)
set_target_properties(jemalloc PROPERTIES IMPORTED_LOCATION ${MALLOC_LIB})
if(NOT JEMALLOC_FOUND)
    add_dependencies(jemalloc prj_jemalloc)
endif()


# Prepare the `kmc` library — required by the Cuttlefish algorithm implementation.
message("Build system will fetch and install KMC3")
ExternalProject_Add(prj_kmc
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/refresh-bio/KMC/archive/refs/tags/v3.2.1.tar.gz -o KMC-3.2.1.tar.gz &&
                        tar -xzf KMC-3.2.1.tar.gz &&
                        rm KMC-3.2.1.tar.gz
    PATCH_COMMAND       patch --strip 1 < ${CMAKE_SOURCE_DIR}/patches/kmc_patch.diff
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/KMC-3.2.1
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external/
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       make -j CC=${CMAKE_CXX_COMPILER} kmc
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

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3")


# Help `conda` build for OSX through circumventing some of its old SDK-based checks.
if(CONDA_BUILD)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D_LIBCPP_DISABLE_AVAILABILITY")
endif()


# Add subdirectory `src` to the build; CMake will open `src/CMakeLists.txt` for such.
add_subdirectory(src)
