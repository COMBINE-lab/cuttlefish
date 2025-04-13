
# =============================================================================

# Specify the minimum version of CMake to use.
# CMake can have different behaviors (policies) based on its version used.
cmake_minimum_required(VERSION 3.17)


# Name the project, its version, and languages used in it.
set(PROJECT_NAME cuttlefish)
project(${PROJECT_NAME}
        VERSION 2.2.0
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

# Fix minimum compiler versions.
set(GCC_VERSION_MIN "9.1")
set(CLANG_VERSION_MIN "9.0")

if(ASAN OR MSAN OR UBSAN)
    set(CLANG TRUE)
endif()

if(CLANG)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
    # set(CMAKE_LINKER_TYPE LLD)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_VERSION_MIN)
        message(FATAL_ERROR "${PROJECT_NAME} requires GCC version to be at least ${GCC_VERSION_MIN}."
                            " Available version is ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS CLANG_VERSION_MIN)
        message(FATAL_ERROR "${PROJECT_NAME} requires Clang version to be at least ${CLANG_VERSION_MIN}."
                            " Available version is ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
endif()


# Bundle the warning flags that we want to pass on to the compiler.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
set(WARNING_FLAGS -Wall -Wextra -Wpedantic)
set(SUPPRESS_WARNING_FLAGS -Wno-unknown-pragmas -Wno-pragmas -Wno-unused-function)

# Warning flags to be suppressed with clang.
set(CLANG_SUPPRESS_WARNING_FLAGS    -Wno-unknown-attributes
                                    -Wno-nested-anon-types -Wno-gnu-anonymous-struct
                                    -Wno-c++20-attribute-extensions
                                    -Wno-gnu-zero-variadic-macro-arguments)

if(CLANG)
    list(APPEND SUPPRESS_WARNING_FLAGS ${CLANG_SUPPRESS_WARNING_FLAGS})
endif()

# Bundle debugging information flags that we want to conditionally
# pass on to the compiler.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html.
set(DEBUG_FLAGS -ggdb)

# Bundle address sanitization flags: https://github.com/google/sanitizers/wiki/AddressSanitizer.
set(ASAN_FLAGS -fsanitize=address -fno-omit-frame-pointer)
# set(ASAN_LINK_FLAGS -static-libasan)

# Bundle memory sanitization flags: https://github.com/google/sanitizers/wiki/MemorySanitizer.
set(MSAN_FLAGS -fsanitize=memory -fPIE -fno-omit-frame-pointer -fsanitize-memory-track-origins)
set(MSAN_LINK_FLAGS -pie)

# Bundle undefined behavior sanitization flags: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html.
set(UBSAN_FLAGS -fsanitize=undefined
                -fsanitize=float-divide-by-zero #-fsanitize=unsigned-integer-overflow
                -fsanitize=implicit-conversion -fsanitize=local-bounds -fsanitize=nullability)
# set(UBSAN_LINK_FLAGS -fuse-ld=lld)

# Bundle the extra optimization flags (not associated with the `-O` levels)
# that we want to pass on to the compiler.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
set(OPTIMIZE_FLAGS -funroll-loops -march=native -mavx -msse4.1 -msse3)

# Bundle miscellaneous required flags.
# `-fconstexpr-ops-limit` for compile-time computations in `rapidgzip`.
if(CLANG)
    set(MISC_FLAGS -fconstexpr-steps=99000100)
else()
    set(MISC_FLAGS -fconstexpr-ops-limit=99000100)
endif()


# Accumulate the compile- and link-flags.
set(COMPILE_FLAGS "")
set(LINK_FLAGS "")

# Specify the warnings and the extra optimization flags to the compiler for the target library.
list(APPEND COMPILE_FLAGS ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS} ${MISC_FLAGS})

if(DEBUG_INFO OR PART_PROFILE)
    list(APPEND COMPILE_FLAGS ${DEBUG_FLAGS})
endif()

if(ASAN)
    list(APPEND COMPILE_FLAGS ${ASAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${ASAN_LINK_FLAGS} ${ASAN_FLAGS})
endif()

if(UBSAN)
    list(APPEND COMPILE_FLAGS ${UBSAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${UBSAN_LINK_FLAGS} ${UBSAN_FLAGS})
endif()

if(MSAN)
    list(APPEND COMPILE_FLAGS ${MSAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${MSAN_LINK_FLAGS} ${MSAN_FLAGS})
endif()


# Add the required preprocessor definitions (`#define`s) to pass on.
add_compile_definitions(PROJECT_VERSION=${CMAKE_PROJECT_VERSION})
add_compile_definitions(FMT_HEADER_ONLY)    # For `spdlog`.
add_compile_definitions(WITH_ISAL)  # For `rapidgzip`.
add_compile_definitions(WITH_RPMALLOC)  # For `rapidzip`.
add_compile_definitions(USE_IGZIP)  # TODO: temporary for `rabbitfx`.
add_compile_definitions(CMAKE_EXPORT_COMPILE_COMMANDS)  # For better `clangd` scanning.
add_compile_definitions(USE_PARLAY) # For `parlayhash`.

if(INSTANCE_COUNT)
    add_compile_definitions(INSTANCE_COUNT=${INSTANCE_COUNT})
endif()

if(FIXED_K)
    add_compile_definitions(FIXED_K=${FIXED_K})
endif()

if(CF_VALIDATION_MODE)
    add_compile_definitions(CF_VALIDATION_MODE)
endif()

if(CF_DEVELOP_MODE)
    add_compile_definitions(CF_DEVELOP_MODE)
endif()

# if(IGZIP_PREFIX)
#     add_compile_definitions(USE_IGZIP)
# endif()

if(PART_PROFILE)
    add_compile_definitions(PART_PROFILE)
endif()


# Here, we have some platform-specific considerations of which we must take
# care.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    execute_process(
        COMMAND getconf LEVEL1_DCACHE_LINESIZE
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    )
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
        # OSX has `getconf`, but doesn't report `LEVEL1_DCACHE_LINESIZE`;
        # so we instead use `sysctl` for the corresponding variable.
        execute_process(
            COMMAND sysctl machdep.cpu.cache.linesize
            COMMAND awk "{print $2}"
            COMMAND tr -d '\n'
            OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
        )
    elseif(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
        execute_process(
            COMMAND sysctl hw.cachelinesize
            COMMAND awk "{print $2}"
            COMMAND tr -d '\n'
            OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
        )
    else()
        message(FATAL_ERROR "Unable to identify Apple Silicon processor architecture")
    endif()

    # Later on, jemalloc will complain if the C compiler  hasn't been properly set.
    if(NOT CONDA_BUILD)
        set(CMAKE_C_COMPILER clang)
        set(CMAKE_CXX_COMPILER clang++)
    endif()
endif()

add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})


# Set path for modules required to search for existing packages in the system.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# Search the file system for the appropriate threads package for this platform, and then set
# the `CMAKE_THREAD_LIBS_INIT` variable (and some other variables as well).
find_package(Threads REQUIRED)

# Search and load setting for the `zlib` library. The library is required to seamlessly adapt
# the `kseq`, the `kmc`, and the `RabbitFX` libraries to gzip-compressed files.
find_package(ZLIB REQUIRED)

# Prepare the `zstd` library. It provides fast lossless compression, targeting real-time
# compression scenarios at zlib-level and better compression ratios.
# find_package(ZSTD REQUIRED)


# Module required to download and install external projects.
include(ExternalProject)

set(EXT_LIB ${CMAKE_SOURCE_DIR}/external/lib/)
set(EXT_INCLUDE ${CMAKE_SOURCE_DIR}/external/include/)

file(MAKE_DIRECTORY ${EXT_LIB} ${EXT_INCLUDE})

# Build instrumented `libc++` to avoid false-positives from memory sanitizers.
# Ref for some of the fixes: https://github.com/google/sanitizers/issues/1685#issuecomment-1734197581
if(MSAN)
    message("Build system will fetch and install instrumented libc++.")
    set(LLVM_RELEASE llvmorg-17.0.1)
    set(LLVM_SRC ${CMAKE_SOURCE_DIR}/external/llvm-project-${LLVM_RELEASE})
    ExternalProject_Add(prj_libcxx
        DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
        DOWNLOAD_COMMAND    curl -k -L https://github.com/llvm/llvm-project/archive/refs/tags/${LLVM_RELEASE}.zip -o llvm.zip &&
                            unzip -q -o llvm.zip &&
                            rm llvm.zip

        SOURCE_DIR          ${LLVM_SRC}
        BUILD_IN_SOURCE     TRUE
        INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external

        CONFIGURE_COMMAND   cmake -G Ninja -S runtimes -B msan
                            -DCMAKE_BUILD_TYPE=Release
                            -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
                            -DLLVM_ENABLE_RUNTIMES=libcxx$<SEMICOLON>libcxxabi
                            -DLLVM_USE_SANITIZER=MemoryWithOrigins
        BUILD_COMMAND       ninja -C msan
        INSTALL_COMMAND     ""
    )

    set(LIBCXX_INCLUDE ${LLVM_SRC}/msan/include/c++/v1/)
    set(LIBCXX ${LLVM_SRC}/msan/lib/)

    set(MSAN_FLAGS ${MSAN_FLAGS} -stdlib=libc++)
    set(MSAN_LINK_FLAGS ${MSAN_LINK_FLAGS} -stdlib=libc++ -lc++abi -lstdc++)
endif()


# Prepare the `jemalloc` library. It provides scalable concurrency support and better avoidance
# of fragmentation.
# TODO: replace with `rpmalloc` or some other more performant library.
set(MALLOC_LIB "")
set(JEMALLOC_MIN_VERSION "5.3.0")
find_package(Jemalloc ${JEMALLOC_MIN_VERSION})
if(JEMALLOC_FOUND)
    message("Found jemalloc (v${JEMALLOC_VERSION}) in the system")
    set(MALLOC_LIB ${JEMALLOC_LIBRARIES})
else()
    message("Build system will fetch and install jemalloc")
    ExternalProject_Add(prj_jemalloc
        DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
        DOWNLOAD_COMMAND    curl -k -L https://github.com/jemalloc/jemalloc/archive/${JEMALLOC_MIN_VERSION}.tar.gz -o jemalloc-${JEMALLOC_MIN_VERSION}.tar.gz &&
                            tar -xzf jemalloc-${JEMALLOC_MIN_VERSION}.tar.gz &&
                            rm jemalloc-${JEMALLOC_MIN_VERSION}.tar.gz

        SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/jemalloc-${JEMALLOC_MIN_VERSION}
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


# Prepare the `lz4` library. It provides extremely fast compression.
find_package(PkgConfig REQUIRED)
pkg_check_modules(LZ4 REQUIRED liblz4)


# Prepare the `lz4` library. It provides extremely fast compression.
find_package(PkgConfig REQUIRED)
pkg_check_modules(LZ4 REQUIRED liblz4)

# Help `conda` build for OSX through circumventing some of its old SDK-based checks.
if(CONDA_BUILD)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D_LIBCPP_DISABLE_AVAILABILITY")
    if (DEFINED ENV{NUM_JOBS})
      set(KMC_SUBBUILD_THREADS $ENV{NUM_JOBS})
    else()
      set(KMC_SUBBUILD_THREADS 1)
    endif()
elseif(NOT DEFINED KMC_SUBBUILD_THREADS)
    cmake_host_system_information(RESULT KMC_SUBBUILD_THREADS QUERY NUMBER_OF_PHYSICAL_CORES)
endif()

message("Setting KMC_SUBBUILD_THREADS = ${KMC_SUBBUILD_THREADS}")

# Prepare the `kmc` library — required by the Cuttlefish algorithm implementation.
# NOTE: do something more intelligent below than the -j4
#[[
set(KMC_VERSION "v3.2.2")
message("Build system will fetch and install KMC3")
ExternalProject_Add(prj_kmc
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    git clone --recurse-submodules https://github.com/refresh-bio/KMC.git --branch ${KMC_VERSION}
    PATCH_COMMAND       patch --strip 1 < ${CMAKE_SOURCE_DIR}/patches/kmc_v3.2.2.patch
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/KMC
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external/
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       make -j ${KMC_SUBBUILD_THREADS} bin/libkmc_core.a
    INSTALL_COMMAND     cp bin/libkmc_core.a ${EXT_LIB} &&
                        cp include/kmc_runner.h ${EXT_INCLUDE}
)

add_library(kmc STATIC IMPORTED)
set_target_properties(kmc PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libkmc_core.a)
add_dependencies(kmc prj_kmc)
]]


# Add / install ISA-L.


# Prepare the `parlay` library—a parallel algorithms programming toolkit.
message("Build system will fetch and install parlaylib")
ExternalProject_Add(prj_parlaylib
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    git clone https://github.com/cmuparlay/parlaylib.git
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/parlaylib
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf include/parlay ${EXT_INCLUDE}
)


# Prepare the `rapidgzip` library—gzip decompressor with multi-cores, and its dependencies.
set(RAPIDGZIP_VERSION "rapidgzip-v0.13.2")
message("Build system will fetch and install rapidgzip")
ExternalProject_Add(prj_rapidgzip
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    git clone --recurse-submodules git@github.com:mxmlnkn/rapidgzip.git --branch ${RAPIDGZIP_VERSION}
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/rapidgzip
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external/
    CONFIGURE_COMMAND   mkdir -p build && cd build &&
                        cmake -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                        -DWITH_RPMALLOC=ON -DWITH_ISAL=ON ..
    BUILD_COMMAND       cd build && make -j
    INSTALL_COMMAND     mkdir -p ${EXT_INCLUDE}/rapidgzip &&
                        cp -rf src/core src/huffman src/indexed_bzip2 src/rapidgzip ${EXT_INCLUDE}/rapidgzip/ &&
                        cp -rf src/external/isa-l/include ${EXT_INCLUDE}/rapidgzip/isa-l &&
                        cp -rf src/external/rpmalloc/rpmalloc ${EXT_INCLUDE}/rapidgzip/ &&
                        cp -rf build/src/libisal_inflate.a build/src/librpmalloc.a ${EXT_LIB}/
)

add_library(isal_inflate STATIC IMPORTED)
set_target_properties(isal_inflate PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libisal_inflate.a)
add_dependencies(isal_inflate prj_rapidgzip)

add_library(rpmalloc STATIC IMPORTED)
set_target_properties(rpmalloc PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/librpmalloc.a)
add_dependencies(rpmalloc prj_rapidgzip)


# Prepare the `boost` library. It provides: support for preprocessor
# metaprogramming and a scalable concurrent hashmap.
set(BOOST_MIN_VERSION "1.85.0")
find_package(Boost ${BOOST_MIN_VERSION})
if(BOOST_FOUND)
    message("Found boost (v${BOOST_VERSION}) in the system")
else()
    message("Build system will fetch and install required Boost libraries")
    ExternalProject_Add(prj_boost
        DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
        DOWNLOAD_COMMAND    curl -k -L https://github.com/boostorg/boost/releases/download/boost-${BOOST_MIN_VERSION}/boost-${BOOST_MIN_VERSION}-cmake.tar.gz -o boost-${BOOST_MIN_VERSION}-cmake.tar.gz &&
                            tar -xzf boost-${BOOST_MIN_VERSION}-cmake.tar.gz &&
                            rm boost-${BOOST_MIN_VERSION}-cmake.tar.gz

        SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/boost-${BOOST_MIN_VERSION}
        BUILD_IN_SOURCE     TRUE
        INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
        CONFIGURE_COMMAND   mkdir -p build && cd build &&
                            cmake -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                            -DBOOST_INCLUDE_LIBRARIES=preprocessor$<SEMICOLON>unordered
                            -DCMAKE_BUILD_TYPE=Release
                            -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/boost-${BOOST_MIN_VERSION}/build/ ..
        BUILD_COMMAND       cd build && cmake --build . --target install --config Release
        INSTALL_COMMAND     cp -rf build/include/boost ${EXT_INCLUDE}/
    )
endif()


# Prepare the `cereal` library. It supports (de)serialization of data
# structures.
set(CEREAL_MIN_VERSION "1.3.2")
message("Build system will fetch and install cereal")
ExternalProject_Add(prj_cereal
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/USCiLab/cereal/archive/refs/tags/v${CEREAL_MIN_VERSION}.tar.gz -o cereal-${CEREAL_MIN_VERSION}.tar.gz &&
                        tar -xzf cereal-${CEREAL_MIN_VERSION}.tar.gz &&
                        rm cereal-${CEREAL_MIN_VERSION}.tar.gz
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/cereal-${CEREAL_MIN_VERSION}
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf include/cereal ${EXT_INCLUDE}
)


# The `Debug` configuration optimizes the program for debugging and enables full debug information.
# The `Release` configuration enables most compiler optimizations for speed and defines `NDEBUG`
# (No Debug) which will remove all traces of the standard library assert calls.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


# Help `conda` build for OSX through circumventing some of its old SDK-based checks.
if(CONDA_BUILD)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D_LIBCPP_DISABLE_AVAILABILITY")
endif()


# Add subdirectory `src` to the build; CMake will open `src/CMakeLists.txt` for such.
add_subdirectory(src)
