
# Bundle the core source files into a variable.
set(PROJECT_SRC
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        xxHash/xxhash.c
        xxHash/xxhsum.c
        Application.cpp
        Parser.cpp
        Async_Logger_Wrapper.cpp
        Thread_Pool.cpp
        DNA_Utility.cpp
        Kmer_u64.cpp
        Vertex.cpp
        State.cpp
        Kmer_Container.cpp
        Kmer_Hash_Table.cpp
        CdBG.cpp
        CdBG_Builder.cpp
        CdBG_Writer.cpp
        CdBG_Plain_Writer.cpp
        CdBG_GFA_Writer.cpp
        CdBG_GFA_Reduced_Writer.cpp
        Read_CdBG.cpp
        Read_CdBG_Constructor.cpp
        Read_CdBG_Extractor.cpp
        Validator.cpp
        Validator_Hash_Table.cpp
        Sequence_Validator.cpp
        Kmers_Validator.cpp
        utility.cpp
        Minimizer_Policy.cpp
    )


# Build a library with the provided source(s).
# `STATIC` will force the creation of a static library (and so no runtime linking).
add_library(core STATIC ${PROJECT_SRC})


# This controls the `-I` (include paths) switch of the compiler invocation.
# Since `include` is declared `PUBLIC` for the library, it is propagated to all the
# consumers of the library.
set(INCLUDE_DIR include)
target_include_directories(core PUBLIC ${CMAKE_SOURCE_DIR}/${INCLUDE_DIR})


# Specify the warnings and the extra optimization flags to the compiler for the target library.
# As the flags are specified as `PRIVATE`, they will only be used when compiling this library.
# These will not be propagated to the targets linking to it.
target_compile_options(core PRIVATE ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})


# Create an executable from the driver program.
add_executable(${PROJECT_NAME} main.cpp)


# Link the core library to the `jemalloc` library, for better `malloc` support.
target_link_libraries(core PRIVATE ${FAST_MALLOC_LIB})

# Link the core library to the `dl` library, required in using dynamic shared object.
# Needed by `jemalloc`.
target_link_libraries(core PRIVATE ${CMAKE_DL_LIBS})

# Link the core library to the threads package in the platform.
target_link_libraries(core PRIVATE ${CMAKE_THREAD_LIBS_INIT})

# Link the core library to the `zlib` library.
target_link_libraries(core PRIVATE ZLIB::ZLIB)


# Link the executable to the required libraries.
target_link_libraries(${PROJECT_NAME} PRIVATE core)


# Specify warnings flags to the compiler for the target executable.
target_compile_options(${PROJECT_NAME} PRIVATE ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})


# Set the parent directory for the target binaries inside the "installation directory". The
# installation directory is set to `/usr/local` by default. If it's set there, i.e. the user
# have not passed any other value for such, then the parent installation directory is set to
# a default (`bin`) directory at the top level of the source tree. To specify a different parent
# installation directory, pass its path as the value of the `-DCMAKE_INSTALL_PREFIX` with the
# `cmake` command.

set(DEFAULT_INSTALL_DIR ${CMAKE_SOURCE_DIR}/bin)        # Set our default installation directory as a `bin` directory inside the top level of the source tree.
set(INSTALL_DIR .)      # Set the installation directory to what CMake has recieved, i.e. either the system default or the user-provided one.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) # If it's set to the system default, then set it to our defined default directory.
        set(INSTALL_DIR ${DEFAULT_INSTALL_DIR})
endif()

install(TARGETS ${PROJECT_NAME}
        DESTINATION ${INSTALL_DIR}
        )


add_executable(test test.cpp)
target_link_libraries(test PRIVATE core)
target_compile_options(test PRIVATE ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})
