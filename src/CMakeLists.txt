
# Bundle the core source files into a variable.
set(PROJECT_SRC
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        kmc_api/Virtual_Prefix_File.cpp
        xxHash/xxhash.c
        xxHash/xxhsum.c
        Build_Params.cpp
        Application.cpp
        Seq_Input.cpp
        Ref_Parser.cpp
        Async_Logger_Wrapper.cpp
        Thread_Pool.cpp
        DNA_Utility.cpp
        Kmer_Utility.cpp
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
        kmer_Enumerator.cpp
        kmer_Enumeration_Stats.cpp
        State_Read_Space.cpp
        Read_CdBG.cpp
        Read_CdBG_Constructor.cpp
        Read_CdBG_Extractor.cpp
        Unitig_Scratch.cpp
        Maximal_Unitig_Scratch.cpp
        Unipaths_Meta_info.cpp
        Data_Logistics.cpp
        dBG_Utilities.cpp
        Character_Buffer_Flusher.cpp
        Progress_Tracker.cpp
        dBG_Info.cpp
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
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
target_include_directories(core PUBLIC ${INCLUDE_DIR} ${EXT_INCLUDE})


# Specify the warnings and the extra optimization flags to the compiler for the target library.
# As the flags are specified as `PRIVATE`, they will only be used when compiling this library.
# These will not be propagated to the targets linking to it.
target_compile_options(core PRIVATE ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})


# Create an executable from the driver program.
add_executable(${PROJECT_NAME} main.cpp)


# Link the core library to the `jemalloc` library, for better `malloc` support.
target_link_libraries(core PRIVATE jemalloc)

# Link the core library to the `dl` library, required in using dynamic shared object.
# Needed by `jemalloc`.
target_link_libraries(core PRIVATE ${CMAKE_DL_LIBS})

# Link the core library to the `kmc` library.
target_link_libraries(core PRIVATE kmc)

# Link the core library to the `zlib` library.
target_link_libraries(core PRIVATE ZLIB::ZLIB)

# Link the core library ti the `bzip2` library.
target_link_libraries(core PRIVATE BZip2::BZip2)

# Link the core library to the threads package in the platform.
target_link_libraries(core PRIVATE ${CMAKE_THREAD_LIBS_INIT})


# Link the executable to the required libraries.
target_link_libraries(${PROJECT_NAME} PRIVATE core)


# Specify warnings flags to the compiler for the target executable.
target_compile_options(${PROJECT_NAME} PRIVATE ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})


# Set the installation directory. The parent directory is `/usr/local/` by default, and the
# user may set it to a different path by passing a value for `CMAKE_INSTALL_PREFIX`. The target
# executable is installed inside a sub-directory named `bin` within the parent directory.
install(TARGETS ${PROJECT_NAME}
        DESTINATION bin
        )


add_executable(test test.cpp)
target_link_libraries(test PRIVATE core)
target_compile_options(test PRIVATE ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})
