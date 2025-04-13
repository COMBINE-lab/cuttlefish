
# =============================================================================

# This controls the `-I` (include paths) switch of the compiler invocation.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Define the targets.

# Bundle the cfcore_static source files into a variable.
set(PROJECT_SRC
        # kmc_api/kmc_file.cpp
        # kmc_api/kmer_api.cpp
        # kmc_api/mmer.cpp
        # kmc_api/Virtual_Prefix_File.cpp
        xxHash/xxhash.c
        Build_Params.cpp
        Application.cpp
        Seq_Input.cpp
        Ref_Parser.cpp
        Async_Logger_Wrapper.cpp
        Thread_Pool.cpp
        DNA_Utility.cpp
        Kmer_Utility.cpp
        Vertex.cpp
        State.cpp
        Kmer_Container.cpp
        Kmer_SPSC_Iterator.cpp
        # Kmer_Hash_Table.cpp
        # CdBG.cpp
        # CdBG_Builder.cpp
        # CdBG_Writer.cpp
        # CdBG_Plain_Writer.cpp
        # CdBG_GFA_Writer.cpp
        # CdBG_GFA_Reduced_Writer.cpp
        # kmer_Enumerator.cpp
        # kmer_Enumeration_Stats.cpp
        # State_Read_Space.cpp
        # Read_CdBG.cpp
        # Read_CdBG_Constructor.cpp
        # Read_CdBG_Extractor.cpp
        Maximal_Unitig_Scratch.cpp
        Unipaths_Meta_info.cpp
        Data_Logistics.cpp
        dBG_Utilities.cpp
        Character_Buffer_Flusher.cpp
        Progress_Tracker.cpp
        dBG_Info.cpp
        # Validator.cpp
        # Validator_Hash_Table.cpp
        # Sequence_Validator.cpp
        # Kmers_Validator.cpp
        utility.cpp
        # Kmer_Index.cpp
        # Kmer_Index_Utility.cpp
        Minimizer_Instance_Iterator.cpp
        Multiway_Merger.cpp
        Discontinuity_Graph_Bootstrap.cpp
        dBG_Contractor.cpp
        Parser.cpp
        Graph_Partitioner.cpp
        Atlas.cpp
        Super_Kmer_Bucket.cpp
        Super_Kmer_Chunk.cpp
        HyperLogLog.cpp
        Discontinuity_Graph.cpp
        Edge_Matrix.cpp
        Unitig_File.cpp
        Subgraphs_Manager.cpp
        Subgraph.cpp
        State_Config.cpp
        Discontinuity_Graph_Contractor.cpp
        Contracted_Graph_Expander.cpp
        Unitig_Collator.cpp
        Unitig_Coord_Bucket.cpp
        Color_Table.cpp
        Color_Repo.cpp
        profile.cpp
        commands.cpp
    )

# Define the core library.
# `STATIC` will force the creation of a static library (and so no runtime linking).
add_library(cfcore_static STATIC ${PROJECT_SRC})

# Define the `RabbitFX` library.
aux_source_directory(${INCLUDE_DIR}/RabbitFX/io RABBITFX_SRC)
add_library(rabbitfx STATIC ${RABBITFX_SRC})

# Define required includes for `parlayhash`.
set(PARLAYHASH_INCLUDE ${INCLUDE_DIR}/parlayhash)
list(APPEND INCLUDE_DIR ${PARLAYHASH_INCLUDE})


# Define an executable from the driver program.
add_executable(${PROJECT_NAME} main.cpp)
add_executable(test test.cpp)


# Define dependencies between the targets.
add_dependencies(rabbitfx prj_rapidgzip)
add_dependencies(cfcore_static prj_parlaylib prj_boost)


# Set the include directories for the targets.

set(RAPIDGZIP_INCLUDE ${EXT_INCLUDE}/rapidgzip ${EXT_INCLUDE}/rapidgzip/core ${EXT_INCLUDE}/rapidgzip/indexed_bzip2 ${EXT_INCLUDE}/rapidgzip/rapidgzip
                      ${EXT_INCLUDE}/rapidgzip/isa-l/ ${EXT_INCLUDE}/rapidgzip/rpmalloc)

target_include_directories(cfcore_static PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(rabbitfx PUBLIC ${RAPIDGZIP_INCLUDE})
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_include_directories(test PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE} ${RAPIDGZIP_INCLUDE})


# Set the compile flags for the targets.

target_compile_options(cfcore_static PRIVATE ${COMPILE_FLAGS})
target_compile_options(rabbitfx PRIVATE ${COMPILE_FLAGS})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_FLAGS})
target_compile_options(test PRIVATE ${COMPILE_FLAGS})


# Link appropriate libraries and targets.

# Link to `jemalloc`, for better `malloc` support.
target_link_libraries(cfcore_static PRIVATE jemalloc)

# Link to `dl`, required in using dynamic shared object.
# Needed by `jemalloc`.
target_link_libraries(cfcore_static PRIVATE ${CMAKE_DL_LIBS})

# Link to `kmc`, for k-mer enumeration in cuttlefish 1/2.
# target_link_libraries(cfcore_static PRIVATE kmc)

# Link to `zlib`, required to support `gzip` files.
target_link_libraries(cfcore_static PRIVATE ZLIB::ZLIB)

# Link to `lz4`, for fast (de)compression of some intermedediate data.
target_link_libraries(cfcore_static INTERFACE ${LZ4_LIBRARIES})

# Link to the threads library of the platform.
target_link_libraries(cfcore_static PRIVATE Threads::Threads)

# Link to `rabbitfx`, for FASTX parsing.
target_link_libraries(cfcore_static PRIVATE rabbitfx)

# When building on conda on OSX with clang, we apparently need to explicitly
# link to the filesystem library; see (https://github.com/bioconda/bioconda-recipes/pull/37269/files).
if("$ENV{CF_INSTALL}" STREQUAL "conda" AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(cfcore_static PRIVATE c++experimental)
endif()

# Link to `zlib`, required to support `gzip` files.
target_link_libraries(rabbitfx PRIVATE ZLIB::ZLIB)

# Link to intel storage acceleration library's decompression core.
target_link_libraries(rabbitfx PRIVATE isal_inflate)

# Link to `rpmalloc` for `rapidgzip`'s memory allocation impl.
target_link_libraries(rabbitfx PRIVATE rpmalloc)

# Link the executable to the core library.
target_link_libraries(${PROJECT_NAME} PRIVATE cfcore_static)

# Link the tester.
target_link_libraries(test PRIVATE cfcore_static rabbitfx)


# Set the link flags for the targets.

target_link_libraries(cfcore_static PRIVATE ${LINK_FLAGS})
target_link_libraries(rabbitfx PRIVATE ${LINK_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LINK_FLAGS})
target_link_libraries(test PRIVATE ${LINK_FLAGS})


# Specialized options for memory-sanitizer.

if(MSAN)
    target_include_directories(rabbitfx PRIVATE ${LIBCXX_INCLUDE})
    target_include_directories(cfcore_static PRIVATE ${LIBCXX_INCLUDE})
    target_include_directories(${PROJECT_NAME} PRIVATE ${LIBCXX_INCLUDE})
    target_include_directories(test PRIVATE ${LIBCXX_INCLUDE})

    target_link_directories(rabbitfx PRIVATE ${LIBCXX})
    target_link_directories(cfcore_static PRIVATE ${LIBCXX})
    target_link_directories(test PRIVATE ${LIBCXX})

    add_dependencies(rabbitfx prj_libcxx)
    add_dependencies(cfcore_static prj_libcxx)
endif()


# Set the installation directory. The parent directory is `/usr/local/` by default, and the
# user may set it to a different path by passing a value for `CMAKE_INSTALL_PREFIX`. The target
# executable is installed inside a sub-directory named `bin` within the parent directory.
install(TARGETS ${PROJECT_NAME}
        DESTINATION bin
        )

install(TARGETS cfcore_static
  DESTINATION lib
)

# Apparently the cfcore_static_static library is not enough itself and tools
# wanting to use this downstream will need to link against `libkmc_core` as
# well.
#[[
install(FILES ${EXT_LIB}/libkmc_core.a
  DESTINATION lib
)
]]
