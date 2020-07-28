
# Bundle the core source files into a variable.
set(PROJECT_SRC
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        Kmer.cpp
        Vertex.cpp
        State.cpp
        Kmer_Container.cpp
        Kmer_Hash_Table.cpp
        CdBG.cpp
        CdBG_Builder.cpp
        CdBG_Writer.cpp
        CdBG_GFA_Writer.cpp
        Validator.cpp
        utility.cpp
    )


# Build a library with the provided source(s).
# `STATIC` will force the creation of a static library (and so no runtime linking).
add_library(core STATIC ${PROJECT_SRC})


# This controls the `-I` (include paths) switch of the compiler invocation.
# Since `include` is declared `PUBLIC` for the library, it is propagated to all the
# consumers of the library.
set(INCLUDE_DIR include)
target_include_directories(core PUBLIC ${CMAKE_SOURCE_DIR}/${INCLUDE_DIR})


# Specify warnings flags to the compiler for the target library. As the warnings are specified as
# `PRIVATE`, they will only be used when compiling this library. The flags will not be propagated to
# the targets linking to it.
target_compile_options(core PRIVATE ${WARNING_FLAGS})


# Create an executable from the driver program.
add_executable(${PROJECT_NAME} main.cpp)


# Link the core library to the threads package in the platform.
target_link_libraries(core PRIVATE ${CMAKE_THREAD_LIBS_INIT})

# Link the core library to the `std::experimental::filesystem` library.
target_link_libraries(core PRIVATE stdc++fs)


# Link the executable to the required libraries.
target_link_libraries(${PROJECT_NAME} PRIVATE core)


# Specify warnings flags to the compiler for the target executable.
target_compile_options(${PROJECT_NAME} PRIVATE ${WARNING_FLAGS})


# Set the parent directory for the target binaries as `bin` directory inside the "installation directory".
# The installation directory is set to `/usr/local` by default. Installing softwares there ensures that
# those are automatically available to users. You may need root access (use `sudo`) for such.
# To specify a different parent installation directory, pass the directory as the value of the
# `-DCMAKE_INSTALL_PREFIX` with the `cmake` command.

# TODO: Employ a better design choice regarding `INSTALL_DIR` (by relating it to `CMAKE_INSTALL_PREFIX`)
set(INSTALL_DIR ${CMAKE_SOURCE_DIR}/bin)
install(TARGETS ${PROJECT_NAME}
        DESTINATION ${INSTALL_DIR}
        )


add_executable(test test.cpp Validator.cpp)
target_link_libraries(test PRIVATE core)
target_compile_options(test PRIVATE ${WARNING_FLAGS})
