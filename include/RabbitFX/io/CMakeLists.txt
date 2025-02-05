cmake_minimum_required(VERSION 3.0)

#openmp
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(EXECUTABLE_OUTPUT_PATH .)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -ffast-math ${CMAKE_CXX_FLAGS}")

aux_source_directory(. SOURCE_LIST)
add_library(rabbitfx ${SOURCE_LIST})
#message(${SOURCE_LIST})
#include_directories(.)

#add_executable (RabbitIO ${SOURCE_LIST})

target_link_libraries(rabbitfx z)
#target_link_libraries(RabbitIO)
