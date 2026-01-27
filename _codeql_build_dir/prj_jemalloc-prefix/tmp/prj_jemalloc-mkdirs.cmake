# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/runner/work/cuttlefish/cuttlefish/external/jemalloc-5.3.0")
  file(MAKE_DIRECTORY "/home/runner/work/cuttlefish/cuttlefish/external/jemalloc-5.3.0")
endif()
file(MAKE_DIRECTORY
  "/home/runner/work/cuttlefish/cuttlefish/_codeql_build_dir/prj_jemalloc-prefix/src/prj_jemalloc-build"
  "/home/runner/work/cuttlefish/cuttlefish/external"
  "/home/runner/work/cuttlefish/cuttlefish/_codeql_build_dir/prj_jemalloc-prefix/tmp"
  "/home/runner/work/cuttlefish/cuttlefish/_codeql_build_dir/prj_jemalloc-prefix/src/prj_jemalloc-stamp"
  "/home/runner/work/cuttlefish/cuttlefish/external"
  "/home/runner/work/cuttlefish/cuttlefish/_codeql_build_dir/prj_jemalloc-prefix/src/prj_jemalloc-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/runner/work/cuttlefish/cuttlefish/_codeql_build_dir/prj_jemalloc-prefix/src/prj_jemalloc-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/runner/work/cuttlefish/cuttlefish/_codeql_build_dir/prj_jemalloc-prefix/src/prj_jemalloc-stamp${cfgdir}") # cfgdir has leading slash
endif()
