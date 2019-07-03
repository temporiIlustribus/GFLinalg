cmake_minimum_required (VERSION 3.8)
project(GFLinalg)

option(RUN_TESTS "compile and run unit tests" OFF)
option(RUN_BENCHMARK "compile and run the benchmarks. Requires google benchmark" OFF)

enable_testing()
add_library(GFLinalg INTERFACE)
target_include_directories(GFLinalg INTERFACE include)

add_subdirectory(tests)
