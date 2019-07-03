cmake_minimum_required (VERSION 3.8)
project(GFLinalg)

enable_testing()
add_library(GFLinalg INTERFACE)
target_include_directories(GFLinalg INTERFACE include)

add_subdirectory(tests)
