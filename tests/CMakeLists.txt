find_package(benchmark REQUIRED)

add_executable(test1 GFtest1.cpp)
target_link_libraries(test1 GFLinalg benchmark::benchmark)
add_test(test1 test1)

