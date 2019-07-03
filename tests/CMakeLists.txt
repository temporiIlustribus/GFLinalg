if(RUN_BENCHMARK)
    find_package(benchmark REQUIRED)
    # target_link_libraries(test1 GFLinalg benchmark::benchmark)
endif()

if(RUN_TESTS)
    add_executable(test1 GFtest1.cpp)
    target_link_libraries(test1 GFLinalg)
    add_test(test1 test1)
endif()

