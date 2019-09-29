if(RUN_BENCHMARK)
    find_package(benchmark REQUIRED)
	add_executable(benchmark1 GFbenchmark1.cc)
    target_link_libraries(benchmark1 GFLinalg benchmark::benchmark)
	add_test(benchmark1 benchmark1)
endif()

if(RUN_TESTS)
    add_executable(test1 GFtest1.cpp GFStorageTest.cpp)
    target_link_libraries(test1 GFLinalg)
    add_test(test1 test1)
endif()

