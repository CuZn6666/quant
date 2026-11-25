cmake_minimum_required(VERSION 3.14)

project(QuantFin LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()
find_package(Threads REQUIRED)

option(QUANT_BUILD_VISUALIZATION "Build the quant executable with Matplot++ visualization" ON)

# Core library: GBM simulation, option payoffs, Monte Carlo pricing
add_library(quant_core STATIC
    src/gbm.cpp
    src/double_call_put_payoff.cpp
    src/monte_carlo_pricing.cpp
    src/parallel_monte_carlo.cpp
)
target_link_libraries(quant_core PUBLIC Threads::Threads)
# Header files
target_include_directories(quant_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(quant_core PUBLIC cxx_std_11)

if(QUANT_BUILD_VISUALIZATION)
    find_package(Matplot++ REQUIRED)

    add_executable(quant
        src/main.cpp
        src/visualization.cpp
    )
    # Matplot++ needs C++17; the pricing library, tests and benchmarks use C++11.
    target_compile_features(quant PRIVATE cxx_std_17)
    target_link_libraries(quant PRIVATE
        quant_core
        Matplot++::matplot
    )
endif()

if(BUILD_TESTING)
    add_executable(test_pricing tests/test_pricing.cpp)
    target_link_libraries(test_pricing PRIVATE quant_core)

    foreach(test_case IN ITEMS
        asian_call
        asian_put
        lookback_call
        lookback_put
        gbm_properties
        monte_carlo_properties
        monte_carlo_zero_volatility
        monte_carlo_one_step
        invalid_inputs
        statistics_raw
        statistics_split_merge
        statistics_unequal_batches
        statistics_empty_singleton
        statistics_discounting
        single_threaded_regression
        statistics_count_overflow
    )
        add_test(NAME ${test_case} COMMAND test_pricing ${test_case})
        set_tests_properties(${test_case} PROPERTIES TIMEOUT 30)
    endforeach()

    add_executable(test_parallel_pricing tests/test_parallel_pricing.cpp)
    target_link_libraries(test_parallel_pricing PRIVATE quant_core)
    foreach(test_case IN ITEMS single_thread reproducibility statistical_consistency
            partition_merge zero_volatility small_samples invalid_parameters overflow)
        add_test(NAME parallel_${test_case} COMMAND test_parallel_pricing ${test_case})
        set_tests_properties(parallel_${test_case} PROPERTIES TIMEOUT 30)
    endforeach()
endif()

# -------- Single-Threaded Benchmark --------

add_executable(benchmark_single
    benchmarks/benchmark.cpp
)

target_link_libraries(benchmark_single PRIVATE
    quant_core
)
add_executable(
    benchmark_parallel
    benchmarks/benchmark_parallel.cpp
)

target_link_libraries(
    benchmark_parallel
    PRIVATE quant_core
)
target_compile_definitions(benchmark_parallel PRIVATE QUANT_BUILD_CONFIG="$<CONFIG>")
