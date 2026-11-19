cmake_minimum_required(VERSION 3.14)

project(QuantFin LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Matplot++ plotting library
find_package(Matplot++ REQUIRED)

# Core library: GBM simulation, option payoffs, Monte Carlo pricing
add_library(quant_core STATIC
    src/gbm.cpp
    src/double_call_put_payoff.cpp
    src/monte_carlo_pricing.cpp
)

# Header files
target_include_directories(quant_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(quant_core PUBLIC cxx_std_17)

# Main executable
add_executable(quant
    src/main.cpp
    src/visualization.cpp
)

# Link libraries
target_link_libraries(quant PRIVATE
    quant_core
    Matplot++::matplot
)

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
    )
        add_test(NAME ${test_case} COMMAND test_pricing ${test_case})
        set_tests_properties(${test_case} PROPERTIES TIMEOUT 30)
    endforeach()
endif()
