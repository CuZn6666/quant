cmake_minimum_required(VERSION 3.14)

project(QuantFin LANGUAGES CXX)

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