cmake_minimum_required(VERSION 3.14)

project(QuantFin LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)

add_library(quant_core STATIC
    src/gbm.cpp
    src/payoff.cpp
    src/pricing.cpp
)

target_include_directories(quant_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(quant_core PUBLIC cxx_std_17)

add_executable(quant
    src/main.cpp
    src/visualization.cpp
)

target_link_libraries(quant PRIVATE quant_core Matplot++::matplot)
