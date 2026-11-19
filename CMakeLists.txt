cmake_minimum_required(VERSION 3.14)
#This is to meet the build environment requirements specified in the official Matplot++ documentation.

project(QuantFin LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Matplot++
find_package(Matplot++ REQUIRED)

# Create executable
add_executable(quant src/main.cpp)

# Link Matplot++
target_link_libraries(quant PRIVATE Matplot++::matplot)