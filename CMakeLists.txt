cmake_minimum_required(VERSION 3.8)
project(aig CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-g -O3 -Wall)

add_library(aig ${CMAKE_CURRENT_SOURCE_DIR}/aig.cpp)
target_include_directories(aig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/aig.hpp)
