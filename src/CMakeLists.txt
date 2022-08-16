cmake_minimum_required(VERSION 3.1.0)
project(matrix)

include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(matrix Matrix.cpp)