cmake_minimum_required(VERSION 3.1.0)
project(test)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)

add_executable(test test.cpp)
target_link_libraries(test matrix gtest)