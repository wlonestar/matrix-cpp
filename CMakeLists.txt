cmake_minimum_required(VERSION 3.21)
project(matrix)

set(CMAKE_CXX_STANDARD 20)

include_directories(${CMAKE_SOURCE_DIR})

add_executable(matrix main.cpp)
