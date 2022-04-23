cmake_minimum_required(VERSION 3.13)
project(murmur_hash)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


file(GLOB murmur_files "*.cpp")
add_library(murmur_lib ${murmur_files})
target_include_directories(murmur_lib INTERFACE .)
