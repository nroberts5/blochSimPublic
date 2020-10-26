cmake_minimum_required(VERSION 3.7)
project(blochSimDemo)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP)

add_subdirectory(blochSim)
include_directories(blochSim)

add_executable(demo main.cpp params.h)
target_link_libraries(demo blochSim)
if(OpenMP_CXX_FOUND)
    target_link_libraries(demo OpenMP::OpenMP_CXX)
endif()