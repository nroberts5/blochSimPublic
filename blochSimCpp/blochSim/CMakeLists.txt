cmake_minimum_required(VERSION 3.7)
project(blochSim)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP)


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_library(blochSim STATIC blochSim.cpp)
target_link_libraries(blochSim ${ARMADILLO_LIBRARIES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(blochSim OpenMP::OpenMP_CXX)
endif()