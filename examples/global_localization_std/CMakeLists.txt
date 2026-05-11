# Standalone build:
#   cd examples/global_localization_std && mkdir build && cd build && cmake .. && make
#
# Needs PCL + Ceres 2.1+ (for Manifold API). System apt's libceres-dev 2.0 won't compile;
# install Ceres 2.1+ to /usr/local (CMake will pick it up first).

cmake_minimum_required(VERSION 3.10)
project(global_localization_std)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
find_package(PCL 1.9 REQUIRED COMPONENTS common io kdtree)
find_package(Ceres REQUIRED)
find_package(OpenMP)

add_executable(std_eval_22submaps
    std_eval_22submaps.cc
    third_party/std/STDesc.cpp
)
target_include_directories(std_eval_22submaps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(std_eval_22submaps
    Eigen3::Eigen
    ${PCL_LIBRARIES}
    Ceres::ceres
)
if (OpenMP_CXX_FOUND)
    target_link_libraries(std_eval_22submaps OpenMP::OpenMP_CXX)
endif()
