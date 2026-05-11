# Standalone build (after `make install` of TEASER++ with BUILD_TEASER_FPFH=ON):
#   cd examples/global_localization_sc && mkdir build && cd build && cmake .. && make
#
# When built from the top-level (BUILD_EXAMPLES=ON + BUILD_TEASER_FPFH=ON), this CMakeLists is
# bypassed and examples/CMakeLists.txt wires the target directly against in-tree libs.

cmake_minimum_required(VERSION 3.10)
project(global_localization_sc)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)
find_package(teaserpp REQUIRED)
find_package(PCL 1.9 REQUIRED COMPONENTS common features kdtree filters)
find_package(OpenMP)

add_executable(global_localization_sc
    global_localization_sc.cc
    third_party/scancontext/Scancontext.cpp
)
target_include_directories(global_localization_sc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(global_localization_sc
    Eigen3::Eigen
    teaserpp::teaser_registration
    teaserpp::teaser_features
    teaserpp::teaser_io
    ${PCL_LIBRARIES}
)
if (OpenMP_CXX_FOUND)
    target_link_libraries(global_localization_sc OpenMP::OpenMP_CXX)
endif()
