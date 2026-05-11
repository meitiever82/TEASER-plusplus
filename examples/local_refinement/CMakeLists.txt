# Standalone build (after `make install` of TEASER++):
#   cd examples/local_refinement && mkdir build && cd build && cmake .. && make
#
# When built from the top-level (BUILD_EXAMPLES=ON + BUILD_TEASER_FPFH=ON), this CMakeLists is
# bypassed and examples/CMakeLists.txt wires the target directly against in-tree libs.

cmake_minimum_required(VERSION 3.10)
project(local_refinement)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)
find_package(teaserpp REQUIRED)
find_package(PCL 1.9 REQUIRED COMPONENTS common features kdtree)
find_package(OpenMP)

add_executable(local_refinement local_refinement.cc)
target_link_libraries(local_refinement
    Eigen3::Eigen
    teaserpp::teaser_registration
    teaserpp::teaser_features
    teaserpp::teaser_io
    ${PCL_LIBRARIES}
)
if (OpenMP_CXX_FOUND)
    target_link_libraries(local_refinement OpenMP::OpenMP_CXX)
endif()

# Copy bunny PLY next to the executable so demo mode works out of the box.
file(COPY ../example_data/
    DESTINATION ./example_data/
    FILES_MATCHING PATTERN *.ply)
