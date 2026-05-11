# Build examples against the in-tree teaserpp targets.
# The per-example subdirectories (teaser_cpp_ply/, teaser_cpp_fpfh/) keep their
# own standalone CMakeLists.txt that use find_package(teaserpp) — those are
# meant to be built after `make install`, not from this top-level build.
#
# Each executable is placed in its own subdirectory under the build tree so the
# layout matches the standalone build (./example_data/ next to the binary), and
# so the target name can't collide with the source subdirectory of the same
# name.

find_package(OpenMP)

set(_example_data_dir "${CMAKE_CURRENT_SOURCE_DIR}/example_data")

# teaser_cpp_ply
add_executable(teaser_cpp_ply teaser_cpp_ply/teaser_cpp_ply.cc)
set_target_properties(teaser_cpp_ply PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/teaser_cpp_ply")
target_link_libraries(teaser_cpp_ply
    PRIVATE
        Eigen3::Eigen
        teaserpp::teaser_registration
        teaserpp::teaser_io
)
if (OpenMP_CXX_FOUND)
    target_link_libraries(teaser_cpp_ply PRIVATE OpenMP::OpenMP_CXX)
endif()
file(COPY ${_example_data_dir}/
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/teaser_cpp_ply/example_data/
    FILES_MATCHING PATTERN *.ply)

# FPFH-based examples need PCL and the optional teaser_features library.
if (BUILD_TEASER_FPFH)
    add_executable(teaser_cpp_fpfh teaser_cpp_fpfh/teaser_cpp_fpfh.cc)
    set_target_properties(teaser_cpp_fpfh PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/teaser_cpp_fpfh")
    target_link_libraries(teaser_cpp_fpfh
        PRIVATE
            Eigen3::Eigen
            teaserpp::teaser_registration
            teaserpp::teaser_features
            teaserpp::teaser_io
            ${PCL_LIBRARIES}
    )

    add_executable(quatro_cpp_fpfh teaser_cpp_fpfh/quatro_cpp_fpfh.cc)
    set_target_properties(quatro_cpp_fpfh PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/teaser_cpp_fpfh")
    target_link_libraries(quatro_cpp_fpfh
        PRIVATE
            Eigen3::Eigen
            teaserpp::teaser_registration
            teaserpp::teaser_features
            teaserpp::teaser_io
            ${PCL_LIBRARIES}
    )

    if (OpenMP_CXX_FOUND)
        target_link_libraries(teaser_cpp_fpfh PRIVATE OpenMP::OpenMP_CXX)
        target_link_libraries(quatro_cpp_fpfh PRIVATE OpenMP::OpenMP_CXX)
    endif()

    file(COPY ${_example_data_dir}/
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/teaser_cpp_fpfh/example_data/
        FILES_MATCHING PATTERN *.ply)

    # local_refinement: rough-pose -> precise pose using FPFH + TEASER/Quatro.
    add_executable(local_refinement local_refinement/local_refinement.cc)
    set_target_properties(local_refinement PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/local_refinement")
    target_link_libraries(local_refinement
        PRIVATE
            Eigen3::Eigen
            teaserpp::teaser_registration
            teaserpp::teaser_features
            teaserpp::teaser_io
            ${PCL_LIBRARIES}
    )
    if (OpenMP_CXX_FOUND)
        target_link_libraries(local_refinement PRIVATE OpenMP::OpenMP_CXX)
    endif()
    file(COPY ${_example_data_dir}/
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/local_refinement/example_data/
        FILES_MATCHING PATTERN *.ply)

    # global_localization_sc: Scan Context retrieval + per-candidate TEASER refinement.
    add_executable(global_localization_sc
        global_localization_sc/global_localization_sc.cc
        global_localization_sc/third_party/scancontext/Scancontext.cpp
    )
    set_target_properties(global_localization_sc PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/global_localization_sc")
    target_include_directories(global_localization_sc
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/global_localization_sc)
    target_link_libraries(global_localization_sc
        PRIVATE
            Eigen3::Eigen
            teaserpp::teaser_registration
            teaserpp::teaser_features
            teaserpp::teaser_io
            ${PCL_LIBRARIES}
    )
    if (OpenMP_CXX_FOUND)
        target_link_libraries(global_localization_sc PRIVATE OpenMP::OpenMP_CXX)
    endif()
    file(COPY ${_example_data_dir}/
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/global_localization_sc/example_data/
        FILES_MATCHING PATTERN *.ply)
else()
    message(STATUS "BUILD_TEASER_FPFH=OFF: skipping teaser_cpp_fpfh / quatro_cpp_fpfh / local_refinement / global_localization_sc examples.")
endif()

# global_localization_std: STD-based retrieval evaluator (Phase 3 SC vs STD test).
# Gated on Ceres availability rather than BUILD_TEASER_FPFH (it doesn't need FPFH).
# Needs Ceres 2.1+ for the Manifold API; system apt's 2.0 won't work.
find_package(Ceres QUIET)
if (Ceres_FOUND)
    find_package(PCL QUIET COMPONENTS common io kdtree)
    if (PCL_FOUND)
        add_executable(std_eval_22submaps
            global_localization_std/std_eval_22submaps.cc
            global_localization_std/third_party/std/STDesc.cpp
        )
        set_target_properties(std_eval_22submaps PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/global_localization_std")
        set_property(TARGET std_eval_22submaps PROPERTY CXX_STANDARD 17)
        target_include_directories(std_eval_22submaps
            PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/global_localization_std)
        target_link_libraries(std_eval_22submaps
            PRIVATE Eigen3::Eigen ${PCL_LIBRARIES} Ceres::ceres)
        if (OpenMP_CXX_FOUND)
            target_link_libraries(std_eval_22submaps PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endif()
else()
    message(STATUS "Ceres not found: skipping global_localization_std (need Ceres 2.1+ for STD).")
endif()

# global_localization_bbs: 3D-BBS coarse pose + TEASER refinement (Phase 5).
# 3D-BBS is vendored in-tree at global_localization_bbs/third_party (KOKIAOKI/3d_bbs upstream).
# We add_subdirectory() it directly so cpu_bbs3d / gpu_bbs3d become available as targets.
# EXCLUDE_FROM_ALL: only built when bbs_smoke_test / bbs_localize is built.
set(_bbs3d_root "${CMAKE_CURRENT_SOURCE_DIR}/global_localization_bbs/third_party")
if (EXISTS "${_bbs3d_root}/CMakeLists.txt")
    find_package(CUDAToolkit QUIET)
    find_package(Boost QUIET COMPONENTS filesystem system)
    if (CUDAToolkit_FOUND AND Boost_FOUND)
        add_subdirectory(${_bbs3d_root}
                         ${CMAKE_CURRENT_BINARY_DIR}/global_localization_bbs/_bbs3d_build
                         EXCLUDE_FROM_ALL)
        set(_bbs3d_includes "${_bbs3d_root}/bbs3d/include")

        add_executable(bbs_smoke_test
            global_localization_bbs/bbs_smoke_test.cc
        )
        set_target_properties(bbs_smoke_test PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/global_localization_bbs")
        target_include_directories(bbs_smoke_test PRIVATE ${_bbs3d_includes})
        target_link_libraries(bbs_smoke_test
            PRIVATE
                Eigen3::Eigen
                gpu_bbs3d
                CUDA::cudart
                Boost::filesystem
                Boost::system
                teaserpp::teaser_io
        )
        file(COPY ${_example_data_dir}/
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/global_localization_bbs/example_data/
            FILES_MATCHING PATTERN *.ply)

        # bbs_localize: BBS coarse + TEASER refine end-to-end. Requires FPFH (PCL + teaser_features).
        if (BUILD_TEASER_FPFH)
            add_executable(bbs_localize
                global_localization_bbs/bbs_localize.cc
            )
            set_target_properties(bbs_localize PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/global_localization_bbs")
            set_property(TARGET bbs_localize PROPERTY CXX_STANDARD 17)
            target_include_directories(bbs_localize
                PRIVATE
                    ${_bbs3d_includes}
                    ${PCL_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}/global_localization_bbs)
            target_link_libraries(bbs_localize
                PRIVATE
                    Eigen3::Eigen
                    gpu_bbs3d
                    CUDA::cudart
                    Boost::filesystem
                    Boost::system
                    teaserpp::teaser_registration
                    teaserpp::teaser_features
                    teaserpp::teaser_io
                    ${PCL_LIBRARIES}
            )
            if (OpenMP_CXX_FOUND)
                target_link_libraries(bbs_localize PRIVATE OpenMP::OpenMP_CXX)
            endif()
        endif()
    else()
        message(STATUS "CUDAToolkit or Boost not found: skipping global_localization_bbs.")
    endif()
else()
    message(STATUS "global_localization_bbs/third_party/ not present: skipping Phase 5 build.")
endif()
