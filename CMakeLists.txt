cmake_minimum_required(VERSION 3.6)

# ------------ general setup -----------
project(arap)

# add our cmake modules under cmake/
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Include CPM dependency manager
include(CPM)

# Set C++ flags
set(CMAKE_CXX_STANDARD 14)

CPMAddPackage(
        NAME eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 3.3.9
        DOWNLOAD_ONLY YES  # Eigen's CMakelists are not intended for library use
)
add_library(Eigen3::Eigen INTERFACE IMPORTED)
target_include_directories(Eigen3::Eigen INTERFACE ${eigen3_SOURCE_DIR})

CPMAddPackage(
        NAME libigl
        GIT_REPOSITORY http://github.com/libigl/libigl/
        GIT_TAG master
)

CPMAddPackage(
        NAME ceres
        GIT_REPOSITORY http://github.com/ceres-solver/ceres-solver/
        GIT_TAG 2.0.0
        DOWNLOAD_ONLY YES
)
add_library(ceres INTERFACE IMPORTED)
target_include_directories(ceres INTERFACE ${ceres_SOURCE_DIR})

# Set files to be compiled
set(HEADER_FILES
        arap/VertexSelectionPlugin.h
        arap/ARAP_Compute.h)

set(SOURCE_FILES
        arap/VertexSelectionPlugin.cpp
        arap/ARAP_Compute.cpp)

add_executable(arap arap/main.cpp ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(arap igl::core igl::opengl igl::opengl_glfw_imgui igl::opengl_glfw ceres)
