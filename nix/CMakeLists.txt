cmake_minimum_required(VERSION 3.6)

# ------------ general setup -----------
project(arap)

# add our cmake modules under cmake/
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake" ${CMAKE_MODULE_PATH})

# Include CPM dependency manager
include(CPM)


# Options
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../libs CACHE PATH "Path to lib folder")
set(Eigen3_INCLUDE_DIR  $ENV{EIGEN_DIR}/include/eigen3 CACHE PATH "Path to Eigen source folder")
set(glog_INCLUDE_DIR ${GLOG_DIR}/include/glog-lib/include CACHE PATH "Path to glog source folder")
set(Ceres_DIR $ENV{CERES_DIR}/lib/Ceres/CMake CACHE PATH "Path to Ceres build folder")
set(Ceres_INCLUDE_DIR $ENV{CERES_DIR}/include)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${Ceres_DIR})


# Set C++ flags
set(CMAKE_CXX_STANDARD 14)


CPMAddPackage(
        NAME libigl
        GIT_REPOSITORY http://github.com/libigl/libigl/
        GIT_TAG master
)

find_package(Ceres REQUIRED)
find_package(glog REQUIRED ${glog_DIR_HINT})


#add_library(ceres INTERFACE IMPORTED)
target_include_directories(ceres INTERFACE ${ceres_SOURCE_DIR})


add_executable(arap 
    ../arap/main.cpp
    ../arap/VertexSelectionPlugin.cpp
    ../arap/VertexSelectionPlugin.h
)
target_link_libraries(arap igl::core igl::opengl_glfw ceres)

target_include_directories(arap PUBLIC ${Eigen3_INCLUDE_DIR} ${Ceres_INCLUDE_DIR} ${glog_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
