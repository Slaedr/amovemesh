cmake_minimum_required(VERSION 2.8)

project (amocurve)

# set compiler
set (CMAKE_C_COMPILER $ENV{GC})
set (CMAKE_CXX_COMPILER $ENV{GPP})

# un-comment to compile with OpenMP
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# un-comment to compile with Eigen 3
set (EIGEN_DIR ../eigen-3.2.8) 

# set compile options
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ggdb")


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -DDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Debug mode")
else()
	message(STATUS "Release mode")
endif()

if (DEFINED EIGEN_DIR)
	message (STATUS "Will build with Eigen from ${EIGEN_DIR}")
	add_definitions (-DEIGEN_LIBRARY=1)
	include_directories (${EIGEN_DIR})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# libraries to be compiled

add_library(amesh3d amesh3d.cpp)
target_link_libraries(amesh3d amatrix2 adatastructures)

add_library(alinalg alinalg.cpp)
target_link_libraries(alinalg amatrix2)

add_library(ageometry3d ageometry3d.cpp)
target_link_libraries(ageometry3d alinalg amesh3d amatrix2 adatastructures)

add_library(arbf arbf.cpp)
target_link_libraries(arbf alinalg amatrix2)

add_library(ageometryh ageometryh.cpp)
target_link_libraries(ageometryh alinalg amatrix2)

add_library(amesh2dh amesh2dh.cpp)
target_link_libraries(amesh2dh amatrix2 adatastructures)

add_library(amatrix2 amatrix2.cpp)

add_library(adatastructures adatastructures.cpp)

# for the final executable(s)
add_subdirectory(curved-mesh-gen-splines)
add_subdirectory(utilities)