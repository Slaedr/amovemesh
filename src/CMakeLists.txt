cmake_minimum_required(VERSION 2.8)

project (amocurve)

# to compile with OpenMP
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# set compile options
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

# set search path(s) for pastix library here - normally no need of this
#find_library(PASTIX_LIB NAMES pastix PATHS /usr/lib64 /usr/lib/)

# PastiX needs a single-threaded BLAS, so set search path(s) for that here
#find_library(BLAS_LIB NAMES blas PATHS /usr/lib64 /usr/lib/)

# ---------------------------------------------------------------------------- #

# flags and stuff

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mavx")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Debug mode")
else()
	message(STATUS "Release mode")
endif()

if (DEFINED ENV{EIGEN_DIR})
	message (STATUS "Will build with Eigen from $ENV{EIGEN_DIR}")
	add_definitions (-DEIGEN_LIBRARY=1)
	include_directories ($ENV{EIGEN_DIR})
endif()

if (DEFINED PASTIX_LIB AND DEFINED BLAS_LIB)
	message(STATUS "Found PastiX at ${PASTIX_LIB} and BLAS at ${BLAS_LIB}")
	add_definitions(-DPASTIX_LIBRARY=1)
else()
	if(NOT DEFINED BLAS_LIB)
		message(STATUS "BLAS library not found!")
	endif()
	if(NOT DEFINED PASTIX_LIB)
		message(STATUS "PastiX library not found! ")
	endif()
	message(STATUS "Compiling without PastiX solver.")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# libraries to be compiled

add_library(amesh3d amesh3d.cpp)
target_link_libraries(amesh3d amatrix adatastructures)

add_library(alinalg alinalg.cpp asparsematrix.cpp)
target_link_libraries(alinalg amatrix)
if(DEFINED PASTIX_LIB)
	target_link_libraries(alinalg ${PASTIX_LIB} ${BLAS_LIB})
endif()

add_library(ageometry ageometry.cpp)
target_link_libraries(ageometry alinalg amatrix adatastructures)

add_library(ageometry3d ageometry3d.cpp)
target_link_libraries(ageometry3d alinalg amesh3d amatrix adatastructures)

add_library(arbf arbf.cpp)
target_link_libraries(arbf alinalg amatrix)

# add_library(aboundaryinfluence aboundaryinfluencedistance.cpp)
# target_link_libraries(aboundaryinfluence amesh2dh)

add_library(ageometryh ageometryh.cpp)
target_link_libraries(ageometryh alinalg amatrix)

add_library(aoutput aoutput.cpp)
target_link_libraries(aoutput amesh2dh)

add_library(amesh2dh amesh2dh.cpp)
target_link_libraries(amesh2dh adatastructures)

add_library(abowyerwatson3d abowyerwatson3d.cpp)
target_link_libraries(abowyerwatson3d amatrix adatastructures)

add_library(amatrix amatrix.cpp)

add_library(adatastructures adatastructures.cpp)

add_subdirectory(binsort)

# for the final executable(s)
add_subdirectory(DG)
add_subdirectory(curved-mesh-gen-splines)
add_subdirectory(curved-mesh-gen-cad)
add_subdirectory(bouncurve)
add_subdirectory(utilities)
