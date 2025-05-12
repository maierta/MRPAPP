# // Copyright (C) 2024 UT-Battelle, LLC
# // All rights reserved.
# //
# // See LICENSE for terms of usage.
# //

######################################################################
# CMake version and policies
######################################################################
cmake_minimum_required(VERSION 3.23.0)
# Note that cmake_minimum_required affects policy defaults.
# All policies known to the running version of CMake and introduced in
# cmake_minimum_required version or earlier will be set to use NEW behavior

######################################################################
# MRPAPP project
######################################################################
project(
  mrpapp
  VERSION 1.0.0
  LANGUAGES C CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PROJECT_CMAKE ${mrpapp_SOURCE_DIR}/CMake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_CMAKE})

set(MRPAPP_CXX_STANDARD
    17
    CACHE STRING "MRPAPP C++ language standard")
if(NOT MRPAPP_CXX_STANDARD EQUAL 17)
  message(WARNING "C++17 is the only language standard officially supported by this MRPAPP version. "
                  "Using other versions of the C++ standard is unsupported and done entirely at user's own risk.")
endif()

option(USE_MPI "use mpi" ON)
if(USE_MPI)
  find_package(MPI)
endif(USE_MPI)

option(ENABLE_CUDA "Build with GPU support through CUDA" OFF)
option(ENABLE_HIP  "Build with GPU support through HIP" OFF)
set(ENABLE_GPU "$<IF:$<OR:$<BOOL:${ENABLE_CUDA}>,$<BOOL:${ENABLE_CUDA}>>,ON,OFF>")
message ("enable gpu: ${ENABLE_GPU}")
set(MRPAPP_GPU_LIBS "" CACHE INTERNAL "")

if(ENABLE_GPU)
  if(ENABLE_CUDA)
    include(mrpapp_cuda)
  endif(ENABLE_CUDA)
  if(ENABLE_HIP)
    include(mrapp_hip)
  endif(ENABLE_HIP)
  if(MRPAPP_HAVE_CUDA OR MRPAPP_HAVE_HIP)
    include(DetermineDeviceArchitectures)
    message(STATUS "GPU device architectures: ${MRPAPP_GPU_ARCHS}")
  endif()
endif(ENABLE_GPU)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

include(ExternalProject)
find_package(Git REQUIRED)

find_package(Catch2 3)
if (NOT CATCH2_FOUND)
ExternalProject_Add(
  catch2
  PREFIX ${CMAKE_BINARY_DIR}/catch2
  GIT_REPOSITORY git@github.com:catchorg/Catch2.git
  GIT_TAG devel
  TIMEOUT 30
  UPDATE_COMMAND ${GIT_EXECUTABLE} pull
  CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/catch2"
  LOG_DOWNLOAD ON)

set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/catch2)
find_package(Catch2 3 REQUIRED)
endif()

set(MRPAPP_SRC
main.cpp)

set(MRPAPP_MODEL "1BAND" CACHE STRING "RPA model options")
set_property(CACHE MRPAPP_MODEL PROPERTY STRINGS SRRUO SRRUO3D SRRUO3DSUH 1BAND 1BANDWSPIN BILAYER_FESC BILAYER_1BAND ORTHOIIBILAYER BSCCOBILAYER BILAYER_FESC BAFEAS KFE2SE2 FOURORBITAL TBFILE COUPLEDLADDERS NDNIO2 MODELFROMFILESO KAGOME 1BANDABWSPIN 1BANDALTERMAGNET 1BANDAB)

include(mrpapp_defines)
mrpapp_write_definitions_file()

include(mrpapp_linking)
add_subdirectory(linalg)

add_executable(mrpapp ${MRPAPP_SRC})
target_include_directories(mrpapp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
  ${PROJECT_SOURCE_DIR}/PartialPsimag 
	${PROJECT_SOURCE_DIR}/PsimagLite 
	${PROJECT_SOURCE_DIR}/PsimagLite/src )

set(HAVE_MPI "$<BOOL:${MPI_FOUND}>")
target_compile_definitions(mrpapp PRIVATE "$<${HAVE_MPI}:USE_MPI>" USE_SCGAP3D "USE_${MRPAPP_MODEL}")

set(mpi_targs MPI::MPI_C MPI::MPI_CXX)
set(MPI_TARGETS "$<$<BOOL:${MPI_FOUND}>:${mpi_targs}>")
target_link_libraries(mrpapp ${MRPAPP_GPU_LIBS} LAPACK::LAPACK BLAS::BLAS ${MPI_TARGETS} )

add_custom_target(genexdebug COMMAND ${CMAKE_COMMAND} -E echo "$<IF:$<OR:$<BOOL:${ENABLE_CUDA}>,$<BOOL:${ENABLE_CUDA}>>,ON,OFF>")
#$<$<BOOL:${MPI_FOUND}>:USE_MPI>")

include(mrpapp_tests)
