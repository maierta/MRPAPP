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
  LANGUAGES CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PROJECT_CMAKE ${mrpapp_SOURCE_DIR}/CMake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_CMAKE})

set(QMC_CXX_STANDARD
    17
    CACHE STRING "MRPAPP C++ language standard")
if(NOT QMC_CXX_STANDARD EQUAL 17)
  message(WARNING "C++17 is the only language standard officially supported by this MRPAPP version. "
                  "Using other versions of the C++ standard is unsupported and done entirely at user's own risk.")
endif()

find_package(MPI)

option(ENABLE_CUDA "Build with GPU support through CUDA" OFF)
option(ENABLE_HIP  "Build with GPU support through HIP" OFF)
set(ENABLE_GPU "$<IF:$<OR:$<BOOL:${ENABLE_CUDA}>,$<BOOL:${ENABLE_CUDA}>>,ON,OFF>")

if(ENABLE_GPU)
  if(ENABLE_CUDA)
    include(mrpapp_cuda)
  endif(ENABLE_CUDA)
  if(ENABLE_HIP)
    include(mrapp_hip)
  endif(ENABLE_HIP)

  if(MRPAPP_HAVE_CUDA OR MRPAPP_HAVE_HIP)
    include(DetermineDeviceArchitectures)
    message(STATUS "GPU device architectures: ${QMC_GPU_ARCHS}")
  endif()
endif(ENABLE_GPU)


