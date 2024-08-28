################################################################################
# Author: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
#
# Checks for HIP and and accordingly sets MRPAPP_HAVE_HIP

set(ROCM_ROOT
  "/opt/rocm-4.2.0"
    CACHE PATH "Root directory of ROCM")
message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")
list(APPEND CMAKE_MODULE_PATH ${ROCM_ROOT}/hip/cmake)
list(APPEND CMAKE_PREFIX_PATH ${ROCM_ROOT})
find_package(HIP REQUIRED)
find_package(hipblas REQUIRED)
find_package(rocsolver REQUIRED)
# architecture flags
set(CMAKE_HIP_ARCHITECTURES
  "gfx906,gfx908"
  CACHE STRING "HIP architecture gfxXXX")
list(APPEND HIP_HIPCC_FLAGS "-fPIC -ffast-math -O3")
list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=${HIP_ARCH}")
list(APPEND HIP_HIPCC_FLAGS "--gpu-max-threads-per-block=256")
# warning suppression
list(APPEND HIP_HIPCC_FLAGS "-Wno-vla")
list(APPEND HIP_HIPCC_FLAGS "-Wno-deprecated-declarations")
list(APPEND HIP_HIPCC_FLAGS "-Wno-unused-command-line-argument")
list(APPEND HIP_HIPCC_FLAGS "-DHIP_PLATFORM_AMD")

# #-------------------------------------------------------------------
# #  set up ROCM compiler options and libraries
# #-------------------------------------------------------------------
if(MRPAPP_WITH_HIP)
  if(${CMAKE_VERSION} VERSION_LESS "3.21.3")
    message(FATAL_ERROR "Compilation for HIP requires CMake 3.21.3 or later.")
  endif()
  set(ENABLE_HIP 1)
  message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")

#-------------------------------------------------------------------
#  set up HIP compiler options
#-------------------------------------------------------------------
  set(CMAKE_MODULE_PATH "${ROCM_ROOT}/hip/cmake" ${CMAKE_MODULE_PATH})
  find_package(HIP REQUIRED)
  find_package(hipblas REQUIRED)
  find_package(hipsparse REQUIRED)
  find_package(rocsolver REQUIRED)

endif(MRPAPP_WITH_HIP)

#set(CUDA_ARCHITECTURES "sm_60" CACHE STRING "Name of the real architecture to build for.")
set(MAGMA_ROOT "" CACHE PATH "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(MRPAPP_HAVE_HIP FALSE CACHE INTERNAL "")
set(MRPAPP_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(MRPAPP_GPU_LIBS "" CACHE INTERNAL "")

include(CheckLanguage)
check_language(HIP)
if (CMAKE_HIP_COMPILER)
  enable_language(HIP)
  list(APPEND CMAKE_HIP_FLAGS "-fgpu-rdc")
  set(MRPAPP_HAVE_HIP TRUE CACHE INTERNAL "")
  set(MRPAPP_HAVE_GPU TRUE CACHE INTERNAL "")
  # Probably probably these should be public properties of the hip targets
  dca_add_haves_define(MRPAPP_HAVE_HIP)
  dca_add_haves_define(MRPAPP_HAVE_GPU)
  dca_add_haves_define(__HIP_PLATFORM_AMD__)
  list(APPEND MRPAPP_GPU_LIBS hip::host roc::hipblas)
  set(MRPAPP_HIP_PROPERTIES "CMAKE_HIP_ARCHITECTURES gfx906,gfx908")
  set(CMAKE_HIP_STANDARD 17)
  list(APPEND HIP_HIPCC_FLAGS "-fPIC")
  # doesn't appear to work
  set(CMAKE_HIP_SOURCE_FILE_EXTENSIONS cu)
  # NOTE: this is solved by dca_linking.cmake: dca_gpu_device_link()
  # alternative method (same issue)
  #file(GLOB_RECURSE CUDA_KERNELS_SRC ${PROJECT_SOURCE_DIR} *.cu)
  #set_source_files_properties(${CUDA_KERNELS_SRC} PROPERTIES LANGUAGE HIP)
# -ffast-math -O3")
# list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=${HIP_ARCH}")
# list(APPEND HIP_HIPCC_FLAGS "--gpu-max-threads-per-block=256")
# # warning suppression
# list(APPEND HIP_HIPCC_FLAGS "-Wno-vla")
# list(APPEND HIP_HIPCC_FLAGS "-Wno-deprecated-declarations")
# list(APPEND HIP_HIPCC_FLAGS "-Wno-unused-command-line-argument")
# list(APPEND HIP_HIPCC_FLAGS "-DHIP_PLATFORM_AMD")
endif()
