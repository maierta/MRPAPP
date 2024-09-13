# // Copyright (C) 2024 UT-Battelle, LLC
# // All rights reserved.
# //
# // See LICENSE for terms of usage.
# //

# Checks for CUDA and accordingly sets MRPAPP_HAVE_CUDA
# In addition, set MRPAPP_GPU_LIBS.
message("checking CUDA environment")
set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING "Name of the real architecture to build for.")

set(MRPAPP_HAVE_CUDA FALSE CACHE INTERNAL "")
set(MRPAPP_GPU_LIBS "" CACHE INTERNAL "")

include(mrpapp_defines)
include(CheckLanguage)

if(NOT CMAKE_CUDA_FLAGS MATCHES "allow-unsupported-compiler")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --allow-unsupported-compiler")
endif()

set(CMAKE_CUDA_HOST_COMPILER
  ${CMAKE_CXX_COMPILER}
  CACHE STRING "nvcc host compiler passed via -ccbin")

# Find CUDA.
find_package(CUDAToolkit REQUIRED)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
  message("Found CUDA compiler!")
  enable_language(CUDA)
  set(MRPAPP_HAVE_CUDA TRUE CACHE INTERNAL "")
  set(MRPAPP_HAVE_GPU TRUE CACHE INTERNAL "")
  mrpapp_add_define(MRPAPP_HAVE_CUDA)
  mrpapp_add_define(MRPAPP_HAVE_GPU)

  list(APPEND MRPAPP_GPU_LIBS CUDA::cudart CUDA::cublas)
  set(MRPAPP_CUDA_PROPERTIES "CMAKE_CUDA_ARCHITECTURES 70")
  list(APPEND CUDAFLAGS "--expt-relaxed-constexpr" ${MRPAPP_CUDA_OPTIONS})
  set(CMAKE_CUDA_STANDARD 17)
endif()
