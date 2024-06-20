# // Copyright (C) 2023 UT-Battelle, LLC
# // All rights reserved.
# //
# // See LICENSE for terms of usage.
# //

# Checks for CUDA and accordingly sets MRPAPP_HAVE_CUDA
# In addition, set MRPAPP_GPU_LIBS.
set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING "Name of the real architecture to build for.")

set(MRPAPP_HAVE_CUDA FALSE CACHE INTERNAL "")
set(MRPAPP_GPU_LIBS "" CACHE INTERNAL "")

# Find CUDA.
include(CheckLanguage)

find_package(CUDAToolkit REQUIRED)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  set(MRPAPP_HAVE_CUDA TRUE CACHE INTERNAL "")
  set(MRPAPP_HAVE_GPU TRUE CACHE INTERNAL "")
  dca_add_haves_define(MRPAPP_HAVE_CUDA)
  dca_add_haves_define(MRPAPP_HAVE_GPU)

  list(APPEND MRPAPP_GPU_LIBS CUDA::cudart CUDA::cublas)
  set(MRPAPP_CUDA_PROPERTIES "CMAKE_CUDA_ARCHITECTURES 70")
  list(APPEND CUDAFLAGS "--expt-relaxed-constexpr" ${MRPAPP_CUDA_OPTIONS})
  set(CMAKE_CUDA_STANDARD 17)
endif()
