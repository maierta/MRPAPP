# cuda and cublas utils

set(srcs allocators/aligned_allocator.cpp)
add_library(platform STATIC ${srcs})
target_include_directories(platform PUBLIC
"${PROJECT_SOURCE_DIR};${PROJECT_BINARY_DIR};${PROJECT_SOURCE_DIR}/platform;")

if(MRPAPP_HAVE_GPU)
  set(srcs util_gpublas.cpp error_gpuBLAS.cpp info_gpu.cpp error_gpu.cpp)
  add_library(gpu_utils STATIC ${srcs})
  message("project source dir in ${PROJECT_SOURCE_DIR}")
  target_include_directories(gpu_utils PUBLIC ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
  mrpapp_gpu_runtime_link(gpu_utils)
  mrpapp_gpu_blas_link(gpu_utils)
  get_property(util_include_dirs TARGET gpu_utils PROPERTY INCLUDE_DIRECTORIES)
  message("gpu_utils includes: ${util_include_dirs}")
  get_property(interface_util_include_dirs TARGET gpu_utils PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  message("gpu_utils interface includes: ${interface_util_include_dirs}")
endif()

