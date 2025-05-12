add_subdirectory(util)

set(srcs vector.cpp)
add_library(linalg STATIC ${srcs})
target_include_directories(linalg PUBLIC
"${PROJECT_SOURCE_DIR};${PROJECT_BINARY_DIR};${PROJECT_SOURCE_DIR}/linalg;${PROJECT_SOURCE_DIR}/linalg/util")


if(MRPAPP_HAVE_GPU)
   target_link_libraries(linalg PRIVATE gpu_utils)
endif(MRPAPP_HAVE_GPU)

get_property(linalg_include_dirs TARGET linalg PROPERTY INCLUDE_DIRECTORIES)
message("linalg target exports includes: ${linalg_include_dirs}")


add_subdirectory(tests)
