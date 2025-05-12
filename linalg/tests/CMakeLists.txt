# linalg unit tests

set(test_srcs vector_cpu_test.cpp)
if(MRPAPP_HAVE_GPU)
set(test_srcs ${test_srcs} vector_gpu_test.cpp vector_cpu_gpu_test.cpp)
endif()

add_executable(test_vector ${test_srcs})
target_link_libraries(test_vector PUBLIC platform linalg Catch2::Catch2WithMain)
mrpapp_gpu_runtime_link(test_vector)




