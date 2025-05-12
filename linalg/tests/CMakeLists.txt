# linalg unit tests

add_executable(test_vector vector_cpu_gpu_test.cpp)
target_link_libraries(test_vector PUBLIC linalg Catch2::Catch2WithMain)
mrpapp_gpu_runtime_link(test_vector)




