// Copyright (C) 2025 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests the interaction between Vector<CPU> and Vector<GPU>.

#include "linalg/vector.hpp"
#include <complex>
#include <string>
#include <utility>
#include "catch2/catch_test_macros.hpp"
#include "gpu_test_util.hpp"

namespace mrpapp {
TEST_CASE("VectorCPUTest::PointerMemoryType", "[linalg]") {
  size_t size = 3;
  size_t capacity = 11;
  std::string name("vector name");

  // Tests all the constructors.
  {
    mrpapp::Vector<float, DeviceType::CPU> vec(name, size, capacity);
    CHECK(testing::isHostPointer(vec.ptr()));
  }
  {
    Vector<int, DeviceType::CPU> vec(size);
    CHECK(testing::isHostPointer(vec.ptr()));
  }
  {
    Vector<std::complex<double>, DeviceType::CPU> vec(size, capacity);
    CHECK(testing::isHostPointer(vec.ptr()));
  }
  {
    Vector<int, DeviceType::CPU> vec(name, size);
    CHECK(testing::isHostPointer(vec.ptr()));
  }
}

TEST_CASE("VectorCPUGPUTest::Constructors", "[linalg]") {
  size_t size = 3;

  Vector<float, DeviceType::CPU> vec("name", size);
  // Set the elements.
  for (int i = 0; i < vec.size(); ++i) {
    float el = 3 * i - 2;
    vec[i] = el;
  }

  Vector<float, DeviceType::GPU> vec_copy(vec);
  CHECK(vec.size() == vec_copy.size());
  CHECK(vec.size() <= vec_copy.capacity());
  CHECK(testing::isDevicePointer(vec_copy.ptr()));

  Vector<float, DeviceType::CPU> vec_copy_copy(vec_copy);
  CHECK(vec.size() == vec_copy_copy.size());
  CHECK(vec.size() <= vec_copy_copy.capacity());
  CHECK(testing::isHostPointer(vec_copy_copy.ptr()));

  for (int i = 0; i < vec.size(); ++i) {
    CHECK(vec[i] == vec_copy_copy[i]);
    CHECK(vec.ptr(i) != vec_copy_copy.ptr(i));
  }
}

TEST_CASE("VectorCPUGPUTest::Set", "[linalg]") {
  SECTION("Assignment fits capacity")
  {
    // Assign a vector that fits into the capacity.
    size_t size = 3;

    Vector<float, DeviceType::GPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();
    Vector<float, DeviceType::CPU> vec_copy_copy(6);
    auto old_ptr_2 = vec_copy_copy.ptr();
    auto capacity_2 = vec_copy_copy.capacity();

    Vector<float, DeviceType::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    CHECK(vec.size() ==  vec_copy.size());
    CHECK(capacity ==  vec_copy.capacity());
    CHECK(old_ptr ==  vec_copy.ptr());
    CHECK(testing::isDevicePointer(vec_copy.ptr()));

    vec_copy_copy.set(vec_copy, 0, 1);
    CHECK(vec.size() ==  vec_copy_copy.size());
    CHECK(capacity_2 ==  vec_copy_copy.capacity());
    CHECK(old_ptr_2 ==  vec_copy_copy.ptr());
    CHECK(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(vec[i] ==  vec_copy_copy[i]);
      CHECK(vec.ptr(i) != vec_copy_copy.ptr(i));
    }
  }
  SECTION("Assignment exceeds capacity")
  {
    // Assign a vector that doesn't fit into the capacity.
    Vector<float, DeviceType::GPU> vec_copy(10);
    Vector<float, DeviceType::CPU> vec_copy_copy(6);
    size_t size = std::max(vec_copy.capacity(), vec_copy_copy.capacity()) + 1;

    Vector<float, DeviceType::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    CHECK(vec.size() ==  vec_copy.size());
    CHECK(vec.size() <= vec_copy.capacity());
    CHECK(!(testing::isHostPointer(vec_copy.ptr())));

    vec_copy_copy.set(vec_copy, 0, 1);
    CHECK(vec.size() ==  vec_copy_copy.size());
    CHECK(vec.size() <= vec_copy_copy.capacity());
    CHECK(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(vec[i] ==  vec_copy_copy[i]);
      CHECK(vec.ptr(i) != vec_copy_copy.ptr(i));
    }
  }
}

TEST_CASE("VectorCPUTest::setAsync", "[linalg]") {
  std::vector<int> vec(4, 1);

  Vector<int, DeviceType::GPU> vec_copy;
  Vector<int, DeviceType::CPU> vec_copy_copy;

  GpuStream stream;

  vec_copy.setAsync(vec, stream);
  vec_copy_copy.setAsync(vec_copy, stream);
  stream.sync();

  CHECK(vec.size() ==  vec_copy_copy.size());
  for (int i = 0; i < vec.size(); ++i)
    CHECK(vec[i] ==  vec_copy_copy[i]);
}

}
