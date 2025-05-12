// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the Vector<GPU> class.

#include "vector.hpp"
#include <complex>
#include <string>
#include <utility>
#include "catch2/catch_test_macros.hpp"
#include "gpu_test_util.hpp"

namespace mrpapp {
TEST_CASE("VectorGPUTest::Constructors", "[linalg]"){
  size_t size = 3;
  size_t capacity = 11;
  std::string name("vector name");

  // Tests all the constructors.
  {
    Vector<float, DeviceType::GPU> vec(name, size, capacity);
    CHECK(name == vec.get_name());
    CHECK(size == vec.size());
    CHECK(capacity <= vec.capacity());
    CHECK(nullptr != vec.ptr());
    CHECK(testing::isDevicePointer(vec.ptr()));
  }
  {
    Vector<double, DeviceType::GPU> vec;
    CHECK(0 == vec.size());
    CHECK(0 <= vec.capacity());
  }
  {
    Vector<int, DeviceType::GPU> vec(size);
    CHECK(size == vec.size());
    CHECK(size <= vec.capacity());
    CHECK(nullptr != vec.ptr());
    CHECK(testing::isDevicePointer(vec.ptr()));
  }
  {
    Vector<std::complex<double>, DeviceType::GPU> vec(size, capacity);
    CHECK(size == vec.size());
    CHECK(capacity <= vec.capacity());
    CHECK(nullptr != vec.ptr());
    CHECK(testing::isDevicePointer(vec.ptr()));
  }
  {
    Vector<double, DeviceType::GPU> vec(name);
    CHECK(name == vec.get_name());
    CHECK(0 == vec.size());
    CHECK(0 <= vec.capacity());
  }
  {
    Vector<int, DeviceType::GPU> vec(name, size);
    CHECK(name == vec.get_name());
    CHECK(size == vec.size());
    CHECK(size <= vec.capacity());
    CHECK(nullptr != vec.ptr());
    CHECK(testing::isDevicePointer(vec.ptr()));
  }
}

TEST_CASE("VectorGPUTest::ElementPointers", "[linalg]"){
  // Check if the pointers are computed correctly.
  size_t size = 5;

  Vector<int, DeviceType::GPU> vec(size);
  const Vector<int, DeviceType::GPU>& vec_const_ref(vec);
  for (int i = 0; i < vec.size(); ++i) {
    int* ptr = vec.ptr();
    CHECK(i == vec.ptr(i) - ptr);
    CHECK(vec.ptr(i) == vec_const_ref.ptr(i));
  }
}

TEST_CASE("VectorGPUTest::CopyConstructor", "[linalg]"){
  size_t size = 4;

  Vector<float, DeviceType::GPU> vec("name", size);
  // Set the elements.
  for (int i = 0; i < vec.size(); ++i) {
    float el = 3 * i - 2;
    testing::setOnDevice(vec.ptr(i), el);
  }

  Vector<float, DeviceType::GPU> vec_copy(vec);
  CHECK(vec.size() == vec_copy.size());
  CHECK(vec.size() <= vec_copy.capacity());

  for (int i = 0; i < vec.size(); ++i) {
    CHECK(testing::getFromDevice(vec.ptr(i)) == testing::getFromDevice(vec_copy.ptr(i)));
    CHECK(vec.ptr(i) != vec_copy.ptr(i));
  }
}

TEST_CASE("VectorGPUTest::Assignement", "[linalg]"){
  {
    // Assign a vector that fits into the capacity.
    size_t size = 4;

    Vector<float, DeviceType::GPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();

    Vector<float, DeviceType::GPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      testing::setOnDevice(vec.ptr(i), el);
    }

    vec_copy = vec;
    CHECK(vec.size() == vec_copy.size());
    CHECK(capacity == vec_copy.capacity());
    CHECK(old_ptr == vec_copy.ptr());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(testing::getFromDevice(vec.ptr(i)) == testing::getFromDevice(vec_copy.ptr(i)));
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
  {
    // Assign a vector that does not fit into the capacity.
    Vector<float, DeviceType::GPU> vec_copy(10);
    auto size = vec_copy.capacity();
    ++size;

    Vector<float, DeviceType::GPU> vec("name", size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      testing::setOnDevice(vec.ptr(i), el);
    }

    vec_copy = vec;
    CHECK(vec.size() == vec_copy.size());
    CHECK(vec.size() <= vec_copy.capacity());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(testing::getFromDevice(vec.ptr(i)) == testing::getFromDevice(vec_copy.ptr(i)));
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
}

TEST_CASE("VectorGPUTest::Set", "[linalg]"){
  {
    // Assign a vector that fits into the capacity.
    size_t size = 4;

    Vector<float, DeviceType::GPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();

    Vector<float, DeviceType::GPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      testing::setOnDevice(vec.ptr(i), el);
    }

    vec_copy.set(vec, 0, 1);
    CHECK(vec.size() == vec_copy.size());
    CHECK(capacity == vec_copy.capacity());
    CHECK(old_ptr == vec_copy.ptr());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(testing::getFromDevice(vec.ptr(i)) == testing::getFromDevice(vec_copy.ptr(i)));
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
  {
    // Assign a vector that does not fit into the capacity.
    Vector<float, DeviceType::GPU> vec_copy(10);
    auto size = vec_copy.capacity();
    ++size;

    Vector<float, DeviceType::GPU> vec("name", size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      testing::setOnDevice(vec.ptr(i), el);
    }

    vec_copy.set(vec, 0, 1);
    CHECK(vec.size() == vec_copy.size());
    CHECK(vec.size() <= vec_copy.capacity());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(testing::getFromDevice(vec.ptr(i)) == testing::getFromDevice(vec_copy.ptr(i)));
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
}

TEST_CASE("VectorGPUTest::Resize", "[linalg]"){
  {
    size_t size = 4;

    Vector<long, DeviceType::GPU> vec(size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      testing::setOnDevice(vec.ptr(i), el);
    }

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    int new_size = capacity;
    vec.resize(new_size);
    CHECK(new_size == vec.size());
    CHECK(capacity == vec.capacity());
    CHECK(old_ptr == vec.ptr());

    // Check the value of the elements.
    for (int i = 0; i < size; ++i) {
      long el = 1 + 3 * i;
      CHECK(el == testing::getFromDevice(vec.ptr(i)));
    }
  }
  {
    size_t size = 5;

    Vector<long, DeviceType::GPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      testing::setOnDevice(vec.ptr(i), el);
    }

    // Shrink the vector. No reallocation has to take place.
    int new_size = 2;
    vec.resize(new_size);
    CHECK(new_size == vec.size());
    CHECK(capacity == vec.capacity());
    CHECK(old_ptr == vec.ptr());

    // Check the value of the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      CHECK(el == testing::getFromDevice(vec.ptr(i)));
    }
  }
  {
    size_t size = 3;

    Vector<long, DeviceType::GPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      testing::setOnDevice(vec.ptr(i), el);
    }

    // New size is larger than capacity.
    // Reallocation has to take place.
    int new_size = capacity + 1;
    vec.resize(new_size);
    CHECK(new_size == vec.size());
    CHECK(new_size <= vec.capacity());
    CHECK(old_ptr != vec.ptr());

    // Check the value of the elements.
    for (int i = 0; i < size; ++i) {
      long el = 1 + 3 * i;
      CHECK(el == testing::getFromDevice(vec.ptr(i)));
    }
  }
}

TEST_CASE("VectorGPUTest::ResizeNoCopy", "[linalg]"){
  {
    size_t size = 4;

    Vector<long, DeviceType::GPU> vec(size);

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    size_t new_size = capacity;
    vec.resizeNoCopy(new_size);
    CHECK(new_size == vec.size());
    CHECK(capacity == vec.capacity());
    CHECK(old_ptr == vec.ptr());
  }
  {
    size_t size = 5;

    Vector<long, DeviceType::GPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();

    // Shrink the vector. No reallocation has to take place.
    size_t new_size = 2;
    vec.resizeNoCopy(new_size);
    CHECK(new_size == vec.size());
    CHECK(capacity == vec.capacity());
    CHECK(old_ptr == vec.ptr());
  }
  {
    size_t size = 3;

    Vector<long, DeviceType::GPU> vec(size);
    auto capacity = vec.capacity();

    // New size is larger than capacity.
    // Reallocation has to take place.
    size_t new_size = capacity + 1;
    vec.resizeNoCopy(new_size);
    CHECK(new_size == vec.size());
    CHECK(new_size <= vec.capacity());
  }
}

TEST_CASE("VectorGPUTest::Clear", "[linalg]"){
  Vector<double, DeviceType::GPU> vec(42);

  CHECK(42 == vec.size());
  vec.clear();
  CHECK(0 == vec.size());
  CHECK(0 == vec.capacity());
}
}
