// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the Vector<CPU> class.

#include "vector.hpp"
#include <complex>
#include <string>
#include <utility>
#include "catch2/catch_test_macros.hpp"

namespace mrpapp {

TEST_CASE("VectorCPUTest::Constructors", "[linalg]"){
  size_t size = 3;
  size_t capacity = 11;
  std::string name("vector name");

  // Tests all the constructors.
  {
    Vector<float, DeviceType::CPU> vec(name, size, capacity);
    CHECK(name ==  vec.get_name());
    CHECK(size ==  vec.size());
    CHECK(capacity <= vec.capacity());
    CHECK(nullptr != vec.ptr());
  }
  {
    Vector<double, DeviceType::CPU> vec;
    CHECK(0 == vec.size());
    CHECK(0 <= vec.capacity());
  }
  {
    Vector<int, DeviceType::CPU> vec(size);
    CHECK(size == vec.size());
    CHECK(size <= vec.capacity());
    CHECK(nullptr != vec.ptr());
  }
  {
    Vector<std::complex<double>, DeviceType::CPU> vec(size, capacity);
    CHECK(size == vec.size());
    CHECK(capacity <= vec.capacity());
    CHECK(nullptr != vec.ptr());
  }
  {
    Vector<double, DeviceType::CPU> vec(name);
    CHECK(name == vec.get_name());
    CHECK(0 == vec.size());
    CHECK(0 <= vec.capacity());
  }
  {
    Vector<int, DeviceType::CPU> vec(name, size);
    CHECK(name == vec.get_name());
    CHECK(size == vec.size());
    CHECK(size <= vec.capacity());
    CHECK(nullptr != vec.ptr());
  }
}

TEST_CASE("VectorCPUTest::ElementPointers", "[linalg]"){
  // Check if the pointers are computed correctly.
  size_t size = 5;

  Vector<int, DeviceType::CPU> vec(size);
  const Vector<int, DeviceType::CPU>& vec_const_ref(vec);
  for (int i = 0; i < vec.size(); ++i) {
    int* ptr = vec.ptr();
    CHECK(i == vec.ptr(i) - ptr);
    CHECK(vec.ptr(i) == vec_const_ref.ptr(i));
    CHECK(vec.ptr(i) == &vec[i]);
    CHECK(vec.ptr(i) == &vec_const_ref[i]);
  }
}

TEST_CASE("VectorCPUTest::ElementAccess", "[linalg]"){
  // Check if the different element accesses return the same value.
  size_t size = 4;

  Vector<int, DeviceType::CPU> vec(size);
  const Vector<int, DeviceType::CPU>& vec_const_ref(vec);
  for (int i = 0; i < vec.size(); ++i) {
    int el = 3 * i - 2;
    vec[i] = el;
    CHECK(el == vec[i]);
    CHECK(el == vec_const_ref[i]);
    CHECK(el == *(vec.ptr(i)));
    CHECK(el == *(vec_const_ref.ptr(i)));
  }
}

TEST_CASE("VectorCPUTest::CopyConstructor", "[linalg]"){
  size_t size = 4;

  Vector<float, DeviceType::CPU> vec("name", size);
  // Set the elements.
  for (int i = 0; i < vec.size(); ++i) {
    float el = 3 * i - 2;
    vec[i] = el;
  }

  Vector<float, DeviceType::CPU> vec_copy(vec, "another name");
  CHECK("another name" == vec_copy.get_name());
  CHECK(vec.size() == vec_copy.size());
  CHECK(vec.size() <= vec_copy.capacity());

  for (int i = 0; i < vec.size(); ++i) {
    CHECK(vec[i] == vec_copy[i]);
    CHECK(vec.ptr(i) != vec_copy.ptr(i));
  }
}

TEST_CASE("VectorCPUTest::Assignement", "[linalg]"){
  {
    // Assign a vector that fits into the capacity.
    size_t size = 4;

    Vector<float, DeviceType::CPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();

    Vector<float, DeviceType::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy = vec;
    CHECK(vec.size() == vec_copy.size());
    CHECK(capacity == vec_copy.capacity());
    CHECK(old_ptr == vec_copy.ptr());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(vec[i] == vec_copy[i]);
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
  {
    // Assign a vector that does not fit into the capacity.
    Vector<float, DeviceType::CPU> vec_copy(10);
    auto size = vec_copy.capacity();
    ++size;

    Vector<float, DeviceType::CPU> vec("name", size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy = vec;
    CHECK(vec.size() == vec_copy.size());
    CHECK(vec.size() <= vec_copy.capacity());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(vec[i] == vec_copy[i]);
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
}

TEST_CASE("VectorCPUTest::Set", "[linalg]"){
  {
    // Assign a vector that fits into the capacity.
    size_t size = 4;

    Vector<float, DeviceType::CPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();

    Vector<float, DeviceType::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    CHECK(vec.size() == vec_copy.size());
    CHECK(capacity == vec_copy.capacity());
    CHECK(old_ptr == vec_copy.ptr());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(vec[i] == vec_copy[i]);
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
  {
    // Assign a vector that does not fit into the capacity.
    Vector<float, DeviceType::CPU> vec_copy(10);
    auto size = vec_copy.capacity();
    ++size;

    Vector<float, DeviceType::CPU> vec("name", size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    CHECK(vec.size() == vec_copy.size());
    CHECK(vec.size() <= vec_copy.capacity());

    for (int i = 0; i < vec.size(); ++i) {
      CHECK(vec[i] == vec_copy[i]);
      CHECK(vec.ptr(i) != vec_copy.ptr(i));
    }
  }
}

TEST_CASE("VectorCPUTest::Resize", "[linalg]"){
  {
    size_t size = 4;

    Vector<long, DeviceType::CPU> vec(size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      vec[i] = el;
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
      CHECK(el == vec[i]);
    }
  }
  {
    size_t size = 5;

    Vector<long, DeviceType::CPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      vec[i] = el;
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
      CHECK(el == vec[i]);
    }
  }
  {
    size_t size = 3;

    Vector<long, DeviceType::CPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      vec[i] = el;
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
      CHECK(el == vec[i]);
    }
  }
}

TEST_CASE("VectorCPUTest::ResizeNoCopy", "[linalg]"){
  {
    size_t size = 4;

    Vector<long, DeviceType::CPU> vec(size);

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

    Vector<long, DeviceType::CPU> vec(size);
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

    Vector<long, DeviceType::CPU> vec(size);
    auto capacity = vec.capacity();

    // New size is larger than capacity.
    // Reallocation has to take place.
    size_t new_size = capacity + 1;
    vec.resizeNoCopy(new_size);
    CHECK(new_size == vec.size());
    CHECK(new_size <= vec.capacity());
  }
}

TEST_CASE("VectorCPUTest::Clear", "[linalg]"){
  Vector<double, DeviceType::CPU> vec(42);

  CHECK(42 == vec.size());
  vec.clear();
  CHECK(0 == vec.size());
  CHECK(0 == vec.capacity());
}
}
