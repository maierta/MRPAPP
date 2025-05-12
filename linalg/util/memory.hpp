// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides memory related utility:
// - setToZero.

#ifndef MRPAPP_LINALG_UTIL_MEMORY_HPP
#define MRPAPP_LINALG_UTIL_MEMORY_HPP

#include <cassert>
#include <complex>
#include <cstring>
#include <stdexcept>

#include "platform/mrpapp_gpu.h"
#include "util/type_help.hpp"
#include "device_type.hpp"
#include "gpu_stream.hpp"
#include "util/ignore.hpp"

#ifdef MRPAPP_HAVE_GPU
#include "gpu_type_mapping.hpp"
#endif

namespace mrpapp {
// mrpapp::

template <DeviceType device_name>
struct Memory {};

template <>
struct Memory<CPU> {
  // Sets the elements to 0. Only defined for arithmetic types and
  // std::complex of aritmetic types.
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      ScalarType* ptr, size_t size) {
    std::memset(ptr, 0, size * sizeof(ScalarType));
  }
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      std::complex<ScalarType>* ptr, size_t size) {
    std::complex<ScalarType> c_zero{0.0, 0.0};
    std::fill_n(
        ptr, size,
        c_zero);  // memset(static_cast<void*>(ptr), 0, size * sizeof(std::complex<ScalarType>));
  }
  template <typename ScalarType>
  static void setToZeroAsync(ScalarType* ptr, size_t size, const GpuStream& /*s*/) {
    setToZero(ptr, size);
  }
  template <typename ScalarType>
  static void setToZero(ScalarType* ptr, size_t size, const GpuStream& /*s*/) {
    setToZero(ptr, size);
  }
#ifdef MRPAPP_HAVE_GPU
  template <typename Scalar>
  static std::enable_if_t<util::IsCUDAComplex_t<Scalar>::value == true, void> setToZero(
      Scalar* ptr, size_t size) {
    std::memset(ptr, 0, sizeof(Scalar) * size);
  }

#ifdef MRPAPP_HAVE_HIP
  template <typename Scalar>
  static std::enable_if_t<util::IsMagmaComplex_t<Scalar>::value == true, void> setToZero(
      Scalar* ptr, size_t size) {
    std::memset(ptr, 0, sizeof(Scalar) * size);
  }
#endif
#endif
};

#ifdef MRPAPP_HAVE_GPU
template <>
struct Memory<GPU> {
  // Sets the elements to 0. Only defined for arithmetic types and
  // std::complex of aritmetic types.

  /// Specialization for float2, double2, cuComplex, cuDoubleComplex
  template <typename ScalarType>
  static std::enable_if_t<dca::util::IsCUDAComplex_t<ScalarType>::value == true, void> setToZero(ScalarType ptr, size_t size) {
    checkRC(cudaMemset(ptr, 0, size * sizeof(ScalarType)));
  }

  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      ScalarType* ptr, size_t size) {
    checkRC(cudaMemset(ptr, 0, size * sizeof(ScalarType)));
  }
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      std::complex<ScalarType>* ptr, size_t size) {
    checkRC(cudaMemset(ptr, 0, size * sizeof(std::complex<ScalarType>)));
  }

  template <typename Scalar>
  static std::enable_if_t<dca::util::IsCUDAComplex_t<Scalar>::value == true, void> setToZero(
      Scalar* ptr, size_t size) {
    checkRC(cudaMemset(ptr, 0, size * sizeof(Scalar)));
  }


  

  // Do nothing for non arithmetic types.
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == false, void> setToZero(
      ScalarType /*ptr*/, size_t /*size*/) {}

  template <typename ScalarType>
  static void setToZeroAsync(ScalarType* ptr, size_t size, const GpuStream& stream) {
    checkRC(cudaMemsetAsync(ptr, 0, size * sizeof(ScalarType), stream));
  }

  template <typename ScalarType>
  static void setToZero(ScalarType* ptr, size_t size, const GpuStream& stream) {
    checkRC(cudaMemsetAsync(ptr, 0, size * sizeof(ScalarType), stream));
    cudaEvent_t zero_event;
    checkRC(cudaEventCreateWithFlags(&zero_event, cudaEventBlockingSync));
    checkRC(cudaEventRecord(zero_event, stream));
    checkRC(cudaEventSynchronize(zero_event));
  }
};
#endif  // MRPAPP_HAVE_GPU

}  // namespace mrpapp

#endif  // MRPAPP_LINALG_UTIL_MEMORY_HPP
