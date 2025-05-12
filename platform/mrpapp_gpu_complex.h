// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides working vender or magma complex gpu headers.
 *  Order of includes is quite brittle for complex types and operators,
 *  check all platforms when making any changes.
 */
#ifndef MRPAPP_GPU_COMPLEX_H
#define MRPAPP_GPU_COMPLEX_H
#include <type_traits>
#include <complex>

#ifdef MRPAPP_HAVE_CUDA
#include <cuComplex.h>
#endif

#include "util/type_fundamentals.hpp"

#if defined(MRPAPP_HAVE_HIP)
#include "platform/cuda2hip.h"

namespace linalg {

template <typename T>
__device__ __host__ inline void assign(T& a, const T b) {
  a = b;
}

}  // namespace linalg
#endif

namespace linalg {
#ifdef MRPAPP_HAVE_GPU
// The contents of the cast come from en.cppreference.com/w/cpp/numeric/complex
template <typename T>
__device__ __host__ inline void assign(std::complex<T>& a, const T b) {
  a = {b, 0.0};
}

__device__ __host__ inline void assign(double2& a, const int8_t b) {
  a.x = static_cast<double>(b);
  a.y = 0.0;
}

__device__ __host__ inline void assign(float2& a, const int8_t b) {
  a.x = static_cast<float>(b);
  a.y = 0.0;
}
#endif
}  // namespace linalg

#endif
