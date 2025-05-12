// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides better type mapping between host and gpu types
 */

#ifndef MRPAPP_GPU_TYPE_MAPPING_HPP
#define MRPAPP_GPU_TYPE_MAPPING_HPP

#include <complex>
#include <type_traits>
#include <string>
#include <memory>
#include "util/type_mapping.hpp"
#include "util/type_fundamentals.hpp"
#include "platform/mrpapp_gpu_complex.h"

namespace mrpapp {

#ifdef MRPAPP_HAVE_GPU

/** Type maps to handle conversion of complex types from GPU to std C++ representation
 *  representations.
 */
template <typename T>
using HOSTTypeMap =
    typename std::disjunction<OnTypesEqual<T, float, float>, OnTypesEqual<T, double, double>,
                              OnTypesEqual<T, const float, const float>,
                              OnTypesEqual<T, const double, const double>,
                              OnTypesEqual<T, float2, std::complex<float>>,
                              OnTypesEqual<T, double2, std::complex<double>>, default_type<void>>::type;

/** Type maps to handle cast from of complex type pointers from GPU to std C++ representation
 *  representations.
 */
template <typename T>
using HOSTPointerMap = typename std::disjunction<

    OnTypesEqual<T, float, float>, OnTypesEqual<T, double, double>, OnTypesEqual<T, float*, float*>,
    OnTypesEqual<T, double*, double*>, OnTypesEqual<T, const float*, const float*>,
    OnTypesEqual<T, const double*, const double*>, OnTypesEqual<T, float**, float**>,
    OnTypesEqual<T, double**, double**>, OnTypesEqual<T, const float**, const float**>,
    OnTypesEqual<T, const double**, const double**>, OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
    OnTypesEqual<T, float2*, std::complex<float>*>, OnTypesEqual<T, double2*, std::complex<double>*>,
    OnTypesEqual<T, float2**, std::complex<float>**>, OnTypesEqual<T, double2**, std::complex<double>**>,
    OnTypesEqual<T, const float2*, const std::complex<float>*>,
    OnTypesEqual<T, const double2*, const std::complex<double>*>,
    OnTypesEqual<T, const float2**, const std::complex<float>**>,
    OnTypesEqual<T, const double2**, const std::complex<double>**>,
    OnTypesEqual<T, std::complex<float>*, std::complex<float>*>,
    OnTypesEqual<T, std::complex<double>*, std::complex<double>*>,
    OnTypesEqual<T, std::complex<float>**, std::complex<float>**>,
    OnTypesEqual<T, std::complex<double>**, std::complex<double>**>,
    OnTypesEqual<T, const std::complex<float>*, const std::complex<float>*>,
    OnTypesEqual<T, const std::complex<double>*, const std::complex<double>*>,
    OnTypesEqual<T, const std::complex<float>**, const std::complex<float>**>,
    OnTypesEqual<T, const std::complex<double>**, const std::complex<double>**>,
    default_type<void>>::type;

/** Type maps to handle conversion of complex types from std c++ to GPU representation
 *  representations.
 */
template <typename T>
using CUDATypeMap = typename std::disjunction<
    OnTypesEqual<T, float, float>, OnTypesEqual<T, double, double>, OnTypesEqual<T, float*, float*>,
    OnTypesEqual<T, double*, double*>, OnTypesEqual<T, const float*, const float*>,
    OnTypesEqual<T, const double*, const double*>, OnTypesEqual<T, float**, float**>,
    OnTypesEqual<T, const float**, const float**>, OnTypesEqual<T, double**, double**>,
    OnTypesEqual<T, const double**, const double**>, OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
    OnTypesEqual<T, std::complex<float>, cuComplex>, OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
    OnTypesEqual<T, std::complex<double>*, cuDoubleComplex*>,
    OnTypesEqual<T, std::complex<float>**, cuComplex**>,
    OnTypesEqual<T, std::complex<double>**, cuDoubleComplex**>,
    OnTypesEqual<T, std::complex<float>*, cuComplex*>, OnTypesEqual<T, float2, cuComplex>,
    OnTypesEqual<T, double2, cuDoubleComplex>,
    OnTypesEqual<T, const std::complex<double>*, const cuDoubleComplex*>,
    OnTypesEqual<T, const std::complex<float>*, const cuComplex*>,
    OnTypesEqual<T, const std::complex<double>&, const cuDoubleComplex&>,
    OnTypesEqual<T, const std::complex<float>&, const cuComplex&>,
    OnTypesEqual<T, const std::complex<float>**, const cuComplex**>,
    OnTypesEqual<T, const std::complex<double>**, const cuDoubleComplex**>,
    OnTypesEqual<T, const std::complex<float>* const*, const cuComplex* const*>,
    OnTypesEqual<T, const std::complex<double>* const*, const cuDoubleComplex* const*>,
    default_type<void>>::type;

template <typename T>
__device__ __host__ CUDATypeMap<T> castGPUType(T var) {
  return reinterpret_cast<CUDATypeMap<T>>(var);
}

template <typename T>
using CUDARealAliasMap =
    typename std::disjunction<OnTypesEqual<T, float2, float>, OnTypesEqual<T, double2, double>,
                              OnTypesEqual<T, cuDoubleComplex, double>,
                              OnTypesEqual<T, cuComplex, float>, default_type<void>>::type;

template <typename T>
CUDARealAliasMap<T> realAliasGPU(T var) {
  return reinterpret_cast<CUDARealAliasMap<T>>(var);
}

template <typename T>
struct IsCUDAComplex_t
    : std::disjunction<std::is_same<float2, T>, std::is_same<double2, T>, std::false_type> {};
  
template <typename T>
using IsCUDAComplex = std::enable_if_t<IsCUDAComplex_t<std::decay_t<T>>::value, bool>;
  
template <typename Real>
struct Real2CudaComplex;

template <>
struct Real2CudaComplex<double> {
  using type = cuDoubleComplex;
};
template <>
struct Real2CudaComplex<float> {
  using type = cuComplex;
};

#ifdef MRPAPP_HAVE_CUDA
template <typename Real>
using GPUComplex = typename Real2CudaComplex<Real>::type;
#endif

template <typename Real>
using CUDAComplex = typename Real2CudaComplex<Real>::type;

#endif

#ifdef MRPAPP_HAVE_GPU
template <typename T>
__device__ __host__ HOSTPointerMap<T> castHostType(T var) {
  if constexpr (std::is_same_v<HOSTPointerMap<T>, T>)
    return var;
  else if constexpr (std::is_pointer_v<T>)
    return reinterpret_cast<HOSTPointerMap<T>>(var);
}
#endif

}  // namespace mrpapp

#endif
