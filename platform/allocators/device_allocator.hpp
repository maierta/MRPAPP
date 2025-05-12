// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides allocators with pinned or mapped memory usable with std::vector.

#ifndef MRPAPP_LINALG_UTIL_ALLOCATORS_DEVICE_ALLOCATOR_HPP
#define MRPAPP_LINALG_UTIL_ALLOCATORS_DEVICE_ALLOCATOR_HPP

#ifdef MRPAPP_HAVE_GPU
#include "platform/mrpapp_gpu.h"
#else
#pragma error "This file requires GPU support."
#endif

namespace mrpapp {
// mrpapp::

template <typename T>
class DeviceAllocator {
protected:
  T* allocate(std::size_t n) {
    if (n == 0)
      return nullptr;
    T* ptr;
    cudaError_t ret = cudaMalloc((void**)&ptr, n * sizeof(T));
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__,
                        "\t DEVICE size requested : " + std::to_string(n * sizeof(T)));
      throw(std::bad_alloc());
    }
    return ptr;
  }

  void deallocate(T*& ptr, std::size_t /*n*/ = 0) noexcept {
    cudaError_t ret = cudaFree(ptr);
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__);
      std::terminate();
    }
    ptr = nullptr;
  }

public:
  // SFINAE method for setting managed memory stream.
  void setStream(const cudaStream_t /*stream*/) const {}
};

}

#endif  // MRPAPP_LINALG_UTIL_ALLOCATORS_DEVICE_ALLOCATOR_HPP
