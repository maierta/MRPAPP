// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides include all types of allocators, and provides a default selector.

#ifndef MRPAPP_LINALG_UTIL_ALLOCATORS_HPP
#define MRPAPP_LINALG_UTIL_ALLOCATORS_HPP

#include <stdexcept>
#include "aligned_allocator.hpp"
#include "device_type.hpp"
#ifdef MRPAPP_HAVE_GPU
#include "device_allocator.hpp"
#include "managed_allocator.hpp"
#include "pinned_allocator.hpp"
#endif  // MRPAPP_HAVE_GPU

namespace mrpapp {
namespace selector {
// mrpapp::selector::
template <typename T, DeviceType device>
struct DefaultAllocator;

#ifdef MRPAPP_HAVE_GPU
template <typename T>
struct DefaultAllocator<T, DeviceType::CPU> {
  using type = PinnedAllocator<T>;
};

template <typename T>
struct DefaultAllocator<T, DeviceType::GPU> {
  using type = DeviceAllocator<T>;
};
#else

template <typename T>
struct DefaultAllocator<T, DeviceType::CPU> {
  using type = AlignedAllocator<T>;
};

template <typename T>
struct DefaultAllocator<T, DeviceType::GPU> {
  struct UnusedAllocator {
    T* allocate(std::size_t) {
      throw(std::logic_error("GPU not available."));
    }
    void deallocate(T*& /*ptr*/, std::size_t /*n*/ = 0) {}
  };
  using type = UnusedAllocator;
};

#endif  // MRPAPP_HAVE_GPU

}  // selector

template <typename T, DeviceType device>
using DefaultAllocator = typename selector::DefaultAllocator<T, device>::type;

} 

#endif  // MRPAPP_LINALG_UTIL_ALLOCATORS_HPP
