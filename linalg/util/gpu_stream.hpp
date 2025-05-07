// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// RAII wrapper for gpu stream.

#ifndef MRPAPP_GPU_STREAM_HPP
#define MRPAPP_GPU_STREAM_HPP

#include <iostream>
#include "defines.hpp"
#include "platform/mrpapp_gpu.h"


namespace linalg {
namespace util {

#ifdef MRPAPP_HAVE_GPU

// dca::linalg::util::

class GpuStream {
public:
  GpuStream() {
    checkRC(cudaStreamCreate(&stream_));
    owning_ = true;
  }

  GpuStream(const cudaStream_t& stream) { 
    stream_ = stream;
    owning_ = false;
  }

  GpuStream(const GpuStream& other) {
    stream_ = other.stream_;
    owning_ = false;
  }

  /** simple assignment does not take possesion of the cuda stream
   */
  GpuStream& operator=(const GpuStream& other)
  {
    if (owning_ && stream_)
      checkRC(cudaStreamDestroy(stream_));
    stream_ = other.stream_;
    owning_ = false;
    return *this;
  }

  GpuStream(GpuStream&& other) noexcept {
    swap(other);
  }

  // clang at least can't do the GpuStream_t() conversion
  cudaStream_t streamActually() const {
    return stream_;
  }

  GpuStream& operator=(GpuStream&& other) noexcept {
    swap(other);
    return *this;
  }

  void sync() const {
    try {
    checkRC(cudaStreamSynchronize(stream_));
    } catch(...) {
      std::cout << "exception thrown from StreamSynchronize.\n";
    }
  }

  ~GpuStream() {
    if (owning_ && stream_)
      checkRC(cudaStreamDestroy(stream_));
  }

  operator cudaStream_t() const {
    return stream_;
  }

  void swap(GpuStream& other) noexcept {
    std::swap(stream_, other.stream_);
  }

private:
  cudaStream_t stream_ = nullptr;
  bool owning_ = false;
};

#else  // MRPAPP_HAVE_GPU

// Mock object.
class GpuStream {
public:
  GpuStream() = default;

  void sync() const {}

  // clang at least can't do the GpuStream_t() conversion
  auto streamActually(){
    return 0;
  }
};

#endif  // MRPAPP_HAVE_GPU

}  // namespace util
}  // namespace linalg

#endif  // MRPAPP_CUDA_STREAM_HPP
