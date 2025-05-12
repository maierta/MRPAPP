// Copyright (C) 2024 ETH Zurich
// Copyright (C) 2024 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file implements cublas related utilities.

#include "platform/mrpapp_gpu.h"
#include "platform/mrpapp_gpu_blas.h"
#include "util_gpublas.hpp"
#include "handle_functions.hpp"

namespace mrpapp {
// mrpapp::

#if defined(MRPAPP_HAVE_GPU)

#if defined(MRPAPP_HAVE_CUDA)
int getGpuBLASVersion() {
  int version = 0;
  cublasStatus_t ret = cublasGetVersion(getHandle(0), &version);
  checkRC(ret);
  return version;
}
#elif defined(MRPAPP_HAVE_HIP)
int getGpuBLASVersion() {
  int version = hipblasVersionMajor;
  return version;
} 
#endif
  
#else

#endif  // MRPAPP_HAVE_CUDA

}  // namespace mrpapp
