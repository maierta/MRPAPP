// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides vender independent basic gpu headers.
 */
#ifndef DCA_GPU_BLAS_H
#define DCA_GPU_BLAS_H

#include "defines.hpp"
#if defined(MRPAPP_HAVE_CUDA)
#include <cublas_v2.h>
#include "error_cuda.hpp"
#elif defined(MRPAPP_HAVE_HIP)
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hip/hip_complex.h>
#include "platform/cuda2hip.h"
#include "error_hip.hpp"
#endif
#include "platform/error_gpuBLAS.hpp"

#endif
