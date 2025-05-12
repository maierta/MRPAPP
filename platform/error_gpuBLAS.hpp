// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file provides cublas related utilities to
// - return code checking,
// - error message printing.

#ifndef MRPAPP_ERROR_GPUBLAS_HPP
#define MRPAPP_ERROR_GPUBLAS_HPP

#include "platform/mrpapp_gpu.h"
#include "platform/mrpapp_gpu_blas.h"

#if defined(MRPAPP_HAVE_CUDA)
#include <cublas_v2.h>
#include "error_cuda.hpp"
#elif defined(MRPAPP_HAVE_HIP)
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hip/hip_complex.h>
#include "cuda2hip.h"
#include "error_hip.hpp"
#endif
#include <stdexcept>
#include <string>

namespace mrpapp {

// Returns the error string related to error.
std::string errorStringCublas(cublasStatus_t error);

// Prints an error message containing function_name, file_name, line, and the error message and
// error code related to error.
void printErrorMessage(cublasStatus_t error, std::string function_name, std::string file_name,
                       int line, std::string extra_error_string = "");

// Prints an error message and throws a std::logic_error if the return code of a cuda function is
// not CUBLAS_STATUS_SUCCESS.
// This function can be invoked with the macros checkRC and checkRCMsg (defined in error_cuda.hpp)
// that automatically include the function name, the filename, and the line to the function call.
inline void checkRCInternal(cublasStatus_t return_code, std::string function_name,
                            std::string file_name, int line, std::string extra_error_string = "") {
  if (return_code != CUBLAS_STATUS_SUCCESS) {
    printErrorMessage(return_code, function_name, file_name, line, extra_error_string);
    throw std::logic_error(function_name);
  }
}

}

#endif  // MRPAPP_ERROR_CUBLAS_HPP
