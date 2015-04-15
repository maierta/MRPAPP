//-*-C++-*-

#ifndef PSIMAG_MatIdentity_H
#define PSIMAG_MatIdentity_H


/** \file MatIdentity.h
 *
 * \brief Contains template functions compute the trace of a matrix.
 *
 * \warn The use of these functions and operators will be faster than
 *       their loop-coded counterparts but they will also have a much
 *       larger code size. Because of this they may not be suitable
 *       for really huge matrix-like objects.
 *
 */
 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"
#include "MatForEachDiagonal.h"

namespace psimag {

  //======================================================================
  //======================================================================
  /**
   * \brief The primary template class for the TRACE<...>:EXEC
   *        function.
   *
   * \param MatType:  The type of the matrix.
   * \param Traits:   A type (e.g. ColMajorTraits<T,NROW,NCOL>) describng the matrix.
   */
  template<typename MatType, 
	   typename Traits>
  class MakeIDENTITY {
  public:

    typedef typename Traits::ElType  T;

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    class IdentityReducerType {
    public:
      
      typedef typename Traits::ElType  T;
      
      static void OPERATE(T& el) {
	el = T(1);
      }
    };

    static void EXEC(MatType& m) { 
      FOREACH_DIAG<IdentityReducerType, MatType, Traits>::EXEC(m);
    }
  };

  template<typename T, size_t NROW, size_t NCOL,
	   template<typename, size_t, size_t>           class TraitsTemplate,
	   template<typename, size_t, size_t, typename> class MatType>
  void MakeIdentity(MatType<T,NROW,NCOL,TraitsTemplate<T,NROW,NCOL> >& m) {
    MakeIDENTITY<MatType<T,NROW,NCOL,TraitsTemplate<T,NROW,NCOL> > ,TraitsTemplate<T,NROW,NCOL> >::EXEC(m);
  }
} // psimag
#endif
