//-*-C++-*-

#ifndef PSIMAG_MatTrace_H
#define PSIMAG_MatTrace_H


/** \file MatTrace.h
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
#include "MatReduceDiagonal.h"

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
  class TRACE {
  public:

    typedef typename Traits::ElType  T;

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    class TraceReducerType {
    public:
      
      typedef typename Traits::ElType  T;
      
      static T COMBINE(T lhs, T rhs) {
	return lhs + rhs;
      }
      
      static T PRODUCE(T el) {
	return el;
      }
    };

    static T EXEC(const MatType& m) { 
      return REDUCE_DIAG<TraceReducerType, MatType, Traits>::EXEC(m);
    }
  };


  //======================================================================
  // Trace Template Functions
  //======================================================================

  /**
   * \brief Calculates the trace of an object of type MatType.
   */
  template<typename T, size_t NROW, size_t NCOL,
	   typename Traits, template<typename, size_t, size_t, typename> class MatType>
inline  T Trace(const MatType<T,NROW,NCOL,Traits>& m) {
    return TRACE<MatType<T,NROW, NCOL, Traits>,Traits>::EXEC(m); 
  }
  
  /**  \brief Overload to calculate the 2x2 trace of a double array. */
inline  double Trace(ColMajorTraits<double,2>::ConstRefType m) { 
    return TRACE<ColMajorTraits<double,2>::Type, ColMajorTraits<double,2> >::EXEC(m); }

  /**  \brief Overload to calculate the 3x3 trace of a double array. */
inline  double Trace(ColMajorTraits<double,3>::ConstRefType m) { 
    return TRACE<ColMajorTraits<double,3>::Type, ColMajorTraits<double, 3> >::EXEC(m); 
  }

  /**  \brief Overload to calculate the 2x2 trace of an int array. */
inline  int Trace(ColMajorTraits<int,2>::ConstRefType m) { 
    return TRACE<ColMajorTraits<int,2>::Type,ColMajorTraits<int,2> >::EXEC(m); }

  /**  \brief Overload to calculate the 3x3 trace of an int array. */
inline  int Trace(ColMajorTraits<int,3>::ConstRefType m) { 
    return TRACE<ColMajorTraits<int,3>::Type,ColMajorTraits<int,3> >::EXEC(m); 
  }


} // psimag

#endif
