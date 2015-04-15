//-*-C++-*-

#ifndef PSIMAG_MatMax_H
#define PSIMAG_MatMax_H


/** \file MatMax.h
 *
 * \brief Contains template functions to test the maxity of two
 *        matrix objects.
 *
 * The matrix objects to be compared may be of different types. The
 * types are specified by the MatType and MatTypeO template
 * arguments. These type need only provide a [] operator.
 *
 * The order of elements returned by operator[] given by the Traits
 * and TraitsO classes respectively.
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
#include "MatTraits.h"
#include "MatReduce2.h"
#include "NullType.h"

namespace psimag {


  //======================================================================
  //======================================================================
  /**
   * \brief The primary template class for the MAX<...>:EXEC
   *        function, which compares the maxity of two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  MAX<MatType, MatTypeO, Traits, TriatsO>::EXEC(m,other)
   * 
   * Forms a single large Boolean operation by assembling Boolean
   * expressions for each row from the bottom row upward and from the
   * left column to the right.
   *
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatTypeO: The type of the matrix used on the RHS.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param TraitO:   A type (e.g. MatTraits<NROW,NCOL>) describng the other matrix.
   *
   */
  template<typename MatType, 
	   typename MatTypeO, 
	   typename Traits,
	   typename TraitsO>
  class MAX {
  public:

    typedef NullType ArgType;

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };

    typedef typename Traits::ElType  T;

    class MaxReducerType {
    public:

      typedef typename Traits::ElType  T;
      typedef typename Traits::ElType  EType;
      typedef typename TraitsO::ElType OType;

      static T COMBINE(T lhs, T rhs, ArgType& arg) {
	if (lhs > rhs) 
	  return lhs;
	else
	  return rhs;
      }

      static T PRODUCE(const EType el, const OType otherEl, ArgType& arg) {
	return el;
      }
    };

    static T EXEC(const MatType& m, const MatTypeO& other) { 
      
      NullType ARG;
      return REDUCE2<MaxReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG);
    }
  };

  //======================================================================
  //======================================================================
  //======================================================================

  //======================================================================
  /**
   * \brief An matrix max operation
   *
   * \param NROW:     The number of rows in the two matricies being compared.
   * \param NCOL:     The number of columns in the two matricies being compared.
   *
   */
  template<typename T,
	   template<typename, size_t, size_t> class MatType,
	   size_t NROW, size_t NCOL >
  T Max(const MatType<T,NCOL,NROW>&  m) {
    return 
      MAX<MatType<T,NCOL,NROW>, 
      MatType<T,NCOL,NROW>, 
      typename MatType<T,NCOL,NROW>::Traits, 
      typename MatType<T,NCOL,NROW>::Traits >::EXEC(m,m);
  }


} /* namespace psimag */

#endif /* PSIMAG_Mat_Equal_H */
