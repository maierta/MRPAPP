//-*-C++-*-

#ifndef PSIMAG_MatMagnitude_H
#define PSIMAG_MatMagnitude_H


/** \file MatMagnitude.h
 *
 * \brief Contains template functions to test the magnitudeity of two
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
   * \brief The primary template class for the MAGNITUDE<...>:EXEC
   *        function, which compares the magnitudeity of two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  MAGNITUDE<MatType, MatTypeO, Traits, TriatsO>::EXEC(m,other)
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
  class MAGNITUDE {
  public:

    typedef NullType ArgType;

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };
    
    class MagnitudeReducerType {
    public:

      typedef typename Traits::ElType  T;
      typedef typename Traits::ElType  EType;
      typedef typename TraitsO::ElType OType;

      static T COMBINE(T lhs, T rhs, ArgType& arg) {

	return lhs + rhs;
      }

      static T PRODUCE(const EType el, const OType otherEl, ArgType& arg) {

	return static_cast<T>(sqrt(static_cast<double>(el * otherEl)));
      }
    };
    
    static typename Traits::ElType EXEC(const MatType& m, const MatTypeO& other) { 
      ArgType ARG;
      return REDUCE2<MagnitudeReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG);
      //typename Traits::ElType result = REDUCE2<MagnitudeReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG);

      //return result;
    }
  };

  //======================================================================
  //======================================================================
  /**
   * \brief The primary template class for the MAGNITUDE<...>:EXEC
   *        function, which compares the magnitudeity of two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  MAGNITUDE<MatType, MatTypeO, Traits, TriatsO>::EXEC(m,other)
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
  class L2NORM {
  public:

    typedef NullType ArgType;

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };
    
    class MagnitudeReducerType {
    public:

      typedef typename Traits::ElType  T;
      typedef typename Traits::ElType  EType;
      typedef typename TraitsO::ElType OType;

      static T COMBINE(T lhs, T rhs, ArgType& arg) {

	return lhs + rhs;
      }

      static T PRODUCE(const EType el, const OType otherEl, ArgType& arg) {

	return static_cast<T>(static_cast<double>(el * otherEl));
      }
    };
    
    static typename Traits::ElType EXEC(const MatType& m, const MatTypeO& other) { 
      ArgType ARG;
      return sqrt(REDUCE2<MagnitudeReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG));
      //typename Traits::ElType result = REDUCE2<MagnitudeReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG);

      //return result;
    }
  };

  //======================================================================
  //======================================================================
  //======================================================================

  //======================================================================
  /**
   * \brief An matrix magnitude operation
   *
   * \param NROW:     The number of rows in the two matricies being compared.
   * \param NCOL:     The number of columns in the two matricies being compared.
   *
   */
  template<typename T, size_t NROW, size_t NCOL,
	   template<typename, size_t, size_t> class TraitsTemplate,
	   template<typename, size_t, size_t, typename> class MatType>
  T Magnitude(const MatType<T,NCOL,NROW, TraitsTemplate<T, NCOL,NROW> >&  m) {
    T result =  
      MAGNITUDE<MatType<T,NCOL,NROW, TraitsTemplate<T, NCOL,NROW> >, 
      MatType<T,NCOL,NROW, TraitsTemplate<T, NCOL,NROW> >, 
      TraitsTemplate<T, NCOL,NROW>, 
      TraitsTemplate<T, NCOL,NROW> >::EXEC(m,m);

    return result;
  }


} /* namespace psimag */

#endif /* PSIMAG_Mat_Equal_H */
