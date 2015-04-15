//-*-C++-*-

#ifndef PSIMAG_MatForEachDiagonal_H
#define PSIMAG_MatForEachDiagonal_H


/** \file MatForEachDiagonal.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         ForEach Diagonals operation. This operation applies a given
 *         function to every diagonal element of the given matrix. It
 *         then collects the results into a single value using another
 *         function.
 *
 * The matrix objects to be operated on may be of different types. The
 * types are specified by template arguments like MatType. These types
 * need only be consistant with the given MatTraits class.
 *
 * \param The MatTraits class which must provide a templated class
 *        REF<MatType, Row,Col> with a static function member
 *        GET(MatType& m). This method defines a code fragment for
 *        returning a reference to the given ROW and COL.
 *
 * Using this metaprogramming technique we construct loop-less matrix
 * diagonal reduce operation for any size matricies.
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

namespace psimag {


  //======================================================================
  /**
   * \brief The primary template class for assembling the
   *        FOREACH_DIAG_<...>:EXEC function.
   *
   * It is initially invoked like: 
   *
   *  FOREACH_DIAG_<FunctorType, MatType, Traits, NROW-1>::EXEC(m,other)
   * 
   * \param FunctorType:  The type which provides the functor which provides a return type and two functions.
   * \param MatType:      The type of the matrix used on the LHS.
   * \param Traits:       A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param ROW:          The current row (zero indexed) to assemble the boolean expression for.
   *
   */
  template<typename Functor, 
	   typename MatType, 
	   typename Traits,
	   size_t   ROW>
  class FOREACH_DIAG_{
  public:
    
    typedef typename Traits::template  REF<MatType,ROW,ROW> REF_;

    static void EXEC(MatType& m) {

      Functor::OPERATE(REF_::GET(m));
      FOREACH_DIAG_ <Functor, MatType, Traits, ROW-1>::EXEC(m);
    }

    static void EXEC(const MatType& m) {

      Functor::OPERATE(REF_::GETCONST(m));
      FOREACH_DIAG_ <Functor, MatType, Traits, ROW-1>::EXEC(m);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH_DIAG_ELS<...>:EXEC.
   *
   * The FOREACH_DIAG_ELS<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename Functor, 
	   typename MatType, 
	   typename Traits>
  class FOREACH_DIAG_<Functor, MatType, Traits, 0> {
  public:

    typedef typename Functor::T                     ReturnType;
    typedef typename Traits::template  REF<MatType,0,0> REF_;
  
    static void EXEC(MatType& m) {

      Functor::OPERATE(REF_::GET(m));
    }

    static void EXEC(const MatType& m) {

      Functor::OPERATE(REF_::GETCONST(m));
    }
  };
  
  //======================================================================
  /**
   * \brief The ForEach template that calls FOREACH_DIAG_ to start the
   *        recursive computation of the aggregate reduce operation at
   *        row NROW-1.
   *
   */
  template<typename Functor,
	   typename MatType, 
	   typename Traits>
  class FOREACH_DIAG {
  public:
    
    typedef typename Functor::T ReturnType; 

    enum { NROW=Traits::NROW };

    static void EXEC(MatType& m) {

      return FOREACH_DIAG_<Functor, MatType, Traits, NROW-1>::EXEC(m);
    }

    static void EXEC(const MatType& m) {

      return FOREACH_DIAG_<Functor, MatType, Traits, NROW-1>::EXEC(m);
    }
  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_ForEach_H */
