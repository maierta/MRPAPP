//-*-C++-*-

#ifndef PSIMAG_MatForEach_H
#define PSIMAG_MatForEach_H


/** \file MatForEach.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         ForEach operation. This operation applies a given function
 *         to every tuple of elements in a given set of matricies.
 *
 * The ForEach operation is constructed as a static method on a
 * templated ForEach class:
 *
 * e.g.  FOREACH_<FunctorType, MatType, MatTypeRHS, NROW, NCOL, ROW>::EXEC(lhs, rhs)
 *
 * This static function is constructed at compile time through
 * recursive invokations of itself (and it's row == 0 specialization) and
 * through invokations of:
 *
 * FOREACH_ELS<FunctorType, MatType, MatTypeRHS, NROW, NCOL, ROW, COL>::EXEC(lhs, rhs)
 *
 * and it's specializations. The specializations are used to halt the
 * recursive construction of the FOREACH_<...>::EXEC function.
 *
 * Using this metaprogramming technique we construct loop-less matrix
 * foreach operation for any size matricies.
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
   *        FOREACH_ELS<...>:EXEC function, which operates on the
   *        corresponding rows of matrix-like objects (not necessarily
   *        of the same type).
   *
   * It is initially invoked like: 
   *
   *  FOREACH_ELS<FunctorType, MatType, MatTypeRHS, Traits, TraitsRHS, ROW, NCOL-1>::EXEC(lhs, rhs)
   * 
   * This operation traverses the matrix in row-major order (for the sake of printing)>
   *
   * \param FunctorType:  The type which provides the operation through it's EXEC statis methos.
   * \param MatType:  The type of the matrix used on the LHS.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param ROW:      The current row (zero indexed) to assemble the boolean expression for.
   * \param COL:      The current row (zero indexed) to assemble the boolean expression for.
   *
   */
  template<typename FunctorType, 
	   typename MatType, 
	   typename Traits,
	   typename ArgType,
	   size_t ROW, size_t COLS_LEFT>
  class FOREACH_ELS {
  public:

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    static const size_t COL = NCOL - COLS_LEFT;

    typedef typename Traits::template REF<MatType,   ROW,COL> REF_LHS;

    static void EXEC(const MatType& lhs, ArgType& arg) {

      FunctorType::EXEC(REF_LHS::GETCONST(lhs), ROW, COL, arg);

      FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, COLS_LEFT-1>::EXEC(lhs, arg);
    }

    static void EXEC( MatType& lhs, ArgType& arg) {

      FunctorType::EXEC(REF_LHS::GETCONST(lhs), ROW, COL, arg);

      FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, COLS_LEFT-1>::EXEC(lhs, arg);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH_ELS<...>:EXEC.
   *
   * The FOREACH_ELS<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename FunctorType, 
	   typename MatType, 
	   typename Traits,
	   typename ArgType,
	   size_t ROW>
  class FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, 1> {
  public:
    
    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    static const size_t COL = NCOL - 1;

    typedef typename Traits::template REF<MatType,   ROW,COL> REF_LHS;
    
    static void EXEC(const MatType& lhs, ArgType& arg) {

      FunctorType::EXEC(REF_LHS::GETCONST(lhs), ROW, COL, arg);
    }

    static void EXEC(MatType& lhs, ArgType& arg) {

      FunctorType::EXEC(REF_LHS::GETCONST(lhs), ROW, COL, arg);
    }
  };
  
  //======================================================================
  /**
   * \brief The primary template class for assembling the
   *        FOREACH_<...>:EXEC function, which operates on
   *        two marix-like objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  FOREACH_<FunctorType,MatType, MatTypeRHS, NROW, NCOL, NROW-1>::EXEC(lhs, rhs)
   * 
   * It recursively forms a single operation by
   * assembling operations for each row from the bottom row
   * upward. The recursion is terminated by a specialization of
   * FOREACH_<...>:EXEC for row zero.
   *
   * \param FunctorType:  The type which provides the operation through it's EXEC statis method.
   * \param MatType:  The type of the matrix used on the LHS.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describing the matrix.
   * \param ROW:      The current row (zero indexed).
   *
   */
  template<typename FunctorType,
	   typename MatType, 
	   typename Traits,
	   typename ArgType,
	   size_t ROWS_LEFT>
  class FOREACH_ {
  public:

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };

    static const size_t ROW = NROW - ROWS_LEFT;
    
    static void EXEC(const MatType& lhs, ArgType& arg) { 

      FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, NCOL>::EXEC(lhs, arg);
      FOREACH_<FunctorType, MatType, Traits, ArgType, ROWS_LEFT-1>::EXEC(lhs, arg);
    }

    static void EXEC(MatType& lhs, ArgType& arg) { 

      FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, NCOL>::EXEC(lhs, arg);
      FOREACH_<FunctorType, MatType, Traits, ArgType, ROWS_LEFT-1>::EXEC(lhs, arg);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH_<...>:EXEC.
   *
   * The FOREACH_<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename FunctorType,
	   typename MatType, 
	   typename Traits,
	   typename ArgType>
  class FOREACH_<FunctorType, MatType, Traits, ArgType, 1> {
  public:
    
    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };

    static const size_t ROW = NROW - 1;

    static void EXEC(const MatType& lhs, ArgType& arg) {

      FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, NCOL>::EXEC(lhs, arg);
    }

    static void EXEC( MatType& lhs, ArgType& arg) {

      FOREACH_ELS<FunctorType, MatType, Traits, ArgType, ROW, NCOL>::EXEC(lhs, arg);
    }
  };

  //======================================================================
  /**
   * \brief The FOREACH template that calls FOREACH_ to start the
   *        recursive computation of the aggregate operation at
   *        ROWS_LEFT =  NROW.
   *
   */
  template<typename FunctorType,
	   typename MatType, 
	   typename Traits,
	   typename ArgType>
  class FOREACH {
  public:
    
    enum { NROW=Traits::NROW };

    static void EXEC(const MatType& lhs, ArgType& arg) {
      FOREACH_<FunctorType, MatType, Traits, ArgType, NROW>::EXEC(lhs, arg);
    }

    static void EXEC(MatType& lhs, ArgType& arg) {
      FOREACH_<FunctorType, MatType, Traits, ArgType, NROW>::EXEC(lhs, arg);
    }
  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_ForEach_H */
