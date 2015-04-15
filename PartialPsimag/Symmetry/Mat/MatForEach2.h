//-*-C++-*-

#ifndef PSIMAG_MatForeach2_H
#define PSIMAG_MatForeach2_H


/** \file MatForeach2.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         Foreach2 operation. This operation applies a given function
 *         to every tuple of elements in a given set of matricies.
 *
 * The Foreach2 operation is constructed as a static method on a
 * templated Foreach2 class:
 *
 * e.g.  FOREACH2_<FunctorType, MatType, MatTypeRHS, NROW, NCOL, ROW>::EXEC(lhs, rhs)
 *
 * This static function is constructed at compile time through
 * recursive invokations of itself (and it's row == 0 specialization) and
 * through invokations of:
 *
 * FOREACH2_ELS<FunctorType, MatType, MatTypeRHS, NROW, NCOL, ROW, COL>::EXEC(lhs, rhs)
 *
 * and it's specializations. The specializations are used to halt the
 * recursive construction of the FOREACH2_<...>::EXEC function.
 *
 * Using this metaprogramming technique we construct loop-less matrix
 * foreach2 operation for any size matricies.
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
   *        FOREACH2_ELS<...>:EXEC function, which operates on the
   *        corresponding rows of matrix-like objects (not necessarily
   *        of the same type).
   *
   * It is initially invoked like: 
   *
   *  FOREACH2_ELS<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, ROW, NCOL-1>::EXEC(lhs, rhs)
   * 
   * This operation traverses the matrix in row-major order (for the sake of printing)>
   *
   * \param FunctorType:  The type which provides the operation through it's EXEC statis methos.
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatType): The type of the matrix used on the RHS.
   * \param TraitsLHS:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param TraitsRHS:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param ROW:      The current row (zero indexed) to assemble the boolean expression for.
   * \param COL:      The current row (zero indexed) to assemble the boolean expression for.
   *
   */
  template<typename FunctorType, 
	   typename MatType, 
	   typename MatTypeRHS, 
	   typename TraitsLHS,
	   typename TraitsRHS,
	   size_t ROW, size_t COLS_LEFT>
  class FOREACH2_ELS {
  public:

    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    static const size_t COL = NCOL - COLS_LEFT;

    typedef typename TraitsLHS::template REF<MatType,   ROW,COL> REF_LHS;
    typedef typename TraitsRHS::template REF<MatTypeRHS,ROW,COL> REF_RHS;
    

    static void EXEC(MatType& lhs, const MatTypeRHS& rhs) {

      FunctorType::EXEC(REF_LHS::GET(lhs), REF_RHS::GETCONST(rhs));

      FOREACH2_ELS<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, ROW, COLS_LEFT-1>::EXEC(lhs, rhs);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH2_ELS<...>:EXEC.
   *
   * The FOREACH2_ELS<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename FunctorType, 
	   typename MatType, 
	   typename MatTypeRHS, 
	   typename TraitsLHS,
	   typename TraitsRHS,
	   size_t ROW>
  class FOREACH2_ELS<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, ROW, 1> {
  public:
    
    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    static const size_t COL = NCOL - 1;

    typedef typename TraitsLHS::template REF<MatType,   ROW,COL> REF_LHS;
    typedef typename TraitsRHS::template REF<MatTypeRHS,ROW,COL> REF_RHS;
    
    static void EXEC(MatType& lhs, const MatTypeRHS& rhs) {

      FunctorType::EXEC(REF_LHS::GET(lhs), REF_RHS::GETCONST(rhs));
    }
  };
  
  //======================================================================
  /**
   * \brief The primary template class for assembling the
   *        FOREACH2_<...>:EXEC function, which operates on
   *        two marix-like objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  FOREACH2_<FunctorType,MatType, MatTypeRHS, NROW, NCOL, NROW-1>::EXEC(lhs, rhs)
   * 
   * It recursively forms a single operation by
   * assembling operations for each row from the bottom row
   * upward. The recursion is terminated by a specialization of
   * FOREACH2_<...>:EXEC for row zero.
   *
   * \param FunctorType:  The type which provides the operation through it's EXEC statis methos.
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatTypeO: The type of the matrix used on the RHS.
   * \param TraitsLHS:   A type (e.g. MatTraits<NROW,NCOL>) describing the matrix.
   * \param TraitsRHS:  A type (e.g. MatTraits<NROW,NCOL>) describing the rhs matrix.
   * \param ROW:      The current row (zero indexed).
   *
   */
  template<typename FunctorType,
	   typename MatType, 
	   typename MatTypeRHS, 
	   typename TraitsLHS,
	   typename TraitsRHS,
	   size_t ROWS_LEFT>
  class FOREACH2_ {
  public:

    enum { NROW=TraitsLHS::NROW, 
	   NCOL=TraitsLHS::NCOL };

    static const size_t ROW = NROW - ROWS_LEFT;
    
    static void EXEC(MatType& lhs, const MatTypeRHS& rhs) { 

      FOREACH2_ELS<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, ROW, NCOL>::EXEC(lhs, rhs);
      FOREACH2_<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, ROWS_LEFT-1>::EXEC(lhs, rhs);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH2_<...>:EXEC.
   *
   * The FOREACH2_<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename FunctorType,
	   typename MatType, 
	   typename MatTypeRHS, 
	   typename TraitsLHS,
	   typename TraitsRHS>
  class FOREACH2_<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, 1> {
  public:
    
    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };

    static const size_t ROW = NROW - 1;

    static void EXEC(MatType& lhs, const MatTypeRHS& rhs) {

      FOREACH2_ELS<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, ROW, NCOL>::EXEC(lhs, rhs);
    }
  };

  //======================================================================
  /**
   * \brief The FOREACH2 template that calls FOREACH2_ to start the
   *        recursive computation of the aggregate operation at
   *        ROWS_LEFT =  NROW.
   *
   */
  template<typename FunctorType,
	   typename MatType, 
	   typename MatTypeRHS, 
	   typename TraitsLHS,
	   typename TraitsRHS>
  class FOREACH2 {
  public:
    
    enum { NROW=TraitsLHS::NROW };

    static void EXEC(MatType& lhs, const MatTypeRHS& rhs) {
      FOREACH2_<FunctorType, MatType, MatTypeRHS, TraitsLHS, TraitsRHS, NROW>::EXEC(lhs, rhs);
    }
  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_Foreach2_H */
