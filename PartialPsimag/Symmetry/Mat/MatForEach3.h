//-*-C++-*-

#ifndef PSIMAG_MatForEach3_H
#define PSIMAG_MatForEach3_H


/** \file MatForEach3.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         ForEach3 operation. This operation applies a given function
 *         to every tuple of elements in a given set of matricies.
 *
 * The matrix objects to be operated on may be of different types. The
 * types are specified by template arguments like MatTypeLHS and
 * MatTypeRHS. These types need only provide a [] operator and be of
 * the same shape.
 *
 * \param The order of elements returned by operator[] for each
 *        MatType is given by a corresponding template argument
 *        Traits.
 *
 * The ForEach3 operation is constructed as a static method on a
 * templated ForEach3 class:
 *
 * e.g.  FOREACH3_<FunctorType, MatTypeLHS, MatTypeRHS, NROW, NCOL, ROW>::EXEC(lhs, rhs)
 *
 * This static function is constructed at compile time through
 * recursive invokations of itself (and it's row == 0 specialization) and
 * through invokations of:
 *
 * FOREACH3_ELS<FunctorType, MatTypeLHS, MatTypeRHS, NROW, NCOL, ROW, COL>::EXEC(lhs, rhs)
 *
 * and it's specializations. The specializations are used to halt the
 * recursive construction of the FOREACH3_<...>::EXEC function.
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
   *        FOREACH3_ELS<...>:EXEC function, which operates on the
   *        corresponding rows of matrix-like objects (not necessarily
   *        of the same type).
   *
   * It is initially invoked like: 
   *
   *  FOREACH3_ELS<FunctorType, MatTypeLHS, MatTypeRHS, TraitsLHS, TraitsRHS, ROW, NCOL-1>::EXEC(lhs, rhs)
   * 
   * It recursively forms a single combined operation for the given ROW
   * by assembling operations for each matrix position in the
   * row from the left column toward the right. The recursion is
   * terminated by a specialization of FOREACH3_ELS<...>:EXEC for column
   * zero.
   *
   * \param FunctorType:  The type which provides the operation through it's EXEC static method.
   * \param MatTypeResult:      The type of the resulting matrix.
   * \param MatTypeLHS:      The type of the matrix used on the LHS.
   * \param MatTypeRHS:   The type of the matrix used on the RHS.
   * \param TraitsResult:    A type (e.g. MatTraits<NROW,NCOL>) describng the resulting matrix.
   * \param TraitsLHS:    A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param TraitsRHS:    A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param ROW:         The current row (zero indexed) to assemble the boolean expression for.
   * \param COL:         The current row (zero indexed) to assemble the boolean expression for.
   *
   */
  template<typename FunctorType, 
	   typename MatTypeResult, 
	   typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename TraitsResult,
	   typename TraitsLHS,
	   typename TraitsRHS,
	   size_t ROW, size_t COL>
  class FOREACH3_ELS {
  public:

    typedef typename TraitsLHS::template REF<MatTypeResult,ROW,COL> REF_Result;
    typedef typename TraitsLHS::template REF<MatTypeLHS,   ROW,COL> REF_LHS;
    typedef typename TraitsRHS::template REF<MatTypeRHS,   ROW,COL> REF_RHS;
    
    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    static void EXEC(const MatTypeLHS& lhs, const MatTypeRHS& rhs, MatTypeResult& result) {

      FunctorType::EXEC(REF_LHS::GETCONST(lhs), REF_RHS::GETCONST(rhs), REF_Result::GET(result));

      FOREACH3_ELS<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, ROW, COL-1>::EXEC(lhs, rhs, result);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH3_ELS<...>:EXEC.
   *
   * The FOREACH3_ELS<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename FunctorType, 
	   typename MatTypeResult, 
	   typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename TraitsResult,
	   typename TraitsLHS,
	   typename TraitsRHS,
	   size_t ROW>
  class FOREACH3_ELS<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, ROW, 0> {
  public:
    
    typedef typename TraitsLHS::template REF<MatTypeResult,ROW,0>   REF_Result;
    typedef typename TraitsLHS::template REF<MatTypeLHS,   ROW,0>   REF_LHS;
    typedef typename TraitsRHS::template REF<MatTypeRHS,   ROW,0>   REF_RHS;
    
    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    static void EXEC(const MatTypeLHS& lhs, const MatTypeRHS& rhs, MatTypeResult& result) {

      FunctorType::EXEC(REF_LHS::GETCONST(lhs), REF_RHS::GETCONST(rhs), REF_Result::GET(result));
    }
  };
  
  //======================================================================
  /**
   * \brief The primary template class for assembling the
   *        FOREACH3_<...>:EXEC function, which operates on
   *        two marix-like objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  FOREACH3_<FunctorType,MatTypeLHS, MatTypeRHS, NROW, NCOL, NROW-1>::EXEC(lhs, rhs)
   * 
   * It recursively forms a single operation by
   * assembling operations for each row from the bottom row
   * upward. The recursion is terminated by a specialization of
   * FOREACH3_<...>:EXEC for row zero.
   *
   * \param FunctorType:  The type which provides the operation through it's EXEC statis methos.
   * \param MatTypeLHS:   The type of the matrix used on the LHS.
   * \param MatTypeRHS:   The type of the matrix used on the RHS.
   * \param TraitsLHS:    A type (e.g. MatTraits<NROW,NCOL>) describing the matrix.
   * \param TraitsRHS:    A type (e.g. MatTraits<NROW,NCOL>) describing the rhs matrix.
   * \param ROW:         The current row (zero indexed).
   *
   */
  template<typename FunctorType,
	   typename MatTypeResult, 
	   typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename TraitsResult,
	   typename TraitsLHS,
	   typename TraitsRHS,
	   size_t ROW>
  class FOREACH3_ {
  public:

    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    static void EXEC(const MatTypeLHS& lhs, const MatTypeRHS& rhs, MatTypeResult& result) { 

      FOREACH3_ELS<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, ROW, NCOL-1>::EXEC(lhs, rhs, result);
      FOREACH3_<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, ROW-1>::EXEC(lhs, rhs, result);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        FOREACH3_<...>:EXEC.
   *
   * The FOREACH3_<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename FunctorType,
	   typename MatTypeResult, 
	   typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename TraitsResult,
	   typename TraitsLHS,
	   typename TraitsRHS>
  class FOREACH3_<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, 0> {
  public:
    
    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };

    static void EXEC(const MatTypeLHS& lhs, const MatTypeRHS& rhs, MatTypeResult& result) {

      FOREACH3_ELS<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, 0, NCOL-1>::EXEC(lhs, rhs, result);
    }
  };

  //======================================================================
  /**
   * \brief The FOREACH3 template that calls FOREACH3_ to start the
   *        recursive computation of the aggregate operation at
   *        row NROW-1.
   *
   */
  template<typename FunctorType,
	   typename MatTypeResult, 
	   typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename TraitsResult,
	   typename TraitsLHS,
	   typename TraitsRHS>
  class FOREACH3 {
  public:
    
    enum { NROW=TraitsLHS::NROW };

    static void EXEC(const MatTypeLHS& lhs, const MatTypeRHS& rhs, MatTypeResult& result) {
      FOREACH3_<FunctorType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS, NROW-1>::EXEC(lhs, rhs, result);
    }
  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_ForEach3_H */
