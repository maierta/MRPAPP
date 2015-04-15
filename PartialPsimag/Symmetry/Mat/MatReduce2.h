//-*-C++-*-

#ifndef PSIMAG_MatReduce2_H
#define PSIMAG_MatReduce2_H


/** \file MatReduce2.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         Reduce operation. This operation applies a given function
 *         to every tuple of elements in a given set of matricies. It
 *         then collects the results into a single value using another
 *         function.
 *
 * The matrix objects to be operated on may be of different types. The
 * types are specified by template arguments like MatType and
 * MatTypeOther. These types need only provide a [] operator and be of
 * the same shape.
 *
 * \param The order of elements returned by operator[] for each
 *        MatType is given by a corresponding template argument
 *        Traits.
 *
 * The Reduce operation is constructed as a static method on a
 * templated Reduce class:
 *
 * e.g.  REDUCE2_<ReducerType, MatType, MatTypeOther, NROW,
 *              NCOL, ROW>::EXEC(m, other)
 *
 * This static function is constructed at compile time through
 * recursive invokations of itself (and it's row == 0 specialization) and
 * through invokations of:
 *
 * REDUCE2_ELS<ReducerType, MatType, MatTypeOther, NROW, NCOL, ROW, COL>::EXEC(m, other)
 *
 * and it's specializations. The specializations are used to halt the
 * recursive construction of the REDUCE2_<...>::EXEC function.
 *
 * Using this metaprogramming technique we construct loop-less matrix
 * reduce operation for any size matricies.
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
   *        REDUCE2_ELS<...>:EXEC function, which operates on the
   *        corresponding rows of matrix-like objects (not necessarily
   *        of the same type).
   *
   * It is initially invoked like: 
   *
   *  REDUCE2_ELS<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, ROW, NCOL-1>::EXEC(m,other)
   * 
   * It recursively forms a single combined operation for the given ROW
   * by assembling operations for each matrix position in the
   * row from the left column toward the right. The recursion is
   * terminated by a specialization of REDUCE2_ELS<...>:EXEC for column
   * zero.
   *
   * \param ReducerType:  The type which provides the reducer which provides a return type and two functions.
   * \param MatType:      The type of the matrix used on the LHS.
   * \param MatTypeO:     The type of the matrix used on the RHS.
   * \param Traits:       A type (e.g. MatTraits<NROW,NCOL>) describing the LHS matrix.
   * \param TraitsO:      A type (e.g. MatTraits<NROW,NCOL>) describing the RHS matrix.
   * \param ArgType:      The type of the auxillary argument.
   * \param ROW:          The current row (zero indexed) to assemble the boolean expression for.
   * \param COL:          The current row (zero indexed) to assemble the boolean expression for.
   *
   */
  template<typename Reducer, 
	   typename MatType, 
	   typename MatTypeOther, 
	   typename Traits,
	   typename TraitsO,
	   typename ArgType,
	   size_t ROW, size_t COL>
  class REDUCE2_ELS {
  public:
    
    typedef typename Reducer::T ReturnType; 
    typedef typename Traits::template  REF<MatType,ROW,COL> REF_LHS;
    typedef typename TraitsO::template REF<MatTypeOther,ROW,COL> REF_RHS;

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    static ReturnType EXEC(const MatType& m, const MatTypeOther& other, ArgType& arg) {

      return Reducer::COMBINE(Reducer::PRODUCE(REF_LHS::GETCONST(m),REF_RHS::GETCONST(other), arg),
			      REDUCE2_ELS<Reducer, MatType, MatTypeOther, Traits, TraitsO, ArgType, ROW, COL-1>::EXEC(m, other, arg),
			      arg);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        REDUCE2_ELS<...>:EXEC.
   *
   * The REDUCE2_ELS<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename ReducerType, 
	   typename MatType, 
	   typename MatTypeOther, 
	   typename Traits,
	   typename TraitsO,
	   typename ArgType,
	   size_t ROW>
  class REDUCE2_ELS<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, ROW, 0> {
  public:

    typedef typename ReducerType::T ReturnType;
    typedef typename Traits::template  REF<MatType,ROW,0> REF_LHS;
    typedef typename TraitsO::template REF<MatTypeOther,ROW,0> REF_RHS;
  
    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    static ReturnType EXEC(const MatType& m, const MatTypeOther& other, ArgType& arg) {

      return ReducerType::PRODUCE(REF_LHS::GETCONST(m), REF_RHS::GETCONST(other), arg);
    }
  };
  
  //======================================================================
  /**
   * \brief The primary template class for assembling the
   *        REDUCE2_<...>:EXEC function, which operates on
   *        two marix-like objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  REDUCE2_<ReducerType,MatType, MatTypeOther, NROW, NCOL, NROW-1>::EXEC(m,other)
   * 
   * It recursively forms a single operation by
   * assembling operations for each row from the bottom row
   * upward. The recursion is terminated by a specialization of
   * REDUCE2_<...>:EXEC for row zero.
   *
   * \param ReducerType:  The type which provides the return type and the combine and procuce operations.
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatTypeO: The type of the matrix used on the RHS.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describing the matrix.
   * \param TraitsO:  A type (e.g. MatTraits<NROW,NCOL>) describing the other matrix.
   * \param ArgType:  The type of the auxillary argument.
   * \param ROW:      The current row (zero indexed).
   *
   */
  template<typename ReducerType, 
	   typename MatType, 
	   typename MatTypeOther, 
	   typename Traits,
	   typename TraitsO,
	   typename ArgType,
	   size_t ROW>
  class REDUCE2_ {
  public:

    typedef typename ReducerType::T ReturnType; 

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    static ReturnType EXEC(const MatType& m, const MatTypeOther& other, ArgType& arg) { 

      return
	ReducerType::COMBINE(REDUCE2_ELS<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, ROW, NCOL-1>::EXEC(m,other, arg),
			     REDUCE2_<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, ROW-1>::EXEC(m,other, arg),
			     arg);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        REDUCE2<...>:EXEC.
   *
   * The REDUCE2<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename ReducerType,
	   typename MatType, 
	   typename MatTypeOther, 
	   typename Traits,
	   typename TraitsO,
	   typename ArgType>
  class REDUCE2_<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, 0> {
  public:
    
    typedef typename ReducerType::T ReturnType; 

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };

    static ReturnType EXEC(const MatType& m, const MatTypeOther& other, ArgType& arg) {
      return REDUCE2_ELS<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, 0, NCOL-1>::EXEC(m, other, arg);
    }
  };

  //======================================================================
  /**
   * \brief The Reduce2 template that calls REDUCE2_ to start the
   *        recursive computation of the aggregate reduce operation at
   *        row NROW-1.
   *
   */
  template<typename ReducerType,
	   typename MatType, 
	   typename MatTypeOther, 
	   typename Traits,
	   typename TraitsO,
	   typename ArgType>
  class REDUCE2 {
  public:
    
    typedef typename ReducerType::T ReturnType; 

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };

    static ReturnType EXEC(const MatType& m, const MatTypeOther& other, ArgType& arg) {
      return REDUCE2_<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, NROW-1>::EXEC(m,other, arg);
    }
  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_Reduce2_H */
