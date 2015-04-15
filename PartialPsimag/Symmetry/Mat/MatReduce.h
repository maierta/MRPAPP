//-*-C++-*-

#ifndef PSIMAG_MatReduce_H
#define PSIMAG_MatReduce_H


/** \file MatReduce.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         Reduce operation. This operation applies a given function
 *         elements of a given matrix and a given argument. It then
 *         collects the results into a single value using another
 *         function.
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
   * \brief The primary template typename for assembling the
   *        REDUCE_ELS<...>:EXEC function, which operates on the
   *        corresponding rows of matrix-like objects (not necessarily
   *        of the same type).
   *
   * It is initially invoked like: 
   *
   *  REDUCE_ELS<ReducerType, MatType, MatTypeOther, Traits, TraitsO, ArgType, ROW, NCOL-1>::EXEC(m,other)
   * 
   * It recursively forms a single combined operation for the given ROW
   * by assembling operations for each matrix position in the
   * row from the left column toward the right. The recursion is
   * terminated by a specialization of REDUCE_ELS<...>:EXEC for column
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
	   typename Traits,
	   typename ArgType,
	   size_t ROW, size_t COL>
  class REDUCE_ELS {
  public:
    
    typedef typename Reducer::T ReturnType; 
    typedef typename Traits::template  REF<MatType,ROW,COL> REF_LHS;

    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    static ReturnType EXEC(const MatType& m, ArgType& arg) {

      return Reducer::COMBINE(Reducer::PRODUCE(REF_LHS::GETCONST(m), arg),
			      REDUCE_ELS<Reducer, MatType, Traits, ArgType, ROW, COL-1>::EXEC(m,  arg),
			      arg);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        REDUCE_ELS<...>:EXEC.
   *
   * The REDUCE_ELS<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename ReducerType, 
	   typename MatType, 
	   typename Traits,
	   typename ArgType,
	   size_t ROW>
  class REDUCE_ELS<ReducerType, MatType, Traits, ArgType, ROW, 0> {
  public:

    typedef typename ReducerType::T ReturnType;
    typedef typename Traits::template  REF<MatType,ROW,0> REF_LHS;
  
    enum { NROW=Traits::NROW, NCOL=Traits::NCOL };
    
    static ReturnType EXEC(const MatType& m, ArgType& arg) {

      return ReducerType::PRODUCE(REF_LHS::GETCONST(m), arg);
    }
  };
  
  //======================================================================
  /**
   * \brief The primary template class for assembling the
   *        REDUCE_<...>:EXEC function, which operates on
   *        two marix-like objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  REDUCE_<ReducerType,MatType, MatTypeOther, NROW, NCOL, NROW-1>::EXEC(m,other)
   * 
   * It recursively forms a single operation by
   * assembling operations for each row from the bottom row
   * upward. The recursion is terminated by a specialization of
   * REDUCE_<...>:EXEC for row zero.
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
	   typename Traits,
	   typename ArgType,
	   size_t ROW>
  class REDUCE_ {
  public:

    typedef typename ReducerType::T ReturnType; 

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };
    
    static ReturnType EXEC(const MatType& m, ArgType& arg) { 

      return
	ReducerType::COMBINE(REDUCE_ELS<ReducerType, MatType, Traits, ArgType, ROW, NCOL-1>::EXEC(m, arg),
			     REDUCE_<ReducerType, MatType, Traits, ArgType, ROW-1>::EXEC(m, arg),
			     arg);
    }
  };

  //======================================================================
  /**
   * \brief The row zero specialization of the primary template class,
   *        REDUCE<...>:EXEC.
   *
   * The REDUCE<..,0>:EXEC function does the same thing as it's
   * counterpart in the primary template except that is does not
   * recurse.
   *
   */
  template<typename ReducerType,
	   typename MatType, 
	   typename Traits,
	   typename ArgType>
  class REDUCE_<ReducerType, MatType, Traits, ArgType, 0> {
  public:
    
    typedef typename ReducerType::T ReturnType; 

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };

    static ReturnType EXEC(const MatType& m, ArgType& arg) {
      return REDUCE_ELS<ReducerType, MatType, Traits, ArgType, 0, NCOL-1>::EXEC(m,  arg);
    }
  };

  //======================================================================
  /**
   * \brief The Reduce template that calls REDUCE_ to start the
   *        recursive computation of the aggregate reduce operation at
   *        row NROW-1.
   *
   */
  template<typename ReducerType,
	   typename MatType, 
	   typename Traits,
	   typename ArgType>
  class REDUCE {
  public:
    
    typedef typename ReducerType::T ReturnType; 

    enum { NROW=Traits::NROW, 
	   NCOL=Traits::NCOL };

    static ReturnType EXEC(const MatType& m, ArgType& arg) {
      return REDUCE_<ReducerType, MatType, Traits, ArgType, NROW-1>::EXEC(m, arg);
    }
  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_Reduce_H */
