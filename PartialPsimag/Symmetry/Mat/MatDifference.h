//-*-C++-*-

#ifndef PSIMAG_MatDifference_H
#define PSIMAG_MatDifference_H


/** 
 * \file MatDifference.h
 *
 * \brief Contains template functions to difference two matrix objects.
 *
 * The matrix objects to be copied may be of different types. 
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
#include "MatForEach3.h"
#include "FieldConvert.h"
namespace psimag {


  //======================================================================
  //======================================================================
  /**
   * \brief The template class for the DIFFERENCE<...>:EXEC
   *        function, which computes the element-wise difference of two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is invoked like: 
   *
   *  DIFFERENCE<MatType, MatTypeO, MAtTypeR, Traits, TriatsO, TraitsR>::EXEC(m,other)
   * 
   * Forms a single large aggregate operation by assembling assignment
   * expressions for each row from the bottom row upward and from the
   * left column to the right.
   *
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatTypeO: The type of the matrix used on the RHS.
   * \param MatTypeO: The type of the resulting matrix.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param TraitO:   A type (e.g. MatTraits<NROW,NCOL>) describng the other matrix.
   * \param TraitResult:   A type (e.g. MatTraits<NROW,NCOL>) describng the resulting matrix.
   *
   */
  template<typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename MatTypeResult, 
	   typename TraitsLHS,
	   typename TraitsRHS,
	   typename TraitsResult>
  class DIFFERENCE {
  public:

    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    class DifferenceFunctionType {
    public:

      typedef typename TraitsLHS::ElType    LHSType;
      typedef typename TraitsRHS::ElType    RHSType;
      typedef typename TraitsResult::ElType ResultType;

      static void EXEC(const LHSType& lhs, const RHSType& rhs, ResultType& result) {
	result = convert<ResultType>(lhs) - convert<ResultType>(rhs);
      }
    };

    static void EXEC(const MatTypeLHS& lhs, const MatTypeRHS& rhs, MatTypeResult& result) { 
      FOREACH3<DifferenceFunctionType, MatTypeResult, MatTypeLHS, MatTypeRHS, TraitsResult, TraitsLHS, TraitsRHS>::EXEC(lhs, rhs, result);
    }
  };

  //======================================================================
  //======================================================================
  //======================================================================

  //======================================================================
  /**
   * \brief An - for matrix types of the form,
   *        template<T,NROW,NCOL> MatType.
   *
   * \param MatType: The <i> template template </i> parameter used to
   *        specify the LHS matrix of the form template<T,NROW,NCOL>
   *        MatType.
   *
   * \param MatTypeO: The <i> template template </i> parameter used to
   *        specify the RHS matrix of the form template<T,NROW,NCOL>
   *        MatTypeO.
   *
   */
  template<typename T, size_t NROW, size_t NCOL ,
	   template<typename, size_t, size_t> class TraitsTemplate,
 	   template<typename, size_t, size_t, typename> class MatType,
	   template<typename, size_t, size_t> class TraitsTemplateO,
 	   template<typename, size_t, size_t, typename> class MatTypeO>

  MatType<T,NCOL,NROW, TraitsTemplate<T,NCOL,NROW> > 
  operator - (const MatType <T,NCOL,NROW, TraitsTemplate <T,NCOL,NROW> >&  m, 
	      const MatTypeO<T,NCOL,NROW, TraitsTemplateO<T,NCOL,NROW> >& other) {

    MatType<T,NCOL,NROW,TraitsTemplate<T,NCOL,NROW> > result;

    DIFFERENCE<MatType<T,NROW,NCOL,TraitsTemplate<T,NCOL,NROW> >, 
      MatTypeO<T,NROW,NCOL,TraitsTemplateO<T,NCOL,NROW> >, 
      MatType<T,NROW,NCOL,TraitsTemplate<T,NCOL,NROW> >,
      TraitsTemplate<T,NCOL,NROW>,
      TraitsTemplateO<T,NCOL,NROW>,
      TraitsTemplate<T,NCOL,NROW> >::EXEC(m, other, result);
    return result;
  }

} /* namespace psimag */

#endif /* PSIMAG_Mat_Equal_H */
