//-*-C++-*-

#ifndef PSIMAG_MatCopy_H
#define PSIMAG_MatCopy_H


/** 
 * \file MatCopy.h
 *
 * \brief Contains template functions to copy two matrix objects.
 *
 * The matrix objects to be copied may be of different types. The
 * types are specified by the MatT`ype and MatTypeO template
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
#include "MatForEach2.h"
#include "FieldConvert.h"

namespace psimag {


  //======================================================================
  //======================================================================
  /**
   * \brief The primary template class for the COPY<...>:EXEC
   *        function, which copies two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is invoked like: 
   *
   *  COPY<MatType, MatTypeO, Traits, TriatsO>::EXEC(m,other)
   * 
   * Forms a single large aggregate operation by assembling assignment
   * expressions for each row from the bottom row upward and from the
   * left column to the right.
   *
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatTypeO: The type of the matrix used on the RHS.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param TraitO:   A type (e.g. MatTraits<NROW,NCOL>) describng the other matrix.
   *
   */
  template<typename MatTypeLHS, 
	   typename MatTypeRHS, 
	   typename TraitsLHS,
	   typename TraitsRHS>
  class COPY {
  public:

    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    class CopyFunctionType {
    public:

      typedef typename TraitsLHS::ElType LHSType;
      typedef typename TraitsRHS::ElType RHSType;

      static void EXEC(LHSType& lhs, const RHSType& rhs) {
	lhs = convert<LHSType>(rhs);
      }
    };
    class CopyPlusFunctionType {
    public:
      
      typedef typename TraitsLHS::ElType LHSType;
      typedef typename TraitsRHS::ElType RHSType;
      
      static void EXEC(LHSType& lhs, const RHSType& rhs) {
	lhs += convert<LHSType>(rhs);
      }
    };
    class CopyMinusFunctionType {
    public:
      
      typedef typename TraitsLHS::ElType LHSType;
      typedef typename TraitsRHS::ElType RHSType;
      
      static void EXEC(LHSType& lhs, const RHSType& rhs) {
	lhs -= convert<LHSType>(rhs);
      }
    };
    class CopyTimesFunctionType {
    public:
      
      typedef typename TraitsLHS::ElType LHSType;
      typedef typename TraitsRHS::ElType RHSType;
      
      static void EXEC(LHSType& lhs, const RHSType& rhs) {
	lhs *= convert<LHSType>(rhs);
      };
    };

    static void EXEC(MatTypeLHS& lhs, const MatTypeRHS& rhs) { 
      return FOREACH2<CopyFunctionType, MatTypeLHS, MatTypeRHS, TraitsLHS, TraitsRHS>::EXEC(lhs, rhs);
    }

    static void EXEC_PLUS(MatTypeLHS& lhs, const MatTypeRHS& rhs) { 
      return FOREACH2<CopyPlusFunctionType, MatTypeLHS, MatTypeRHS, TraitsLHS, TraitsRHS>::EXEC(lhs, rhs);
    }

    static void EXEC_MINUS(MatTypeLHS& lhs, const MatTypeRHS& rhs) { 
      return FOREACH2<CopyMinusFunctionType, MatTypeLHS, MatTypeRHS, TraitsLHS, TraitsRHS>::EXEC(lhs, rhs);
    }

    static void EXEC_TIMES(MatTypeLHS& lhs, const MatTypeRHS& rhs) { 
      return FOREACH2<CopyTimesFunctionType, MatTypeLHS, MatTypeRHS, TraitsLHS, TraitsRHS>::EXEC(lhs, rhs);
    }

  };

  //======================================================================
  //======================================================================
  //======================================================================

  //======================================================================
  /**
   * \brief An = operator for assigning matrix types of the from,
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
   * Since there is no partial specialization of operators in c++,
   * operators of this sort must be selected by matching the argument
   * types. Since the arguments are of the proscribed form the
   * compiler then knows NROW, NCOL, and the element type T. This
   * allows it to select the right EQUAL_ELS<...>:EXEC function.
   *
   * \param NROW:     The number of rows in the two matricies being compared.
   * \param NCOL:     The number of columns in the two matricies being compared.
   *
   */
    template<typename MatType, typename MatTypeO>
    void Copy (MatType&  m, MatTypeO& other) {
      return COPY<MatType, MatTypeO, typename MatType::Traits, typename MatTypeO::Traits>::EXEC(m,other);
    }


  //======================================================================
  /**
   * \brief An preprocessor macro for generating = operators which
   *        work with square arrays represented by C types like T[4]
   *        and T[9] etc.
   *
   * C type like T[9] are constructed by the ColMajorTraits<T,DIM> class.
   *
   * One would think that the template system could be used instead of
   * this macro but it proved to be problematic.
   *
   * \param T:     The type of the C array T[DIM*DIM]
   * \param DIM:   The dimension of the square array represented by the given C array.
   *
   * \note This is still under development. It is tricky since passing
   *       int (&) [9] as a function argument forgets the dimension
   *       resulting in int* (&).
   *
   */
  /*
   #define Mat_Equal_Template(T, DIM)					\
     bool  operator == (ColMajorTypeTraits<T,DIM, DIM>::ConstRefType m,	\
   		     ColMajorTypeTraits<T,DIM, DIM>::ConstRefType other) {	\
       return EQUAL<ColMajorTypeTraits<T,DIM, DIM>::Type,			\
         ColMajorTraits<T,DIM,DIM>::Type, ColMajorTraits<DIM, DIM>, DIM-1>::EXEC(m, other); }
  */
  /**
   * Note the following does not work. 
   */
  //   bool  operator == (ColMajorTypeTraits<int,3,3>::Type m, ColMajorTypeTraits<int,3,3>::Type other) {
  //     for(size_t i =0; i< 9; i++) 
  //       if (m[i] != other[i]) return false;
  //     return true;
  //   }
  /**  \brief Overload to equal a  2x2 double array. */
  //Mat_Equal_Template(double,2);
  
  /**  \brief Overload to equal a 3x3 determinant of a double array. */
  //Mat_Equal_Template(double,3);

  /**  \brief Overload to equal a 2x2 determinant of an int array. */
  //Mat_Equal_Template(int,2);

  /**  \brief Overload to equal a 3x3 determinant of an int array. */
  //Mat_Equal_Template(int,3);


} /* namespace psimag */

#endif /* PSIMAG_Mat_Equal_H */
