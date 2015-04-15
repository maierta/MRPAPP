//-*-C++-*-

#ifndef PSIMAG_MatEqual_H
#define PSIMAG_MatEqual_H


/** \file MatEqual.h
 *
 * \brief Contains template functions to test the equality of two
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
#include "MatReduce.h"
#include "MatReduce2.h"
#include "NullType.h"

namespace psimag {


  //======================================================================
  //======================================================================
  /**
   * \brief The primary template class for the EQUAL<...>:EXEC
   *        function, which compares the equality of two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  EQUAL<MatType, MatTypeO, Traits, TriatsO>::EXEC(m,other)
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
  class EQUAL {
  public:

    typedef NullType ArgType;
    
    enum { NROW = Traits::NROW, 
	   NCOL = Traits::NCOL };
    
    class EqualReducerType {
    public:

      typedef bool T;
      typedef typename Traits::ElType  EType;
      typedef typename TraitsO::ElType OType;

      static T COMBINE(T lhs, T rhs, ArgType& arg) {
	return lhs && rhs;
      }

      static T PRODUCE(EType el, OType otherEl, ArgType& arg) {
	return el == otherEl;
      }
    };

    static bool EXEC(const MatType& m, const MatTypeO& other) { 
      
      ArgType ARG;
      return REDUCE2<EqualReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG);
    }
  };

  //======================================================================
  // 

  template<typename MatType, 
	   typename MatTypeO, 
	   typename Traits,
	   typename TraitsO,
	   typename Algorithms>
  class CLOSE {
  public:

    typedef NullType ArgType;
    
    enum { NROW = Traits::NROW, 
	   NCOL = Traits::NCOL };
    
    class CloseReducerType {
    public:

      typedef bool T;
      typedef typename Traits::ElType  EType;
      typedef typename TraitsO::ElType OType;

      static T COMBINE(T lhs, T rhs, ArgType& arg) {
	return lhs && rhs;
      }

      static T PRODUCE(EType el, OType otherEl, ArgType& arg) {
	return Algorithms::close(el,otherEl);
      }
    };

    static bool EXEC(const MatType& m, const MatTypeO& other) { 
      
      ArgType ARG;
      return REDUCE2<CloseReducerType, MatType, MatTypeO, Traits, TraitsO, ArgType>::EXEC(m,other, ARG);
    }
  };

  //======================================================================
  // 

  template<typename T, 
	   template<typename, size_t> class MatType,
	   template<typename, size_t> class MatTypeO,
	   size_t NROW,
	   typename Algorithms>
  bool Close (const MatType<T,NROW>&  m, 
	      const MatTypeO<T,NROW>& other) {
    return CLOSE<MatType <T,NROW>, 
                 MatTypeO<T,NROW>, 
                 typename MatType<T,NROW>::Traits,  
      typename MatTypeO<T,NROW>::Traits, Algorithms >::EXEC(m,other); 
  }
 
  //======================================================================
  //======================================================================
  //======================================================================

  //======================================================================
  /**
   * \brief An == operator comparing matrix types of the from,
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
  template<typename T, size_t NROW, size_t NCOL,
	   template<typename, size_t, size_t>           class TraitsTemplate,
	   template<typename, size_t, size_t>           class TraitsTemplateO,
	   template<typename, size_t, size_t, typename> class MatType,
	   template<typename, size_t, size_t, typename> class MatTypeO>
  bool operator == (const MatType <T,NCOL,NROW, TraitsTemplate <T,NCOL,NROW> >& m, 
		    const MatTypeO<T,NCOL,NROW, TraitsTemplateO<T,NCOL,NROW> >& other) {
    return EQUAL<MatType<T,NCOL,NROW, TraitsTemplate<T,NCOL,NROW> >, 
      MatTypeO<T,NCOL,NROW, TraitsTemplateO<T,NCOL,NROW> >, 
      TraitsTemplate<T,NCOL,NROW>,
      TraitsTemplateO<T,NCOL,NROW> >::EXEC(m,other); 
  }


//======================================================================
/**
 * \brief An preprocessor macro for generating == operators which
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

  //======================================================================
  //======================================================================
  /**
   * \brief The primary template class for the EQUAL<...>:EXEC
   *        function, which compares the equality of two marix-like
   *        objects (not necessarily of the same type).
   *
   * It is initially invoked like: 
   *
   *  EQUAL<MatType, MatTypeO, Traits, TriatsO>::EXEC(m,other)
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
	   typename Traits,
	   typename ArgType>
  class EQUAL_VAL {
  public:

    enum { NROW = Traits::NROW, 
	   NCOL = Traits::NCOL };
    
    class EqualReducerType {
    public:

      typedef bool T;
      typedef typename Traits::ElType  EType;

      static T COMBINE(T lhs, T rhs, ArgType& arg) {
	return lhs && rhs;
      }

      static T PRODUCE(EType el, ArgType& arg) {
	return el == arg;
      }
    };

    static bool EXEC(const MatType& m, const ArgType& arg) { 
      ArgType copy(arg);
      return REDUCE<EqualReducerType, MatType, Traits, ArgType>::EXEC(m, copy);
    }
  };

  //======================================================================
  //======================================================================
  //======================================================================

  //======================================================================
  /**
   * \brief An == operator comparing matrix types of the from,
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
  template<typename T, 
	   template<typename, size_t, size_t> class MatType,
	   size_t NROW, size_t NCOL >
  bool operator == (const MatType<T,NCOL,NROW>&  m, 
		    const T& other) {
    return EQUAL_VAL<MatType<T,NCOL,NROW>, 
      typename MatType<T,NCOL,NROW>::Traits,  
      T >::EXEC(m,other); 
  }

} /* namespace psimag */

#endif /* PSIMAG_Mat_Equal_H */
