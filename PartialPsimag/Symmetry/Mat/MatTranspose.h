//-*-C++-*-

#ifndef PSIMAG_MatTranspose_H
#define PSIMAG_MatTranspose_H

/** \file MatTranspose.h
 *
 *  \brief Contains template functions providing a transpose operation
 *  on Matrices of any sort.
 *
 */
 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"

namespace psimag {

  /**
   *
   * \note ROW nd COL are zero indexed
   */
  template<typename MatType, 
	   typename MatTypeR, 
	   typename Traits, 
	   typename TraitsR, 
	   size_t ROW, 
	   size_t COL>
  class TRANSPOSE_ELS {
  public:
    
    typedef MatType              M;
    typedef MatTypeR             MR;
    typedef Traits               MT;
    typedef TraitsR              MRT;
    typedef typename MT::ElType  T;

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
      static T& from(MR& im) {
	return MRT::template REF<MR,R,C>::GET(im);
      }
    };
    
    static void EXEC(const MatType& m, MatTypeR& result) {
      get<ROW,COL>::from(result) = get<COL,ROW>::from(m);
      TRANSPOSE_ELS<MatType, MatTypeR, Traits, TraitsR, ROW, COL-1>::EXEC(m, result);
    }
  };

  //======================================================================

  template<typename MatType, 
	   typename MatTypeR, 
	   typename Traits, 
	   typename TraitsR, 
	   size_t ROW>
  class TRANSPOSE_ELS<MatType, MatTypeR, Traits, TraitsR, ROW, 0> {
  public:
    
    typedef MatType              M;
    typedef MatTypeR             MR;
    typedef Traits               MT;
    typedef TraitsR              MRT;
    typedef typename MT::ElType  T;

    enum { COL = 0 };

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
      static T& from(MR& im) {
	return MRT::template REF<MR,R,C>::GET(im);
      }
    };
    static void EXEC(const MatType& m, MatType& result) {
      get<ROW,COL>::from(result) = get<COL,ROW>::from(m);
    }
  };
  
  //======================================================================

  template<typename MatType, 
	   typename MatTypeR, 
	   typename Traits, 
	   typename TraitsR, 
	   size_t ROW>
  class TRANSPOSE {
  public:

    enum { NCOL = Traits::NCOL };

    static void EXEC(const MatType& m, MatTypeR& result) {
      TRANSPOSE_ELS<MatType, MatTypeR, Traits, TraitsR, ROW, NCOL-1>::EXEC(m,result);
      TRANSPOSE<MatType, MatTypeR, Traits, TraitsR, ROW-1>::EXEC(m,result);
    }
  };

  //======================================================================

  template<typename MatType, 
	   typename MatTypeR, 	   
	   typename Traits, 
	   typename TraitsR>
  class TRANSPOSE<MatType, MatTypeR, Traits, TraitsR, 0> {
  public:

    enum { ROW=0, NCOL = Traits::NCOL };
    
    static void EXEC(const MatType& m, MatType& result) {
      TRANSPOSE_ELS<MatType, MatTypeR, Traits, TraitsR, ROW, NCOL-1>::EXEC(m,result);
    }
  };

  //======================================================================
  //======================================================================
  /**
   * \brief Calculates the transpose of an object of type MatType.
   */
  template<typename T, size_t DIM,
	   template<typename, size_t, size_t>           class TraitsTemplate, 
	   template<typename, size_t, size_t, typename> class MatType, 
	   template<typename, size_t, size_t>           class TraitsTemplateR, 
	   template<typename, size_t, size_t, typename> class MatTypeR>

  inline  void Transpose(const MatType <T,DIM,DIM, TraitsTemplate <T,DIM,DIM> >& m, 
			       MatTypeR<T,DIM,DIM, TraitsTemplateR<T,DIM,DIM> >& result) { 

    TRANSPOSE<MatType <T,DIM,DIM, TraitsTemplate <T,DIM,DIM> >, 
      MatTypeR<T,DIM,DIM, TraitsTemplateR <T,DIM,DIM> >, 
      TraitsTemplate <T,DIM,DIM>,
      TraitsTemplateR<T,DIM,DIM>,
      DIM-1>::EXEC(m, result); 
  }

  //======================================================================

  /**  \brief Overload to calculate the DIMxDIM determinant of a T array. */
  //  template<typename T, size_t DIM>
#define TRANSPOSE_ARRAY(T, DIM)				\
 inline  void Transpose(const ColMajorTraits<T,DIM>::Type& m,	\
		 ColMajorTraits<T,DIM>::Type& result)	\
  {							\
    TRANSPOSE<ColMajorTraits<T,DIM>::Type,		\
      ColMajorTraits<T,DIM>::Type,			\
      ColMajorTraits<T,DIM>,				\
      ColMajorTraits<T,DIM>,				\
      DIM-1>::EXEC(m, result);   }
  
  TRANSPOSE_ARRAY(int,2)
  TRANSPOSE_ARRAY(int,3)
  TRANSPOSE_ARRAY(double,2)
  TRANSPOSE_ARRAY(double,3)

} /* namespace psimag */

#endif /* PSIMAG_MatTranspose_H */
