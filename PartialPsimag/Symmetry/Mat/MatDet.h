//-*-C++-*-

#ifndef PSIMAG_MatDet_H
#define PSIMAG_MatDet_H

/** \file MatDet.h
 *
 *  \brief Contains template functions providing det for
 *         standard operations on Matrices of any sort.
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

  template<typename MatType,
	   typename Traits,
	   size_t DIM> 
  class DET 
  {
  public:
    typedef typename Traits::ElType T;
    
    static T of(const MatType& m) {
      ASSERT(false,
	     std::exception("Algorithm not implemented yet."));
    }
  };

  //======================================================================

  template<typename MatType,
	   typename Traits>
  class DET<MatType, Traits, 1> {
  public:

    typedef MatType              M;
    typedef Traits               MT;
    typedef typename MT::ElType   T;

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
    };

    static T of(const MatType& m) {
      return get<0,0>::GETCONST(m);
    }
  };

  //======================================================================

  template<typename MatType,
	   typename Traits>
  class DET<MatType, Traits, 2> {
  public:
    
    typedef MatType              M;
    typedef Traits               MT;
    typedef typename MT::ElType   T;

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
    };

    static T of(const MatType& m) {
      return 
	+ get<0,0>::from(m)*get<1,1>::from(m)
	- get<1,0>::from(m)*get<0,1>::from(m);
    }
  };

  //======================================================================

  template<typename MatType,
	   typename Traits> 
  class DET<MatType, Traits, 3> {
  public:

    typedef MatType              M;
    typedef Traits               MT;
    typedef typename MT::ElType   T;

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
    };

    static T of(const MatType& m) {
      return 
	+ get<0,0>::from(m)*get<1,1>::from(m)*get<2,2>::from(m) 
	+ get<0,1>::from(m)*get<1,2>::from(m)*get<2,0>::from(m) 
	+ get<0,2>::from(m)*get<1,0>::from(m)*get<2,1>::from(m)
	- get<2,0>::from(m)*get<1,1>::from(m)*get<0,2>::from(m) 
	- get<2,1>::from(m)*get<1,2>::from(m)*get<0,0>::from(m) 
	- get<2,2>::from(m)*get<1,0>::from(m)*get<0,1>::from(m);
    }
  };

  //======================================================================
  // Det Template Functions
  //======================================================================

  /**
   * \brief Calculates the determinant of an object of type MatType.
   */
  template<typename T, 
	   template<typename, size_t, size_t,typename> class MatType, 
	   size_t NROW, size_t NCOL,
	   typename Traits>
  inline  T Det(const MatType<T,NROW,NCOL,Traits>& m) {
    return DET<MatType<T,NROW, NCOL,Traits>, Traits, NROW >::of(m); 
  }
  
  /**  \brief Overload to calculate the 2x2 determinant of a double array. */
  inline  double Det(ColMajorTraits<double,2>::ConstRefType m) { 
    return DET<ColMajorTraits<double,2>::Type, ColMajorTraits<double,2>,2>::of(m); }

  /**  \brief Overload to calculate the 3x3 determinant of a double array. */
  inline  double Det(ColMajorTraits<double,3>::ConstRefType m) { 
    return DET<ColMajorTraits<double,3>::Type,ColMajorTraits<double,3>,3>::of(m); 
  }

  /**  \brief Overload to calculate the 2x2 determinant of an int array. */
  inline  int Det(ColMajorTraits<int,2>::ConstRefType m) { 
    return DET<ColMajorTraits<int,2>::Type,ColMajorTraits<int,2>,2>::of(m); }

  /**  \brief Overload to calculate the 3x3 determinant of an int array. */
  inline  int Det(ColMajorTraits<int,3>::ConstRefType m) { 
    return DET<ColMajorTraits<int,3>::Type,ColMajorTraits<int,3>,3>::of(m); 
  }

} /* namespace psimag */

#endif /* PSIMAG_MatDet_H */
