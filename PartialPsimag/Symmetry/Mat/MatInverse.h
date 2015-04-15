//-*-C++-*-

#ifndef PSIMAG_MatInverse_H
#define PSIMAG_MatInverse_H

/** \file MatInverse.h
 *
 *  \brief Contains template functions providing fro matrix inverse
 *         operations.
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
   * \brief Template class to calculate the inverse of a 2x2
   *        matrix of type MatType storing the result in a square
   *        matrix of type MatTypeR.
   */
  template<typename T, 
	   template<typename, size_t, size_t, typename> class MatType, 
	   template<typename, size_t, size_t, typename> class MatTypeR,
	   typename Traits, typename RTraits>
  
  class AssembleInverse2
  {
  public:

    typedef MatType<T,2,2,Traits>   M;
    typedef MatTypeR<T,2,2,RTraits> MR;

    typedef Traits   MT;
    typedef RTraits MRT;

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
      static T& from(MR& im) {
	return MRT::template REF<MR,R,C>::GET(im);
      }
    };

    static void EXEC(const M& m, MR& im)
    {
      
      // Find reciprocal of determinant
      T iDet = Det(m);
      if(iDet==0)
	{
	  std::string message("Mat.h::Inverse(Mat<T,2,2>): ");
	  message += " attempt to invert matrix with zero determinant!";
	  throw std::range_error(message);
	}
      iDet = static_cast<T>(1)/iDet;
      get<0,0>::from(im) =  get<1,1>::from(m) * iDet;
      get<1,0>::from(im) = -get<1,0>::from(m) * iDet;
      get<0,1>::from(im) = -get<0,1>::from(m) * iDet;
      get<1,1>::from(im) =  get<0,0>::from(m) * iDet;
      
      return;
    }
  };

  /**
   * \brief Template class to calculate the inverse of a 2x2
   *        matrix of type MatType storing the result in a square
   *        matrix of type MatTypeR.
   */
  template<typename T, 
	   template<typename, size_t, size_t, typename> class MatType, 
	   template<typename, size_t, size_t, typename> class MatTypeR,
	   typename Traits, typename RTraits>
  class AssembleInverse3
  {
  public:
    
    typedef MatType<T,3,3,Traits>   M;
    typedef MatTypeR<T,3,3,RTraits> MR;

    typedef Traits   MT;
    typedef RTraits MRT;

    template<size_t R, size_t C> class get {
    public:
      static const T& from(const M& m) {
	return MT::template REF<M,R,C>::GETCONST(m);
      }
      static T& from(MR& im) {
	return MRT::template REF<MR,R,C>::GET(im);
      }
    };

    static void EXEC(const MatType<T,3,3,MT>& m, MatTypeR<T,3,3,MRT>& im)
    {
      
      // Find reciprocal of determinant
      T iDet = Det(m);
      if(iDet==0)
	{
	  std::string message("Mat.h::Inverse(Mat<T,3,3>): ");
	  message += " attempt to invert matrix with zero determinant!";
	  throw std::range_error(message);
	}
      iDet = static_cast<T>(1)/iDet;

      get<0,0>::from(im) =  ( get<1,1>::from(m) * get<2,2>::from(m) - get<1,2>::from(m) * get<2,1>::from(m) ) * iDet;
      get<0,1>::from(im) =  ( get<0,2>::from(m) * get<2,1>::from(m) - get<0,1>::from(m) * get<2,2>::from(m) ) * iDet;
      get<0,2>::from(im) =  ( get<0,1>::from(m) * get<1,2>::from(m) - get<0,2>::from(m) * get<1,1>::from(m) ) * iDet;
      
      get<1,0>::from(im) =  ( get<1,2>::from(m) * get<2,0>::from(m) - get<1,0>::from(m) * get<2,2>::from(m) ) * iDet;
      get<1,1>::from(im) =  ( get<0,0>::from(m) * get<2,2>::from(m) - get<0,2>::from(m) * get<2,0>::from(m) ) * iDet;
      get<1,2>::from(im) =  ( get<0,2>::from(m) * get<1,0>::from(m) - get<0,0>::from(m) * get<1,2>::from(m) ) * iDet;
      
      get<2,0>::from(im) =  ( get<1,0>::from(m) * get<2,1>::from(m) - get<1,1>::from(m) * get<2,0>::from(m) ) * iDet;
      get<2,1>::from(im) =  ( get<0,1>::from(m) * get<2,0>::from(m) - get<0,0>::from(m) * get<2,1>::from(m) ) * iDet;
      get<2,2>::from(im) =  ( get<0,0>::from(m) * get<1,1>::from(m) - get<0,1>::from(m) * get<1,0>::from(m) ) * iDet;

    }
  };

  //======================================================================
  /**
   * \brief Calculates the inverse of an object of type MatType and
   *        stores the result in MatTypeR.
   */
  template<typename T, 
	   template<typename, size_t, size_t, typename> class MatType, 
	   template<typename, size_t, size_t, typename> class MatTypeR,
	   typename Traits, typename RTraits,
	   size_t NROW, 
	   size_t NCOL>
  inline void Inverse(const MatType<T,NROW,NCOL,Traits>& m, 
		      MatTypeR<T,NROW,NCOL,RTraits>& im)
  {
    std::ostringstream msg;
    msg << "Did you want to use Inverse<t>(const Mat<T,3,3>& m, Mat<T,3,3>& im) instead??\n";
    msg << " Inverse(const Mat<T,NROW,NCOL>& m, Mat<T,NCOL,NROW>& im):\n";
    msg << " Not implemented Yet!\n";
    throw std::exception(msg.str());
    // Insert appropriate BLAS routine. We'll do this later.
  }

  /**
   * \brief Overload to calculate the inverse of an
   *        object of a square MatType of dimension 3.
   */
  template<typename T, 
	   template<typename, size_t, size_t, typename> class MatType, 
	   template<typename, size_t, size_t, typename> class MatTypeR,
	   typename Traits, typename RTraits>
  inline void Inverse(const MatType<T,2,2,Traits>& m, 
		      MatTypeR<T,2,2,RTraits>& im) {
    AssembleInverse2<T,MatType, MatTypeR,Traits,RTraits>::EXEC(m,im);
  }

  /**
   * \brief Overload to calculate the inverse of an
   *        object of a square MatType of dimension 3.
   */
  template<typename T, 
	   template<typename, size_t, size_t, typename> class MatType, 
	   template<typename, size_t, size_t, typename> class MatTypeR,
	   typename Traits, typename RTraits>
  inline void Inverse(const MatType<T,3,3,Traits>& m, 
		      MatTypeR<T,3,3,RTraits>& im) {
    AssembleInverse3<T,MatType, MatTypeR,Traits,RTraits>::EXEC(m,im);
  }

} /* namespace psimag */

#endif /* PSIMAG_MatInverse_H */

      //im[INDEX(0,0)] =  m[INDEX(1,1)] * iDet;
      //im[INDEX(1,0)] = -m[INDEX(1,0)] * iDet;
      //im[INDEX(0,1)] = -m[INDEX(0,1)] * iDet;
      //im[INDEX(1,1)] =  m[INDEX(0,0)] * iDet;
      //
//     im[INDEX(0,0)] = ( m[INDEX(1,1)]*m[INDEX(2,2)] - m[INDEX(1,2)]*m[INDEX(2,1)] ) * iDet;
//     im[INDEX(0,1)] = ( m[INDEX(0,2)]*m[INDEX(2,1)] - m[INDEX(0,1)]*m[INDEX(2,2)] ) * iDet;
//     im[INDEX(0,2)] = ( m[INDEX(0,1)]*m[INDEX(1,2)] - m[INDEX(0,2)]*m[INDEX(1,1)] ) * iDet;
//     //
//     im[INDEX(1,0)] = ( m[INDEX(1,2)]*m[INDEX(2,0)] - m[INDEX(1,0)]*m[INDEX(2,2)] ) * iDet;
//     im[INDEX(1,1)] = ( m[INDEX(0,0)]*m[INDEX(2,2)] - m[INDEX(0,2)]*m[INDEX(2,0)] ) * iDet;
//     im[INDEX(1,2)] = ( m[INDEX(0,2)]*m[INDEX(1,0)] - m[INDEX(0,0)]*m[INDEX(1,2)] ) * iDet;
//     //
//     im[INDEX(2,0)] = ( m[INDEX(1,0)]*m[INDEX(2,1)] - m[INDEX(1,1)]*m[INDEX(2,0)] ) * iDet;
//     im[INDEX(2,1)] = ( m[INDEX(0,1)]*m[INDEX(2,0)] - m[INDEX(0,0)]*m[INDEX(2,1)] ) * iDet;
//     im[INDEX(2,2)] = ( m[INDEX(0,0)]*m[INDEX(1,1)] - m[INDEX(0,1)]*m[INDEX(1,0)] ) * iDet;
    //
