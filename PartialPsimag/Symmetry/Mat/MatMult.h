//-*-C++-*-

#ifndef PSIMAG_MatMultvec_H
#define PSIMAG_MatMultvec_H

/** \file MatMultvec.h
 *
 *  \brief Contains template functions providing Multvec for
 *         standard operations on Matrices of any sort.
 *
 */

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"
#include "Mat.h"
#include "MatTraits.h"


namespace psimag {

  //======================================================================
  // &*&*&*&* 
  // Generalize so that it can handle matricies of different types
  //======================================================================

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    size_t   LROW,
	    size_t   LCOL,
	    size_t   RCOL>
  class MatMultVecPos {
  public:
    typedef typename Traits::ElType T;
    enum { RROW=LCOL };
      
    template<size_t R, size_t C> class get {
    public:
      static const T& from(const MatType& m) {
	return Traits::template REF<MatType,R,C>::GETCONST(m);
      }
    };

//     template<size_t R, size_t C> class getO {
//     public:
//       static const T& from(const MatTypeO& m) {                    // &*&*&*&* can we have it both ways does it make a difference?
// 	return TraitsO::template REF<MatTypeO,R,C>::GETCONST(m);   // If the field is different then there is an implicit temp
//       }
//     };

    template<size_t R, size_t C> class getO {
    public:
      static T from(const MatTypeO& m) {                    // &*&*&*&* can we have it both ways does it make a difference?
	return TraitsO::template REF<MatTypeO,R,C>::GETCONST(m);   // If the field is different then there is an implicit temp
      }
    };

    static T EXEC(const MatType& lhs, const MatTypeO& rhs) {

      return get<LROW,LCOL>::from(lhs)*getO<RROW,RCOL>::from(rhs)
	+ MatMultVecPos<MatType, Traits, MatTypeO, TraitsO, LROW, LCOL-1, RCOL>::EXEC(lhs,rhs);
    }
  };

  //======================================================================

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    size_t   LROW,
	    size_t   RCOL>
  class MatMultVecPos<MatType, Traits, MatTypeO, TraitsO, LROW, 0, RCOL > {
  public:
    typedef typename Traits::ElType T;
    enum { LCOL=0, RROW=LCOL };
    
    template<size_t R, size_t C> class get {
    public:
      static const T& from(const MatType& m) {
	return Traits::template REF<MatType,R,C>::GETCONST(m);
      }
    };

    template<size_t R, size_t C> class getO {
    public:
      static const T from(const MatTypeO& m) {
	return TraitsO::template REF<MatTypeO,R,C>::GETCONST(m);
      }
    };
    
    static T EXEC(const MatType& lhs, const MatTypeO& rhs) {

      return get<LROW,LCOL>::from(lhs)*getO<RROW,RCOL>::from(rhs);
    }
  };

  //====================================================================== ROW
  //======================================================================

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    typename MatTypeR,
	    typename TraitsR,
	    size_t   LROW,
	    size_t   RCOL>
  class MatMultVecRow {
  public:
    enum { NLCOL=Traits::NCOL};

    typedef typename Traits::ElType  T;
    typedef typename TraitsR::ElType TR;

    enum { LCOL=0, RROW=LCOL };
    
    template<size_t R, size_t C> class get {
    public:
      static TR& from(MatTypeR& m) {
	return TraitsR::template REF<MatTypeR,R,C>::GET(m);
      }
    };

    static void EXEC(const MatType& m, const MatTypeO& v, MatTypeR& result) {

      get<LROW,RCOL>::from(result) = static_cast<TR>(MatMultVecPos<MatType, Traits, MatTypeO, TraitsO, LROW, NLCOL-1, RCOL>::EXEC(m,v));
      MatMultVecRow<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, LROW-1, RCOL>::EXEC(m,v, result);
    }
  };

  //====================================================================== ROW 0

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    typename MatTypeR,
	    typename TraitsR,
	    size_t   RCOL>
  class MatMultVecRow<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, 0, RCOL> {
  public:
    typedef typename Traits::ElType  T;
    typedef typename TraitsR::ElType TR;

    enum { LROW=0, NLCOL=Traits::NCOL};
    
    template<size_t R, size_t C> class get {
    public:
      static TR& from(MatTypeR& m) {
	return TraitsR::template REF<MatTypeR,R,C>::GET(m);
      }
    };

    static void EXEC(const MatType& m, const MatTypeO& v, MatTypeR& result) {

      get<LROW,RCOL>::from(result) = static_cast<TR>(MatMultVecPos<MatType, Traits, MatTypeO, TraitsO, LROW, NLCOL-1, RCOL>::EXEC(m,v));
    }
  };

  //====================================================================== COL

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    typename MatTypeR,
	    typename TraitsR,
	    size_t   RCOL>
  class MatMultVecCol {
  public:
    enum { NLROW=Traits::NROW};
    
    static void EXEC(const MatType& m, const MatTypeO& v, MatTypeR& result) {

      MatMultVecRow<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, NLROW-1, RCOL>::EXEC(m,v,result);
      MatMultVecCol<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, RCOL-1>::EXEC(m,v,result);
    }
  };

  //====================================================================== COL 0


  /**
   *
   * The scope of the multiply operation is controlled by NRRow and NLCOL.
   */

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    typename MatTypeR,
	    typename TraitsR >
  class MatMultVecCol<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, 0> {
  public:
    enum { RCOL= 0, NLROW=Traits::NROW};

    static void EXEC(const MatType& m, const MatTypeO& v, MatTypeR& result) {

      MatMultVecRow<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, NLROW-1, RCOL>::EXEC(m,v,result);
    }
  };

  //====================================================================== Whole Thing

  template< typename MatType, 
	    typename Traits,
	    typename MatTypeO,
	    typename TraitsO,
	    typename MatTypeR,
	    typename TraitsR >
  class MatMultVec {
  public:
    
    enum { NRCOL =TraitsO::NCOL};
    
    static void EXEC(const MatType& m, const MatTypeO& mo, MatTypeR& result) {
      MatMultVecCol<MatType, Traits, MatTypeO, TraitsO, MatTypeR, TraitsR, NRCOL-1>::EXEC(m,mo,result);
    }
  };

  //====================================================================== Functions and Operators
  
 /**
   * \brief Generic Multiply M*V = V.
   */
      template< typename MatType, typename MatTypeO, typename MatTypeR >

      inline void Multiply(const MatType& m,
			   const MatTypeO& v,
			         MatTypeR& result) { 
	
	MatMultVec<MatType, typename MatType::Traits,
	  MatTypeO, typename MatTypeO::Traits,
	  MatTypeR, typename MatTypeR::Traits >::EXEC(m,v,result);
      }

  /**
   * \Brief Generic Multiply M*M3 = M3 operator.
   */
  template<typename MatType,
	   typename MatTypeO,
	   typename MatTypeR>
  inline MatTypeR
  operator * (const MatType&  m, 
	      const MatTypeO& other) {
    MatTypeR result;
    Multiply(m, other, result);
    return result;
  }

} /* namespace psimag */

#endif /* PSIMAG_MatMultvec_H */
