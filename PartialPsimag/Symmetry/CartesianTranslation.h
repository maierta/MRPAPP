//-*-C++-*-
#ifndef  Psimag_Cartesian_Translation
#define  Psimag_Cartesian_Translation

/** \ingroup extendedMatrices */
/*@{*/

/** \file CartesianTranslation.h
 *
 * Contains classes for implementing CartesianTranslation (a subclass of SeitzVector) objects.
 */

#include "SeitzVector.h"
#include "FieldConvert.h"

namespace psimag {
  
  //============================================================
  /** \ingroup extendedMatrices
   *
   * \brief A class indicating that the translation is in the
   *        (implicit) cartesian system.
   */
  template<typename Field, size_t DIM> 
  class CartesianTranslation: public SeitzVector<Field, DIM, 0> {
    
  public:

    typedef          SeitzVector<Field,DIM,0> BaseType;
    typedef typename BaseType::ArrayType      ArrayType;
    typedef CartesianPosition<Field,DIM>      CartesianPositionType;

    // Constructors are never inherited :-(

    /** The Default Constructor produces the zero translation. */
    CartesianTranslation(): SeitzVector<Field, DIM, 0>() {}

    /** Construct a translation whose components are set to the given value. */
    template<typename IN_TYPE> CartesianTranslation(const IN_TYPE& val): SeitzVector<Field, DIM, 0>(val) {}

     /** Construct a translation whose components are set to the given value. */
    CartesianTranslation(const Vec<Field, DIM>& v):  SeitzVector<Field, DIM, 0>(v) {}

    /** Construct a translation whose components are set to the given value. */
    CartesianTranslation(const CartesianTranslation<Field, DIM>& v): SeitzVector<Field, DIM, 0>(v) {}

    /** Construct a translation whose components are set to the given value array. */
    CartesianTranslation(const ArrayType* vals):  SeitzVector<Field, DIM, 0>(vals) {}

    /** 
     * Return the length of the cartesian translation. 
     *
     * Convert data types as necessary.
     */
    Field length() const {
      return sqrt((*this) * (*this));
    }

    CartesianPositionType plusOrigin() const {
      static CartesianPositionType origin = CartesianPositionType();
      return CartesianPositionType((*this) + origin);
    }

//     /** 
//      *
//      * Difference operator, translation concatenation
//      */
//     template<typename Field, size_t DIM>
//     CartesianTranslation<Field,DIM> operator- (const CartesianTranslation<Field,DIM>& lhs, 
// 					       const CartesianTranslation<Field,DIM>& rhs)
//     {
//       return static_cast<SeitzVector<Field, DIM> >(lhs) 
// 	- static_cast<SeitzVector<Field, DIM> >(rhs);
//     }
  
//     /** 
//      *
//      * Plus operator, translation concatenation
//      */
//     template<typename Field, size_t DIM>
//     CartesianTranslation<Field,DIM> operator+ (const CartesianTranslation<Field,DIM>& lhs, 
// 					       const CartesianTranslation<Field,DIM>& rhs)
//     {
//       return static_cast<SeitzVector<Field, DIM> >(lhs) 
// 	- static_cast<SeitzVector<Field, DIM> >(rhs);
//     }
    
  };

  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<typename Field, size_t DIM> 
  Tag toXML(const CartesianTranslation<Field,DIM> pos,
	    std::string name="CartesianTranslation") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << pos[i] ;
    return result;
  }

  template<typename Field, typename IN_TYPE>
  CartesianTranslation<Field,2> cartesianTranslation(IN_TYPE c0, IN_TYPE c1) {
    CartesianTranslation<Field,2> result;
    result[0] = convert<Field>(c0);
    result[1] = convert<Field>(c1);
    return result;
  }

  template<typename Field, typename IN_TYPE>
  CartesianTranslation<Field,3> cartesianTranslation(IN_TYPE c0, IN_TYPE c1, IN_TYPE c2) {
    CartesianTranslation<Field,3> result;
    result[0] = convert<Field>(c0);
    result[1] = convert<Field>(c1);
    result[2] = convert<Field>(c2);
    return result;
  }

}  /* namspace psimag */
#endif // Psimag_Cartesian_Translation
/*@}*/
