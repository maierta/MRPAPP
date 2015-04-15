//-*-C++-*-
#ifndef  Psimag_Cartesian_Position
#define  Psimag_Cartesian_Position

/** \ingroup extendedMatrices */
/*@{*/

/** \file CartesianPosition.h
 *
 * Contains classes for implementing CartesianPosition (a subclass of SeitzPosition) objects.
 */

#include "SeitzVector.h"
#include "SeitzPosition.h"
#include "Tag.h"

namespace psimag {

  //============================================================
  /** \ingroup extendedMatrices
   *
   * \brief A class indicating that the position is in the
   *        (implicit) cartesian system.
   */
  template<typename Field, size_t DIM> 
  class CartesianPosition: public SeitzPosition<Field,DIM> {
    
  public:

    typedef          SeitzPosition<Field,DIM> BaseType;
    typedef typename BaseType::ArrayType      ArrayType;

   /** The Default Constructor produces the zero position. */
    CartesianPosition(): BaseType() {}

    /** Construct a position whose components are set to the given value. */
    CartesianPosition(const Field& val): BaseType(val) {}

    /** Construct a position whose components are set to the given value. */
    //CartesianPosition(const Field& val): BaseType(val) {}

     /** Construct a position whose components are set to the given value. */
    CartesianPosition(const Vec<Field, DIM>& v):  BaseType(v) {}

    /** Construct a position whose components are set to the given value. */
    CartesianPosition(const SeitzVector<Field,DIM,1>& v): BaseType(v) {}

    /** Construct a position whose components are set to the given value. */
    CartesianPosition(const BaseType& v): BaseType(v) {}

    /** Construct a position whose components are set to the given value. */
    CartesianPosition(const CartesianPosition<Field, DIM>& v): BaseType(v) {}

    /** Construct a position whose components are set to the given value array. */
    CartesianPosition(const ArrayType& vals):  BaseType(vals) {}

  };
    
  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<typename Field, size_t DIM> 
  Tag toXML(const CartesianPosition<Field,DIM> pos,
	    std::string name = "CartesianPosition") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << pos[i] ;
    return result;
  }

  //====================================================================== 

  template<typename Field, typename IN_TYPE>
  CartesianPosition<Field,2> cartesianPosition(IN_TYPE c0, IN_TYPE c1) {
    CartesianPosition<Field,2> result;
    result[0] = convert<Field>(c0);
    result[1] = convert<Field>(c1);
    return result;
  }

  template<typename Field, typename IN_TYPE>
  CartesianPosition<Field,3> cartesianPosition(IN_TYPE c0, IN_TYPE c1, IN_TYPE c2) {
    CartesianPosition<Field,3> result;
    result[0] = convert<Field>(c0);
    result[1] = convert<Field>(c1);
    result[2] = convert<Field>(c2);
    return result;
  }

//   /** 
//    *
//    * Difference operator, position concatenation
//    */
//   template<typename Field, size_t DIM>
//   CartesianPosition<Field,DIM> operator- (const CartesianPosition<Field,DIM>& lhs, 
// 					  const CartesianPosition<Field,DIM>& rhs)
//   {
//     return static_cast<SeitzVector<Field, DIM> >(lhs) 
//       - static_cast<SeitzVector<Field, DIM> >(rhs);
//   }
  
}  /* namspace psimag */

#endif // Psimag_Cartesian_Position
/*@}*/
