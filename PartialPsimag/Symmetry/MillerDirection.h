//-*-C++-*-
#ifndef  Psimag_Miller_Direction
#define  Psimag_Miller_Direction

/** \ingroup extendedMatrices */
/*@{*/

/** \file MillerDirection.h
 *
 * Contains a class for implementing MillerDirection (SeitzTranslation subclass) objects.
 */

#include "SeitzVector.h"

namespace psimag {
  
  /** \ingroup extendedMatrices
   *
   * \brief A marker class indicating that the Translation indicates a direction from within a Cell.
   *
   * These can be converted to  and from MillerIndices for a given Cell.
   */
  template<typename Field, size_t DIM> 
  class MillerDirection: public SeitzVector< Field, DIM, 0 > {
    
  public:

    typedef          SeitzVector<Field,DIM,0> BaseType;
    typedef typename BaseType::ArrayType      ArrayType;

    // Constructors are never inherited :-(

    /** The Default Constructor produces the zero translation. */
    MillerDirection(): SeitzVector<Field, DIM, 0>() {}

    /** Construct a translation whose components are set to the given value. */
    template<typename IN_TYPE> MillerDirection(const IN_TYPE& val): SeitzVector<Field, DIM, 0>(val) {}

     /** Construct a translation whose components are set to the given value. */
    MillerDirection(const Vec<Field, DIM>& v):  SeitzVector<Field, DIM, 0>(v) {}

    /** Construct a translation whose components are set to the given value. */
    MillerDirection(const MillerDirection<Field, DIM>& v): SeitzVector<Field, DIM, 0>(v) {}

    /** Construct a translation whose components are set to the given value array. */
    MillerDirection(const ArrayType& vals):  SeitzVector<Field, DIM, 0>(vals) {}

  };

  template<typename Field, typename IN_TYPE> 
  MillerDirection<Field,2> millerDirection(IN_TYPE t0, IN_TYPE t1) {
    MillerDirection<Field,2> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    return result;
  }

  template<typename Field, typename IN_TYPE> 
  MillerDirection<Field,2> millerDirection(IN_TYPE t0, IN_TYPE t1, IN_TYPE t2) {
    MillerDirection<Field,2> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    result[2] = convert<Field>(t2);
    return result;
  }
  

}  /* namspace psimag */

#endif // Psimag_Miller_Direction
/*@}*/
