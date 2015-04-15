//-*-C++-*-
#ifndef  Psimag_Cell_Direction
#define  Psimag_Cell_Direction

/** \ingroup extendedMatrices */
/*@{*/

/** \file CellDirection.h
 *
 * Contains classes for implementing directions in cell coordinates.
 */

#include "Vec.h"
#include <iostream>
#include "FieldConvert.h"
#include "SeitzVectors.h"
#include "Mat.h"

namespace psimag {
  
  /** \ingroup extendedMatrices
   *
   * \brief A class representing a direction which is defined in cell
   * coordinates. Which cell is implicit.
   *
   * If the coordinates of a CellDirection are integers and the
   * implicit cell is primitive then the resulting directions are
   * lattice vectors.
   *
   * If the coordinates are [0,0,1] , [0,1,0] or [1,0,0] then the
   * CellDirection represents the same direction as the
   * corresponding basis vectors.
   *
   */
  template<typename Field, size_t DIM> 
  class CellDirection: public SeitzVector<Field,DIM,0> {
    
  public:

    typedef          SeitzVector<Field,DIM,0> BaseType;
    typedef typename BaseType::ArrayType      ArrayType;

    /** 
     * The Default Constructor produces the zero translation.
     */
    CellDirection(): SeitzVector<Field,DIM,0>() {}

    /** 
     * Construct a translation whose components are set to the given value. 
     *
     * Convert data types as necessary.
     */
    template<typename IN_TYPE>
    CellDirection(const IN_TYPE& val): SeitzVector<Field,DIM,0>(val) {}

    /** 
     * Construct a translation whose components are set to the given value array. 
     *
     * Convert data types as necessary.
     */
    CellDirection(const ArrayType vals): SeitzVector<Field,DIM,0>(vals) { }

  };

  template<typename Field, typename IN_TYPE>
  CellDirection<Field,2> cellDirection(IN_TYPE x, IN_TYPE y) {
    CellDirection<Field,2> result;
    result[0] = convert<Field>(x);
    result[1] = convert<Field>(y);
    return result;
  }

  template<typename Field, typename IN_TYPE>
  CellDirection<Field,3> cellDirection(IN_TYPE x, IN_TYPE y, IN_TYPE z) {
    CellDirection<Field,3> result;
    result[0] = convert<Field>(x);
    result[1] = convert<Field>(y);
    result[2] = convert<Field>(z);
    return result;
  }


}  /* namspace psimag */

#endif
/*@}*/
