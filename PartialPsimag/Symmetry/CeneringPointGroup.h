//-*-C++-*-

#ifndef PSIMAG_Centering2D_H
#define PSIMAG_Centering2D_H


/** \file Centering2D.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         ForEach operation over the Centerings. This operation
 *         applies a given function to centring for the specified
 *         dimension.
 *
 */
 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"
#include "MatTraits.h"

namespace psimag {

  template<typename Field, size_t DIM, size_t CenteringNum, size_t PointGroupNum>
  class CenteringPointGroup {
  public:
  };

  //======================================================================
  /**
   * \brief The 2 Dimensions Centering2D template. 
   */
  template<typename Field>
  class CenteringPointGroup<2,0,0> {
  public:
    static const std::string Name() { return "mmm"; }

    static const ???? searchMap
  };

  //======================================================================
  /**
   * \brief The 2 Dimensions Centering2D template. 
   */
  template<typename Field>
  class CenteringPointGroup<2,1,1> {
  public:

    typedef std::vector<CellDirection<Field,DIM> > DirectionVectorType;
    typedef LatticeTransformation<Field,2>         TransformationType;

    enum { NumPointGroups:2 }

    static const std::string Name() { return "BodyCentered"; }

    static const DirectionVectorType Translations() {
      DirectionVectorType(1) result;
      result[0] = cellDirection("1/2","1/2");
      return result;
    }

    /**
     * \brief Return the transformation from Primitive to BodyCentered
     *        lattices.
     *
     * \note We assume here that the input lattice is reduced and
     *       normalized.
     */
    static const TransformationType PrimitiveToCentered() {
      Vec<int,DIM> direction(2, -1);
      TransformationType result(direction);
      return result;
    }
    /**
     * \brief Return the transformation from BodyCentered To Primitive
     *        lattices.
     *
     * \note We assume here that the input lattice is reduced and
     *       normalized.
     */
    static const TransformationType CenteredToPrimitive() {

      // The inverse of:
      //Vec<int,DIM> direction(2, -1);
      //TransformationType result(direction);
      TransformationType result();
      return result;
    }
  };


} /* namespace psimag */

#endif /* PSIMAG_Centering2D */
