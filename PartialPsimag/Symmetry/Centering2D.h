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

  template<size_t DIM>
  class Centerings {
  public:
    enum{ NumCenterings=0 };
  };

  template<>
  class Centerings<2> {
  public:
    enum{ NumCenterings=2 };
  };

  template<>
  class Centerings<3> {
  public:
    enum{ NumCenterings=6 }; // Probably more check this.
  };

  //======================================================================
  /**
   * \brief The Generic Centering2D template. Not used directly.
   */
  template<typename Field, size_t DIM, size_t NUM>
  class Centering {
  public:
  };

  //======================================================================
  /**
   * \brief The 2 Dimensions Centering2D template. 
   */
  template<typename Field>
  class Centering<2,0> {
  public:

    typedef std::vector<CellDirection<Field,DIM> > DirectionVectorType;

    static const std::string Name() { return "Primitive"; }

    static const DirectionVectorType Translations() {
      DirectionVectorType(0) result;
      return result;
    }

    /**
     * \brief Return the transformation from Primitive to Centered
     *        lattices. That is the identity in this case.
     *
     * \note We assume here that the input lattice is reduced and
     *       normalized.
     */
    static const TransformationType PrimitiveToCentered() {
      TransformationType result();
      return result;
    }
    /**
     * \brief Return the transformation from Centered To Primitive
     *        lattices. In this case the identity.
     *
     * \note We assume here that the input lattice is reduced and
     *       normalized.
     */
    static const TransformationType CenteredToPrimitive() {
      TransformationType result();
      return result;
    }
  };

  //======================================================================
  /**
   * \brief The 2 Dimensions Centering2D template. 
   */
  template<typename Field>
  class Centering<2,1> {
  public:

    typedef std::vector<CellDirection<Field,DIM> > DirectionVectorType;
    typedef LatticeTransformation<Field,2>         TransformationType;

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
