//-*-C++-*-

#ifndef PSIMAG_GlideSymmetry2D_H
#define PSIMAG_GlideSymmetry2D_H

/** \ingroup symmetryConcepts */
/*@{*/

/** \file  Glide.h
 *
 *  Contains the 2D Glide sublass of symmetry element.
 *  
 */
 
#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Mat.h"
#include "Real.h"

#include "SymmetryOperation.h"
#include "CellDirection.h"
#include "CellPosition.h"
#include "CartesianTranslation.h"
#include "FieldConvert.h"

#include "PSIMAGAssert.h"

namespace psimag {


  //============================================================ Fundamental Glide Functions
  
  /** \ingroup symmetryConcepts
   *
   * \brief A class for implementing glide symmetry elements.
   *
   * Symmetry elements store the geometric characteristics of a
   * SymmetryOperation and can produce the corresponding
   * SymmetryOperation or AppliedSymmetryOperation given a lattice or a
   * LatticeWithPattern.
   *
   */

  template<typename Field, typename Algorithms>
  class Glide<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {
  public:

    enum { DIM=2 };

    typedef Glide<Field,DIM,Algorithms>                           GlideType;
    typedef Mirror<Field,DIM,Algorithms>                          MirrorType;
    typedef SymmetryOperation<Field,DIM,Algorithms>               SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>                 SymmetryElementType;
    typedef CellTranslation<Field,DIM>                            CellTranslationType;
    typedef CellPosition<Field,DIM,Algorithms>                    CellPositionType;
    typedef CartesianTranslation<Field,DIM>                       CartesianTranslationType;
    typedef CartesianPosition<Field,DIM>                          CartesianPositionType;
    typedef LatticeCoordinates<DIM>                               LatticeCoordinatesType;
    typedef typename SymmetryOperationType::TranslationType       TranslationType;
    typedef Lattice<Field,DIM,Algorithms>                         LatticeType;
    typedef std::vector<SymmetryOperationType>                    RatedOperationsType;

    //======================================================================
    /**
     * \brief The default constructor produces a b-glide symmetry element.
     */
    Glide():
      SymmetryElementType("glide",
			  latticeCoordinate(0,1),    // The b-axis direction
			  cellPosition<Field,Field>(0,0),       // no offset
			  convert<Field>("1/2") )    // the ususal glide
    {}
    
    //======================================================================
    /**
     * \brief Construct a glide symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    template<typename IN_TYPE>
    Glide(const LatticeCoordinatesType& orientation,
	  const IN_TYPE&                glideFraction,
	  const CellTranslationType&    offset):
      SymmetryElementType("glide", 
			  orientation, 
			  offset, 
			  convert<Field>(glideFraction))
    {}

    //======================================================================
    /** 
     * \brief Construct a glide symmetry element from a direction
     *        given by a lattice coordinate and the zero translation.
     */
    template<typename IN_TYPE>
    Glide(const LatticeCoordinatesType& orientation,
	  const IN_TYPE&                glideFraction):
      SymmetryElementType("glide", 
			  orientation, 
			  cellPosition<Field,Field>(0,0), 
			  convert<Field>(glideFraction))
    {}

    //======================================================================
    /** 
     * \brief Construct a glide symmetry element from a direction
     *        given by a lattice coordinate and an offset vector given
     *        by a lattice coordinate and a offset fraction.
     */
    template<typename IN_TYPE1, typename IN_TYPE2>
    Glide(const LatticeCoordinatesType&  orientation,
	  const IN_TYPE1&                glideFraction,
	  const LatticeCoordinates<DIM>& offset,
	  const IN_TYPE2&                offsetFraction):
      SymmetryElementType("glide",
			  orientation,
			  CellPositionType(offset.cellPosition<Field,Algorithms>() * convert<Field>(offsetFraction)),
			  convert<Field>(glideFraction))
    {}

    //====================================================================== Static a_Glide

    static GlideType* a() { 
      return new GlideType(latticeCoordinate(1,0), "1/2");
    }

    template<typename IN_TYPE>
    static GlideType* a(IN_TYPE shift) { 
      GlideType* result = new GlideType(latticeCoordinate(1,0), "1/2",
					latticeCoordinate(0,1), 
					shift);
      return result;
    }
    
    //---------------------------------------------------------------------- Static b_Glide

    static GlideType* b() { 
      return new GlideType(latticeCoordinate(0,1), "1/2");
    }

    template<typename IN_TYPE>
    static GlideType* b(IN_TYPE shift) { 
      GlideType* result = new GlideType(latticeCoordinate(0,1), "1/2",
					latticeCoordinate(1,0), 
					shift);
      return result;
    }

    //======================================================================
    /** 
     * \brief Return the symmetry operation corresponding to this
     *        element in a given lattice.
     */
    virtual SymmetryOperationType operator () (const Lattice<Field,DIM,Algorithms>& lattice) const {
      return operationFor(*this,lattice);
    }


   //======================================================================
    /** 
     * \brief Return the symmetry operation corresponding to this
     *        element in a given lattice.
     */
    static SymmetryOperationType operationFor(const SymmetryElementType&         element, 
					      const Lattice<Field,2,Algorithms>& lattice) {

      typedef ColMajorTraits<Field,DIM,1> VTraits;
      typedef typename SymmetryOperationType::TranslationType STranslationType;

      //      static bool trace(false);

      CartesianTranslationType glideLine     = lattice.cartesianTranslation(element.netDirection);
      CartesianPositionType    cartPos       = lattice.cartesianPosition   (element.cellPosition);
      CartesianTranslationType cartTransPart = lattice.cartesianTranslation(element.translationPart);

      SymmetryOperationType result = MirrorType::operation(glideLine,cartPos);

      COPY<STranslationType,CartesianTranslationType,VTraits,VTraits>::EXEC_PLUS(result.translation, cartTransPart);
      
      result.name = element.name;
      return result;
    }

  private:

    //====================================================================== Fundamental Glide Functions

    //====================================================================== 
    /**
     * \brief Given a cartesian glide vector and ratio returns the cartrsian
     *        translation part of a cartesian glide symmetry
     *        operation.
     *
     * \param glideVector: The direction of the glide operation.
     *
     * \paran ratio: The proportion of the glide vector that the
     *               pattern is moved after reflection.
     */
    template<typename IN_TYPE>
    static TranslationType getTranslationPart(CartesianTranslationType glideVector, 
					      IN_TYPE ratio) {

      Field l           = glideVector.length();
      Field cos_a       = glideVector[0]/l;
      Field sin_a       = glideVector[1]/l;
      Field glideLength = convert<Field>(ratio) * l; /** The glide length */
      
      return TranslationType( Algorithms::normalize(cos_a*glideLength), 
			      Algorithms::normalize(sin_a*glideLength) ) ;
    }

//     //====================================================================== 
//     /**
//      * \brief Returns the cartrsian symmetry operation for a glide
//      *        about an axis in the direction of the given vector
//      *        through the origin. The glide distance is specified as a
//      *        percentage of the glideLine vector.
//      */
//     template<typename IN_TYPE>
//     static SymmetryOperationType operation(CartesianTranslationType glideLine, 
// 					   IN_TYPE glideRatio) {
      
//       SymmetryOperationType result = MirrorType::operation(glideLine);
//       result.translation          += getTranslationPart(glideLine, glideRatio); 
//       return result;
//     }

//     //------------------------------------------------------------

//     /**
//      * \brief Returns the Cartrsian Symmetry operation for a glide
//      *        about the axis in the direction of the given vector and the
//      *        origin shifted by the given originShift. The glide
//      *        distance is specified as a percentage of the glide vector.
//      */
//     template<typename IN_TYPE>
//     static SymmetryOperationType operation(CartesianTranslationType glideLine, 
// 					   IN_TYPE                  glideRatio,
// 					   CartesianPositionType    originPos) {
      
//       SymmetryOperationType result = MirrorType::operation(glideLine,originPos);
//       result.translation += getTranslationPart(glideLine, glideRatio); 
//       return result;
//     }
  };


}

#endif 
/*@}*/
