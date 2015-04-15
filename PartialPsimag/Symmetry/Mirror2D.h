//-*-C++-*-

#ifndef PSIMAG_Mirror2DSymmetry2D_H
#define PSIMAG_Mirror2DSymmetry2D_H

/** \ingroup symmetryConcepts */
/*@{*/

/** \file  Mirror2D.h
 *
 *  Contains the 2D Mirror sublass of symmetry element.
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
#include "CellTranslation.h"
#include "CellPosition.h"
#include "LatticeCoordinates.h"
#include "CartesianTranslation.h"
#include "FieldConvert.h"


#include "PSIMAGAssert.h"

namespace psimag {


  /** \ingroup symmetryConcepts
   *
   * \brief A class for implementing mirror symmetry elements.
   *
   * Symmetry elements store the geometric characteristics of a
   * SymmetryOperation and can produce the corresponding
   * SymmetryOperation or AppliedSymmetryOperation given a lattice or a
   * LatticeWithPattern.
   *
   */

  template<typename Field, typename Algorithms>
  class Mirror<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {
  public:

    enum { DIM=2 };

    typedef Mirror<Field,DIM,Algorithms>                       MirrorType;
    typedef SymmetryElement<Field,DIM,Algorithms>              SymmetryElementType;

    typedef SymmetryOperation<Field,DIM,Algorithms>            SymmetryOperationType;

    typedef CellTranslation<Field,DIM>                         CellTranslationType;
    typedef CellPosition<Field,DIM,Algorithms>                 CellPositionType;
    typedef LatticeCoordinates<DIM>                            LatticeCoordinatesType;
    typedef CartesianTranslation<Field,DIM>                    CartesianTranslationType;
    typedef CartesianPosition<Field,DIM>                       CartesianPositionType;

    typedef typename SymmetryOperationType::TranslationType    TranslationType;

    //======================================================================
    /**
     * \brief The default constructor produces a b-mirror symmetry element.
     */
    Mirror():
      SymmetryElementType("mirror",
			  latticeCoordinate(0,1),          // The b-axis direction (already normalized)
			  cellTranslation<Field,int>(0,0), // no offset
			  convert<Field>("0") )            // no glide
    {}
    
    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    Mirror(const LatticeCoordinatesType& orientation,
	   const CellTranslationType&    offset,
	   std::string                   tp="Mirror(@?)"):
      SymmetryElementType("mirror", 
			  orientation, 
			  offset, 
			  convert<Field>("0"))
    {}

    //======================================================================
    /** 
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and the zero translation.
     */
    Mirror(const LatticeCoordinatesType& orientation,
	   std::string                   tp="?"):
      SymmetryElementType("mirror", 
			  orientation, 
			  cellPosition<Field,Field>(0,0), 
			  convert<Field>("0"))
    {}

    //======================================================================
    /** 
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and an offset vector given
     *        by a lattice coordinate and a offset fraction.
     */
    template<typename IN_TYPE>
    Mirror(const LatticeCoordinatesType&  orientation,
	   const LatticeCoordinatesType&  offsetDir,
	   const IN_TYPE&                 offsetFraction,
	   std::string                    tp="?"):
      SymmetryElementType("mirror",
			  orientation,
			  CellPositionType(offsetDir.cellPosition<Field,Algorithms>() * convert<Field>(offsetFraction)),
			  convert<Field>("0"))
    {}

    //======================================================================
    //---------------------------------------------------------------------- Static A_Mirror

    static MirrorType* a() { 
      return new MirrorType(latticeCoordinate(1,0), "a-Mirror");
    }

    template<typename IN_TYPE>
    static MirrorType* a(IN_TYPE shift) { 
      MirrorType* result = new MirrorType(latticeCoordinate(1,0),
					  latticeCoordinate(0,1), 
					  shift);
      return result;
    }

    //---------------------------------------------------------------------- Static b_Mirror

    static MirrorType* b() { 
      return new MirrorType(latticeCoordinate(0,1), "b-Mirror");
    }

    template<typename IN_TYPE>
    static MirrorType* b(IN_TYPE shift) { 
      MirrorType* result = new MirrorType(latticeCoordinate(0,1),
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

      CartesianTranslationType mirrorLine = lattice.cartesianTranslation(element.netDirection);
      CartesianPositionType    cartPos    = lattice.cartesianPosition(element.cellPosition);

      SymmetryOperationType    result     = operation(mirrorLine, cartPos);

      result.name = element.name;

      return  result;
    }

    //---------------------------------------------------------------------- Static

    /**
     * \brief Returns the cartrsian SymmetryOperation for a mirror about the
     *        given cartesian axis.
     */
    static SymmetryOperationType operation(const CartesianTranslationType& mirrorLine) {
      
      Field l       = mirrorLine.length();
      Field cos_a   = mirrorLine[0]/l;
      Field sin_a   = mirrorLine[1]/l;
      Field squares = Algorithms::normalize((cos_a*cos_a)-(sin_a*sin_a));
      Field both    = Algorithms::normalize(Field(2) * sin_a * cos_a);
      
      
      SymmetryOperationType  op = symmetryOperation<Field,Field,Algorithms>(squares, both,
									    both,    Field(-1) * squares);
      return op;
    }	

    /**
     * \brief Returns the SymmetryOperation for the mirror about the
     *        axis in the direction of the given cartesian vector through the
     *        origin shifted by the given cartesian origin shift vector.
     */
    static SymmetryOperationType operation(const CartesianTranslationType& mirrorLine, 
					   const CartesianPositionType&    originPos) {

      SymmetryOperationType result = operation(mirrorLine);

      Multiply(result.rotation, originPos, result.translation);

      result.translation *= Field(-1);

      COPY<TranslationType, CartesianTranslationType,
	typename TranslationType::Traits,
	typename TranslationType::Traits>::EXEC_PLUS(result.translation, originPos);

      return result;
    }						   

  };

}

#endif 
/*@}*/
