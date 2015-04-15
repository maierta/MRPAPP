//-*-C++-*-

#ifndef PSIMAG_ThreeFoldNSymmetryElement2D_H
#define PSIMAG_ThreeFoldNSymmetryElement2D_H

#include "LatticeCoordinates.h"
#include "CellTranslation.h"

/** \ingroup symmetryConcepts */

/*@{*/

/** \file  ThreeFoldN2D.h
 *
 *  Contains the two-dimensional ThreeFold Symmetry Element Class.
 *  
 */

namespace psimag {


  /**\ingroup symmetryConcepts
   *
   * \brief A class for implementing negative sense two-dimensional ThreeFold Symmetry Elements.
   *
   */
  template<typename Field, typename Algorithms>
  class ThreeFoldN<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {

  public:

    enum { DIM=2 };

    typedef ThreeFold<Field,DIM,Algorithms>           ThreeFoldType;
    typedef SymmetryOperation<Field,DIM,Algorithms>   SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>     SymmetryElementType;
    typedef CellTranslation<Field,DIM>                CellTranslationType;
    typedef LatticeCoordinates<DIM>                   LatticeCoordinatesType;
 
    //======================================================================
    /**
     * \brief The default constructor produces a three-fold at the origin.
     */
    ThreeFoldN():
      SymmetryElementType("threeFoldN",
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellPosition<Field,Field>(0,0),       // no offset
			  convert<Field>("0") )      // no glide
    {
      this->trace         = Field(-1);
      this->sense         = -1;
      this->rotationAngle = Field(-120);
    }
    
    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    ThreeFoldN(const CellTranslationType&    offset):
      SymmetryElementType("threeFoldN", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  offset, 
			  convert<Field>("0"))
    {

      this->trace         = Field(-1);
      this->sense         = -1;
      this->rotationAngle = Field(-120);
    }

    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    template<typename IN_TYPE>
    ThreeFoldN(IN_TYPE p1, IN_TYPE p2):
      SymmetryElementType("threeFoldN", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellPosition<Field,Field>(convert<Field>(p1), convert<Field>(p2)), 
			  convert<Field>("0"))
    {

      this->trace         = Field(-1);
      this->sense         = -1;
      this->rotationAngle = Field(-120);
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
    static SymmetryOperationType operationFor (const SymmetryElementType&           element,
					       const Lattice<Field,DIM,Algorithms>& lattice)
    {
      static  SymmetryOperationType result = 
	symmetryOperation<Field,std::string,Algorithms>("-1/2",           "1/2*sqrt(3)",
							"-1/2*sqrt(3)",   "-1/2"           );
      result.setFixedPoint(lattice.cartesianPosition(element.cellPosition));
      result.name    = element.name;
      
      return result;
    }

  };
}

#endif 
/*@}*/


