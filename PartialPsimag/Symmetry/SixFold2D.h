//-*-C++-*-

#ifndef PSIMAG_SixFoldSymmetryElement2D_H
#define PSIMAG_SixFoldSymmetryElement2D_H

#include "SymmetryElement.h"

/** \ingroup symmetryConcepts */

/*@{*/

/** \file  SixFold2D.h
 *
 *  Contains the two-dimensional SixFold Symmetry Element Class.
 *  
 */


namespace psimag {

  template<typename Field, typename Algorithms>
  class SixFold<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {

  public:

    enum { DIM=2 };

    typedef SixFold<Field,DIM,Algorithms>             SixFoldType;
    typedef SymmetryOperation<Field,DIM,Algorithms>   SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>     SymmetryElementType;
    typedef CellTranslation<Field,DIM>                CellTranslationType;
    typedef LatticeCoordinates<DIM>                   LatticeCoordinatesType;
 
    //======================================================================
    /**
     * \brief The default constructor produces a six-fold at the origin.
     */
    SixFold():
      SymmetryElementType("sixFold",
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellPosition<Field,int>(0,0),       // no offset
			  convert<Field>("0") )      // no glide
    {
      this->trace         = Field(1);
      this->sense         = 1;
      this->rotationAngle = Field(60);
    }
    
    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    SixFold(const CellTranslationType&    offset):
      SymmetryElementType("sixFold", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  offset, 
			  convert<Field>("0"))
    {
      this->trace         = Field(1);
      this->sense         = 1;
      this->rotationAngle = Field(60);
      //std::ostringstream buff;
      //buff << "SixFold(" << offset[0] << "," << offset[1] << ")";
      //this->name = buff.str();
    }

    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    template<typename IN_TYPE>
    SixFold(IN_TYPE p1, IN_TYPE p2):
      SymmetryElementType("sixFold", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellTranslation(convert<Field>(p1), convert<Field>(p2)), 
			  convert<Field>("0"))
    {
      this->trace         = Field(1);
      this->sense         = 1;
      this->rotationAngle = Field(60);
    }

    //======================================================================
    /** 
     * \brief Return the symmetry operation corresponding to this
     *        element in a given lattice.
     */
    SymmetryOperationType operator () (const Lattice<Field,DIM,Algorithms>& lattice) const {
      return operationFor(*this,lattice);
    }


   //======================================================================
    /** 
     * \brief Return the symmetry operation corresponding to this
     *        element in a given lattice.
     */
    static SymmetryOperationType operationFor(const SymmetryElementType&         element, 
					      const Lattice<Field,2,Algorithms>& lattice) {

      static  SymmetryOperationType result = 
	symmetryOperation<Field,std::string,Algorithms>("1/2",         "-1/2*sqrt(3)",
							"1/2*sqrt(3)", "1/2");
    
      result.setFixedPoint(lattice.cartesianPosition(element.cellPosition));
      result.name    = element.name;
      
      return result;
    }

  };


}

#endif 
/*@}*/


