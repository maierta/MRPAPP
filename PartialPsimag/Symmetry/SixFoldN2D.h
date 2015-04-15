//-*-C++-*-

#ifndef PSIMAG_SixFoldNSymmetryElement2D_H
#define PSIMAG_SixFoldNSymmetryElement2D_H

/** \ingroup symmetryConcepts */

/*@{*/

/** \file  SixFoldN2D.h
 *
 *  Contains the two-dimensional SixFoldN Symmetry Element Class.
 *  
 */


namespace psimag {


  /**\ingroup symmetryConcepts
   *
   * \brief A class for implementing negative-sense, two-dimensional SixFold Symmetry Elements.
   *
   */
  template<typename Field, typename Algorithms>
  class SixFoldN<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {

  public:

    enum { DIM=2 };

    typedef SixFoldN<Field,DIM,Algorithms>            SixFoldNType;
    typedef SymmetryOperation<Field,DIM,Algorithms>   SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>     SymmetryElementType;
    typedef CellTranslation<Field,DIM>                CellTranslationType;
    typedef LatticeCoordinates<DIM>                   LatticeCoordinatesType;
 
    //======================================================================
    /**
     * \brief The default constructor produces a six-fold at the origin.
     */
    SixFoldN():
      SymmetryElementType("sixFoldN",
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellPosition<Field,Field>(0,0),       // no offset
			  convert<Field>("0") )      // no glide
    {
      this->trace = Field(1);
      this->sense = -1;
      this->rotationAngle = Field(-60);
    }
    
    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    SixFoldN(const CellTranslationType&    offset):
      SymmetryElementType("sixFoldN", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  offset, 
			  convert<Field>("0"))
    {
      this->trace = Field(1);
      this->sense = -1;
      this->rotationAngle = Field(-60);
    }

    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    template<typename IN_TYPE>
    SixFoldN(IN_TYPE p1, IN_TYPE p2):
      SymmetryElementType("sixFoldN", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellTranslation(convert<Field>(p1), convert<Field>(p2)), 
			  convert<Field>("0"))
    {
      this->trace = Field(1);
      this->sense = -1;
      this->rotationAngle = Field(-60);
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
	symmetryOperation<Field,std::string,Algorithms>("1/2",         "1/2*sqrt(3)",
							"-1/2*sqrt(3)", "1/2");
    
      result.setFixedPoint(lattice.cartesianPosition(element.cellPosition));
      result.name    = element.name;
      
      return result;
    }

  };
}
#endif 
/*@}*/


