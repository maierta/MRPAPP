//-*-C++-*-

#ifndef PSIMAG_TwoFoldSymmetryElement2D_H
#define PSIMAG_TwoFoldSymmetryElement2D_H

#include "CellTranslation.h"

/** \ingroup symmetryConcepts */

/*@{*/

/** \file  TwoFold2D.h
 *
 *  Contains the two-dimensional TwoFold Symmetry Element Class.
 *  
 */

namespace psimag {


  /**\ingroup symmetryConcepts
   *
   * \brief A class for implementing two-dimensional TwoFold Symmetry Elements.
   *
   */

  template<typename Field, typename Algorithms>
  class TwoFold<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {

  public:

    enum { DIM=2 };

    typedef TwoFold<Field,DIM,Algorithms>             TwoFoldType;
    typedef SymmetryOperation<Field,DIM,Algorithms>   SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>     SymmetryElementType;
    typedef CellTranslation<Field,DIM>                CellTranslationType;
    typedef LatticeCoordinates<DIM>                   LatticeCoordinatesType;
 
    //======================================================================
    /**
     * \brief The default constructor produces a two-fold at the origin.
     */
    TwoFold():
      SymmetryElementType("twoFold",
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellPosition<Field,Field>(0,0),       // no offset
			  convert<Field>("0") )      // no glide
    {}
    
    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    TwoFold(const CellTranslationType&    offset):
      SymmetryElementType("twoFold", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  offset, 
			  convert<Field>("0"))
    {}

    //======================================================================
    /**
     * \brief Construct a mirror symmetry element from a direction
     *        given by a lattice coordinate and a offset given by a
     *        cell translation.
     */
    template<typename IN_TYPE>
    TwoFold(IN_TYPE p1, IN_TYPE p2):
      SymmetryElementType("twoFold", 
			  latticeCoordinate(0,0),    // => the c-axis direction
			  cellPosition<Field,Field>(convert<Field>(p1), convert<Field>(p2)), 
			  convert<Field>("0"))
    {}

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
	symmetryOperation<Field,int,Algorithms>(-1, 0, 0, -1);
      
      result.setFixedPoint(lattice.cartesianPosition(element.cellPosition));
      result.name    = element.name;
      
      return result;
    }
    
  };

}
#endif 
/*@}*/


