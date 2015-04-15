//-*-C++-*-

#ifndef PSIMAG_IdentitySymmetryElement2D_H
#define PSIMAG_IdentitySymmetryElement2D_H

/** \ingroup symmetryConcepts */

/*@{*/

/** \file  Identity2D.h
 *
 *  Contains the two-dimensional Identity Symmetry Element Class.
 *  
 */

namespace psimag {


  /**\ingroup symmetryConcepts
   *
   * \brief A class for implementing two-dimensional Identity Symmetry Elements.
   *
   */

  template<typename Field, typename Algorithms>
  class IdentityElement<Field,2,Algorithms>: 
    public SymmetryElement<Field,2,Algorithms> 
  {

  public:

    enum { DIM=2 };

    typedef IdentityElement<Field,DIM,Algorithms>     IdentityElementType;
    typedef SymmetryOperation<Field,DIM,Algorithms>   SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>     SymmetryElementType;
    typedef CellTranslation<Field,DIM>                CellTranslationType;
    typedef CartesianPosition<Field,DIM>              CartesianPositionType;
    typedef LatticeCoordinates<DIM>                   LatticeCoordinatesType;
 
    //======================================================================
    /**
     * \brief The default constructor produces 
     */
    IdentityElement():
      SymmetryElementType("identity",
			  latticeCoordinate             (0,0),   // => the c-axis direction
			  CellPosition<Field,DIM,Algorithms>::make(Field(0.0),Field(0.0)),   // no offset
			  convert<Field>("0") )           // no glide
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
    static SymmetryOperationType operationFor(const SymmetryElementType&           element,
					      const Lattice<Field,DIM,Algorithms>& lattice) 
    {
      static  SymmetryOperationType result = symmetryOperation<Field,int,Algorithms>(1, 0, 0, 1);
      result.name    = element.name;
      
      return result;
    }

  };

}
#endif 
/*@}*/

