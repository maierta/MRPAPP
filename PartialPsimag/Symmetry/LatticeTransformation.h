//-*-C++-*-

#ifndef PSIMAG_LatticeTransformation_H
#define PSIMAG_LatticeTransformation_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file LatticeTransformation.h Class definition for LatticeTransformation. */

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Mat.h"
#include "Real.h"
#include "Lattice.h"
#include "CellPosition.h"
#include "LatticeWithPattern.h"
#include "SeitzMatrix.h"
#include "SymmetryOperation.h"
#include "SymmetryElement.h"
#include "AppliedSymmetryElement.h"
//#include "Crystal.h"
#include "PSIMAGAssert.h"
#include "InverseLatticeTransformation.h"

namespace psimag {

  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class LatticeWithPattern;

  template<typename Field, size_t DIM>
  class InverseLatticeTransformation;

  /** \ingroup latticeConcepts
   *
   * \brief A SeitzMatrix< Field, DIM > which can be used to performs
   *        mappings related to the change of basis from one lattice to
   *        another.
   *
   *  The process by which a Lattice is characterized and it's symmetry group
   *  identified involves a series of change-of-basis operations. These
   *  operations are represented by \b LatticeTransform objects which are
   *  (usually integer) matrices which transform a given set of lattice
   *  basis vectors into another set. The transpose of this matrix can be
   *  used to transform CellPosition objects from one cell to another.
   *
   *  LatticeTransformation objects transform:
   * 
   *  Lattices: Right multiplying the SeitzMatrix<Field,DIM> produces
   *            a matrix which is the basis of the transformed
   *            lattice. 
   *
   *            If the originonal and transformed Cells are primitive
   *            the transformation matrix must be an integer
   *            matrix. 
   *
   *            If the transformation is from a primitive to a
   *            conventional cell or from one conventional setting to
   *            another the matirix may contain select rationals in
   *            addition to integral values.
   *
   *  SymmetryOperations: Transform a SymmetryOperation (specified in
   *                      the coordinate system of one lattice) to the
   *                      same SymmetryOperation in the coordinate
   *                      system of the transformed lattice.
   *
   *  CellPositions: Transform a CellPosition given in the coordinate system
   *                 of one lattice to the same CellPosition in the
   *                 coordinate system of the transformed lattice).
   *
   *  Patterns: Produce a new Pattern in which all of the pattern's
   *            CellPositions are replaced by the transformed CellPositions.
   *
   *  Crystals: Produces a new Crystal which contains
   *                     a transformed cell, patternPattern, etc.
   *
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   *  \param DIM: Dimensionality of the lattice being
   *              represented. Legal values are 1,2,3.
   */
  template <typename Field, size_t DIM>
  class LatticeTransformation: 
    public SeitzMatrix< Field, DIM >
  {
  public:

    typedef SeitzMatrix<Field,DIM>                  MatrixType;
    typedef CellTranslation<Field,DIM>              CellTranslationType;
    typedef LatticeCoordinates<DIM>                 LatticeCoordinatesType;
    typedef InverseLatticeTransformation<Field,DIM> InverseLatticeTransformationType;

    MatrixType cellOrigToNew;
    MatrixType cellNewToOrig;

    /**
     * \brief Default constructor
     *
     */
    LatticeTransformation<Field,DIM>(): MatrixType() {}

    /**
     * \brief Construct a transform so that the only basis vector that
     *        is changed is the last one and it will be a combination
     *        of the origional basis vectors specified by bcoords.
     *
     * \note bcoords is also the coordinates of the new basis vector
     *       (i.e. the last one) in the origional coordinate system.
     */
    LatticeTransformation<Field,DIM>(const Vec<int,DIM>& bcoords): 
    MatrixType() 
    {
      for( size_t i=0 ; i < DIM ;i++)
	this->rotation(DIM-1, i) = bcoords[i];
      finalize();
    }

    /**
     * \brief Finalize this transform.
     *
     */
    void finalize()   {
      cellNewToOrig.rotation = this->rotation;
      //      Transpose(this->rotation,cellNewToOrig.rotation);
      //Copy(this->rotation,cellNewToOrig.rotation);
      Inverse(cellNewToOrig.rotation,cellOrigToNew.rotation);
    }

    /**
     * \brief Return a close of this object that acts like its inverse.
     *
     */
    InverseLatticeTransformationType inverse()   {
      InverseLatticeTransformationType result(*this);
      return result;
    }

    /**
     * \brief Matrix Element Accessor
     *
     */
    Field& operator () (const size_t& row, const size_t& col) { 
      MatrixType& m(*this);
      return m(row,col);
    }

    //====================================================================== 'Virtuals'
    //
    //---------------------------------------------------------------------- LatticeCoordinatesType
    /** 
     * \brief Return the CellPosition in the transformed cell given
     *        the equivalent position in the origional cell.
     *
     */
    LatticeCoordinatesType operator()(const LatticeCoordinatesType& latticeCoordinates) const
    {
      LatticeCoordinatesType  result;
      Multiply(cellOrigToNew, latticeCoordinates, result);
      return result;
    }

    //---------------------------------------------------------------------- CellTranslationType
    /** 
     * \brief Return the CellPosition in the transformed cell given
     *        the equivalent position in the origional cell.
     *
     */
    CellTranslationType operator()(const CellTranslationType& cellTranslation) const 
    {
      CellTranslationType  result;
      Multiply(cellOrigToNew, cellTranslation, result);
      return result;
    }
   
    //---------------------------------------------------------------------- CellPositionType
    /** 
     * \brief Return the CellPosition in the transformed cell given
     *        the equivalent position in the origional cell.
     *
     */
    template<typename Algorithms>  
    CellPosition<Field,DIM,Algorithms> operator()(const CellPosition<Field,DIM,Algorithms>& cellPosition) const 
    {
      CellPosition<Field,DIM,Algorithms> result;
      Multiply(cellOrigToNew, cellPosition, result);
      return result;
    }
  
    //---------------------------------------------------------------------- LatticeTemplate
    /** 
     * \brief This is the operation that makes the given tgtLattice a
     *        transformed version of the given srcLattice.
     *
     */
    template<typename Algorithms, 
	     template<typename,size_t,typename> class InLatticeTemplate,
	     template<typename,size_t,typename> class OutLatticeTemplate>
    void operator() (const InLatticeTemplate<Field,DIM,Algorithms>& srcLattice,
		     OutLatticeTemplate<Field,DIM,Algorithms>& tgtLattice) const 
    {
      const MatrixType& rhs(*this);
      const MatrixType& lhs(srcLattice);
      MatrixType&       out(tgtLattice);
      
      Multiply(lhs,rhs,out);
      tgtLattice.transform = *this;
      tgtLattice.update();
    }

    //---------------------------------------------------------------------- SymmetryOperationType
    /** 
     * \brief Transform the given symop into the other. 
     *        This is a Similarity transform.
     *
     */
    template<typename Algorithms>
    void operator()(const SymmetryOperation<Field,DIM,Algorithms>& symop,  
		    SymmetryOperation<Field,DIM,Algorithms>& tgtsymop) const
    {
      MatrixType              temp;
      const MatrixType&       symopM(symop);
      
      Multiply(symopM,         cellNewToOrig, temp);
      Multiply(cellOrigToNew,  temp,          tgtsymop);
      
      tgtsymop.name              = symop.name;
    }
    
    //====================================================================== 'Abstract Base'

#include "LatticeTransformationBase.h"

  };

  //======================================================================
  //======================================================================

  /** \ingroup xml
   *
   * XML Output function for 2D CellParameters.
   */
  template<typename Field, size_t DIM>
  Tag toXML(const LatticeTransformation<Field,DIM>& latTrans,
	    std::string name="LatticeTransformation") {
      
    typedef typename LatticeTransformation<Field,DIM>::MatrixType MatType; 
    typedef typename MatType::Traits                              MatTraitsType;

    Tag result(name);
    result.add(toXML(latTrans.cellOrigToNew,"cellOrigToNew"));
    result.add(toXML(latTrans.cellNewToOrig,"cellNewToOrig"));
    
    std::ostringstream buff;
    MAT_PRINT<MatType,MatTraitsType>::JUST_NUMBERS( latTrans, buff);
    result.content << buff.str();
    
    return result;
  }
  
  //======================================================================

  /** \ingroup ostream
   * Lattice output stream operator 
   **/
  template<typename Field, size_t DIM>
  std::ostream& operator << (std::ostream& os, 
			     const LatticeTransformation<Field,DIM>& ct) {

    typedef typename LatticeTransformation< Field, DIM >::MatrixType CTMatType;
    
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    
    os << " ----------------------------------- Lattice Transformation: " << std::endl;
    os << ( static_cast<CTMatType>(ct) ) << std::endl;
    os << " orig CellCoordinates to new CellCoordinates:" << std::endl;
    os << ( static_cast<CTMatType>(ct.cellOrigToNew) ) << std::endl;
    os << " new CellCoordinates to orig CellCoordinates:" << std::endl;
    os << ( static_cast<CTMatType>(ct.cellNewToOrig) ) << std::endl;
    os << "------------------------------------ End Lattice Transformation";
    
    return os;
  }

} /* namespace psimag */

#endif //PSIMAG_LatticeTransformation_H


/*@}*/

