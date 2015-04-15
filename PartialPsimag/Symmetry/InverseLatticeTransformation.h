//-*-C++-*-

#ifndef PSIMAG_InverseLatticeTransformation_H
#define PSIMAG_InverseLatticeTransformation_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file InverseLatticeTransformation.h Class definition for InverseLatticeTransformation. */

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
#include "LatticeTransformation.h"
#include "AppliedSymmetryElement.h"

#include "PSIMAGAssert.h"

namespace psimag {

  template <typename Field, size_t DIM> class LatticeTransformation;

  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class AppliedSymmetryElement;
  
  /** \ingroup latticeConcepts
   *
   * \brief Closure containing a reference to a LatticeTransformation, with the appropriate methods.
   *
   */
  template <typename Field, size_t DIM>
  class InverseLatticeTransformation {

  public:
    typedef LatticeTransformation<Field,DIM>        LatticeTransformationType;
    typedef SeitzMatrix<Field,DIM>                  MatrixType;
    typedef CellTranslation<Field,DIM>              CellTranslationType;
    typedef LatticeCoordinates<DIM>                 LatticeCoordinatesType;

    const LatticeTransformationType& transformation;

    /* \brief constructor
     */
    InverseLatticeTransformation(const LatticeTransformationType& latticeTransformation): 
      transformation(latticeTransformation) 
    {}

    //====================================================================== s

    //---------------------------------------------------------------------- LatticeCoordinatesType
    /** 
     * \brief Return the CellPosition in the transformed cell given
     *        the equivalent position in the origional cell.
     *
     */
    LatticeCoordinatesType operator()(const LatticeCoordinatesType& latticeCoordinates) const 
    {
      LatticeCoordinatesType  result;
      Multiply(transformation.cellNewToOrig, latticeCoordinates, result);
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
      Multiply(transformation.cellNewToOrig, cellTranslation, result);
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
      Multiply(transformation.cellNewToOrig, cellPosition, result);
      return result;
    }
    
    //---------------------------------------------------------------------- LatticeTemplate
    /** 
     * \brief This is the operation that makes the given tgtLattice a
     *        transformed version of the given srcLattice.
     *        &*&*&* check this!
     */
    template<typename Algorithms, 
	     template<typename,size_t,typename> class InLatticeTemplate,
	     template<typename,size_t,typename> class OutLatticeTemplate>
    void operator() (const InLatticeTemplate<Field,DIM,Algorithms>& srcLattice, 
		     OutLatticeTemplate<Field,DIM,Algorithms>& tgtLattice) const
    {      
      //      MatrixType invTransform;
      //      Inverse(transform,invTransform);
      const MatrixType& invTransform = transformation.cellOrigToNew;
      
      const MatrixType& rhs(invTransform);
      const MatrixType& lhs(srcLattice);
      MatrixType&       out(tgtLattice);
      
      Multiply(lhs,rhs,out);
      //      tgtLattice.transformation = *this;
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
      
      Multiply(symopM,                        transformation.cellOrigToNew, temp);
      Multiply(transformation.cellNewToOrig,  temp,                         tgtsymop);
      
      tgtsymop.name              = symop.name;
    }
    
#include "LatticeTransformationBase.h"

};  

} /* namespace psimag */

#endif //PSIMAG_LatticeTransformation_H

/*@}*/

