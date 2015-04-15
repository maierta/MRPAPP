//-*-C++-*-

#ifndef PSIMAG_PatternWithLattice_H
#define PSIMAG_PatternWithLattice_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file PatternWithLattice.h
 *  Contains the class definition for PatternWithLattice objects.
 */

#include <iostream>                                  
#include <ostream>
#include <sstream>                                  
#include <iomanip>

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <functional> 

#include "Vec.h"
#include "Mat.h"

#include "PSIMAGAssert.h"
#include "Lattice.h"
#include "CellPosition.h"
#include "CartesianPosition.h"
#include "CellTranslation.h"
#include "CartesianTranslation.h"
#include "SymmetryOperation.h"
#include "Pattern.h"
#include "TestPattern.h"
#include "PatternData.h"
#include "LatticeTransformation.h"

namespace psimag {

  template<typename Field, size_t DIM, typename> class Lattice;
  template<typename Field, size_t DIM>           class LatticeTransformation;

  /** \ingroup patternConcepts
   *
   * \brief Template class for Patternbase, inherited by Pattern and TestPattern
   *
   * \note Occupants must be comparable
   *  
   */
  template<typename Field, 
	   size_t   DIM,
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class PatternWithLattice: 
    public Pattern<Field,DIM,Occupant,Algorithms>  
  {
  public:

    enum { NROWS=DIM+1};
    
    //====================================================================== typedefs

    typedef Pattern<Field,DIM,Occupant,Algorithms>                 PatternBaseType;
    typedef PatternData<Field,DIM,Occupant,Algorithms>             PatternDataType;
    typedef LatticeTransformation<Field,DIM>                       LatticeTransformationType;
    typedef SymmetryOperation<Field,DIM>                           SymmetryOperationType;

    //====================================================================== members

    /** 
     *  \brief A reference to the Lattice that defines the cell
     *         geometry and cell positions for this pattern.
     */
    const LatticeType&           referenceLattice;
    
    //============================================================ Constructors
    //============================================================

    /** \brief  Construct an PatternWithLattice with an empty pattern.
     */
    PatternWithLattice(const LatticeType& lattice): 
      PatternBaseType(),
      referenceLattice(lattice)
    {

    }

    /** \brief  Construct an PatternWithLattice with an empty pattern.
     */
    PatternWithLattice(const LatticeType& lattice, const Occupant& occupant): 
      PatternBaseType(occupant),
      referenceLattice(lattice)
    {

    }

    /** \brief  Construct an PatternWithLattice with an empty pattern.
     */
    PatternWithLattice(const LatticeType& lattice, PatternDataType& patternData): 
      PatternBaseType(patternData),
      referenceLattice(lattice)
    {
      this->normalizeCellPositions();
      this->generateCartesianPositions();

    }

    /** \brief  Copy Construct a PatternWithLattice. */
    PatternWithLattice(const PatternWithLattice& pat): 
      PatternBaseType(pat),
      referenceLattice(pat.referenceLattice)
    {
    }

    /** \brief  Copy Construct a PatternBase. 
     *
     * \note This is for constructing ReducedPatterns from Patterns of other lattice types.
     **/
    template<typename OtherLatticeType>
    PatternWithLattice(const Pattern<Field,DIM,Occupant,Algorithms>& pat, 
		       OtherLatticeType& lattice): 
      PatternBaseType(pat),
      referenceLattice(lattice)
    {

    }
    
    /** \brief  Construct a transformed Patternbase from the given pattern.  */
    PatternWithLattice(const PatternWithLattice& pat, const LatticeTransformationType& transform): 
      PatternBaseType(pat),
      referenceLattice(pat.referenceLattice)
    {
      transform(*this);

    }
    
    //============================================================ Methods relying on referenceLattice:

    /** 
     * \brief Use this patterns cart. positions and reference lattice
     *        to generate cell positions from this patterns cartesian
     *        positions.
     **/
    void generateCellPositions()  {
      for(size_t pos=0; pos < this->NUMPOS; pos++)
	Multiply(referenceLattice.inverse(),this->cartesianPositions[pos],this->cellPositions[pos]);
    }
    
    //============================================================

    /** 
     * \brief Regenerate the cartesian positions within the cell from
     *        the normalized cell positions.
     **/
    void generateCartesianPositions() {
      for(size_t pos=0; pos < this->NUMPOS; pos++)
	Multiply(referenceLattice,this->cellPositions[pos],this->cartesianPositions[pos]);
    }

    //============================================================

    /**
     * \brief Returns true if the given symmetry operations map this
     *        pattern onto itself.
     */
    bool satisfies(const SymmetryOperationType&  symmetryOperation) {
      if (differenceFor(symmetryOperation) == 0) return true;
      return false;
    }
    
    //============================================================
    
    /**
     * \brief Returns true if the given symmetry operations map this
     *        pattern onto itself.
     */
    bool differenceFor(const SymmetryOperationType&  symmetryOperation) {
      TestPattern<Field,DIM,Occupant,LatticeType,Algorithms> testPattern(*this,symmetryOperation);
      return testPattern.distanceTo(*this);
    }
  
    //============================================================ done elsewhere

//     /**
//      * \brief Returns a vector that maps a CellPosition index into the
//      *        index of the CellPosition produced by applying the given
//      *        symmetry operations.
//      *
//      * Not right rework! &*&*&*
//      */
//     std::vector<int> positionMap(SymmetryOperationType& symOp) {

//       std::vector<size_t> result(this->size());

//       for(size_t i=0; i < this->positions.size(); i++) 
// 	result[i] = indexOf(symOp * this->positions[i]);
    
//       return result;
//     }

  };
}  

#endif

/*@}*/
