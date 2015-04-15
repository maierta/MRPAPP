//-*-C++-*-

#ifndef ORIGIN_LOCATOR_H
#define ORIGIN_LOCATOR_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file OriginLocator.h
 *  Contains the class definition for OriginLocator objects.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <ostream>
#include <map>
#include <algorithm>
#include <functional>

#include "Vec.h"
#include "MAt.h"
#include "PSIMAGAssert.h"
#include "Lattice.h"
#include "CellPosition.h"
#include "CartesianPosition.h"
#include "CartesianTranslation.h"
#include "PatternData.h"
#include "PatternBase.h"
#include "TestPattern.h"
#include "OccupantClosure.h"
#include "LatticeTransformation.h"
#include "SymmetryOperation.h"

namespace psimag {

  template<typename Field, size_t DIM, typename> class Lattice;

  template<typename Field, size_t DIM, size_t NUMPOS, typename Occupant> class OccupantClosure;

  /** \ingroup patternConcepts
   *
   * \brief Template class for Pattern.
   *
   * \note Occupants must be comparable
   *  
   */
  template<typename Field, size_t DIM, size_t NUMPOS, 
	   typename Occupant, 
	   typename Algorithms>
  class OriginLocator: public PatternBase<Field,DIM,NUMPOS,Occupant,Algorithms>  {
    
  public:
    
    enum { PositionMatCols= NUMPOS, PositionMatRows=DIM+1};

    typedef Pattern<Field,DIM,NUMPOS,Occupant,Algorithms>     PatternType;
    typedef PatternBase<Field,DIM,NUMPOS,Occupant,Algorithms> PatternBaseType;

    typedef OccupantClosure<Field,DIM,NUMPOS,Occupant>    OccupantClosureType;
    typedef SymmetryOperation<Field,DIM>                  SymmetryOperationType;

    typedef PatternData<Field,DIM,NUMPOS,Occupant>        PatternDataType;
    typedef typename PatternDataType::CellPositionType    CellPositionType;
    typedef typename PatternDataType::CellPositions       CellPositions;
    typedef typename PatternDataType::OccItr              OccRef;
    typedef typename PatternDataType::CellPosItr          CellPosRef;

  public:
    
    /** 
     *  The object used to collect the pattern data while constructing
     *  the object, prior to finialization. 
     */
    PatternDataType         patternData;
    
    //============================================================
    /** \brief  Construct an empty pattern Pattern */
    Pattern(LatticeType& lattice): 
      PatternBaseType(lattice)
    {

    }
    
    //============================================================
    //============================================================

    /**
     * Return a OccupantClosure object that can be used with the += operator. 
     *
     * This is part of making: 
     *     pattern[occupant] += cellPosition 
     * work.
     *
     * \note OccupantClosure is a left hand side object supporting the
     *       += operator. This operator adds the cell position to
     *       occupant's cellPositions via a call back to patternData's
     *       insert member.
     */
    OccupantClosureType operator [] (const Occupant& occupant) {
      return OccupantClosureType(occupant, patternData);
    }
    
    //============================================================

    /**
     * Use the patternData to initialize the occupant data in this pattern.
     *
     * The Occupant cellPositions are placed in this->positions and
     * the occupant is noted in the coreesponding index of
     * this->occupant. The Occupants themselves are recorded in
     * this->occupants.
     *
     * \note The Occupant positions are the positions [DIM,PositionMatRows).
     */
    void loadCellPositions() {

      size_t     posIndex = DIM;
      size_t     occIndex = 0;
      OccRef     o;
      CellPosRef ocp;

      for(o  = patternData.occupants.begin(); 
	  o != patternData.occupants.end(); 
	  o++) {

	occupants.push_back(*o);
	occupantStart.push_back(posIndex);
	CellPositions& dataPositions    = patternData.occupantPositions[*o]; 
	size_t         numPos           = dataPositions.size();
	occupantNumPositions.push_back(numPos);
	
	for(ocp=dataPositions.begin(); ocp != dataPositions.end(); ocp++) {
	  
	  for(size_t i=0; i<PositionMatRows; i++)
	    positions(i,posIndex) = (*ocp)[i];
	  
	  occupant[posIndex]      = occIndex;
	  
	  posIndex++;
	}
	occIndex++;
      }
    }

    //============================================================

    void finalize() {

      if (patternData.finalized) return;

      loadCellPositions();

      normalizeOccupantPositions();
      
      generateCartesianPositions();

      patternData.finalized = true;
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
      TestPattern<Field,DIM,NUMPOS,Occupant,Algorithms> testPattern(*this,symmetryOperation);
      return testPattern.distanceTo(*this);
    }
  
    //============================================================

    /**
     * \brief Returns a vector that maps a CellPosition index into the
     *        index of the CellPosition produced by applying the given
     *        symmetry operations.
     *
     * Not right rework! &*&*&*
     */
    std::vector<int> positionMap(SymmetryOperationType& symOp) {

      std::vector<size_t> result(this->size());

      for(size_t i=0; i < positions.size(); i++) 
	result[i] = indexOf(symOp * positions[i]);
    
      return result;
    }

  };


  //============================================================

  /** \ingroup ostream
   *
   * \brief the Pattern ostream operator
   *
   */
  template<typename Field, typename Occupant, size_t DIM, size_t NUMPOS,typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const Pattern<Field,DIM,NUMPOS,Occupant,Algorithms>& pat) {
    //os << pat.patternData << std::endl;
    os << ((Patternbase) pat) << std::endl;
    return os;
  }
}  

#endif

/*@}*/
