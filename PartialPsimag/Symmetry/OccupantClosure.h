//-*-C++-*-

#ifndef OCCUPANT_CLOSURE_H
#define OCCUPANT_CLOSURE_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file OccupantClosure.h
 *  Contains the class definition for OccupantClosure objects.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <ostream>

#include "PSIMAGAssert.h"
#include "Pattern.h"
#include "PatternData.h"
#include "CartesianPosition.h"

namespace psimag {

  template<typename Field, 
	   size_t DIM,
	   typename Occupant, 
	   typename Algorithms>
  class Pattern;

  /** \ingroup patternConcepts
   *
   * \brief Template class for OccupantClosure.
   *
   * These object store a reference to a Pattern and an Occupant.
   * so that they can be passed used on the LHS in += operator.
   *
   * This enables expression like pattern[occupant] += cellPosition to
   * work as desired.
   */
  template<typename Field, 
	   size_t DIM, 
	   typename Occupant, 
	   typename LatticeType, 
	   template<typename Field,
		    size_t   DIM,
		    typename Occupant,
		    typename LatticeType,
		    typename Algorithms> PatternTemplate,
	   typename Algorithms>
  class OccupantClosure  {
  public:
    //====================================================================== typedefs

    typedef PatternTemplate<Field,DIM,Occupant,LatticeType,Algorithms> PatternType;
    typedef PatternData<Field,DIM,Occupant,Algorithms>                 PatternDataType;
    typedef CartesianPosition<Field,DIM>                               CartesianPositionType;
    typedef typename PatternDataType::CellPositionType                 CellPositionType;
    typedef typename PatternDataType::CellPositions                    CellPositions;
    
    //====================================================================== members
  public:

    const Occupant&   occupant;
    PatternType&      pattern;
    PatternDataType&  patternData;
    
    //====================================================================== Constructors

    OccupantClosure(const Occupant& occupant_, PatternType& pattern_): 
      occupant(occupant_), pattern(pattern_), patternData(pattern_.patternData) 
    {}
    
    //====================================================================== Methods

    OccupantClosure& operator+= (const CellPositionType& cellPosition) {
      patternData.insert(occupant,cellPosition);
      return (*this);
    }
    
    OccupantClosure& operator+= (const CartesianPositionType& cartPosition) {
      patternData.insert(occupant, pattern.referenceLattice().cellPosition(cartPosition));
      return (*this);
    }
    
    /** Provide a conversion operator to make this object look ike a CellPositions object.  */
    operator CellPositions() const {					\
      return patternData.tempPattern.occupantPositions[occupant]; 
    } 
  };
}

#endif

/*@}*/
