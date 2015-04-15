//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file LatticeWithPattern.h  
 *
 *  Contains a the LatticeWithPattern class.
 */

#ifndef Psimag_LatticeWithPattern_H
#define Psimag_LatticeWithPattern_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "Tag.h"
#include "PSIMAGAssert.h"
#include "Matrix.h"

//#include "Lattice.h"
#include "PatternData.h"
#include "TestPattern.h"
#include "CellPosition.h"
//#include "SpaceGroup.h"
//#include "SpaceGroup2D.h"

namespace psimag {


  /** \ingroup crystallography
   *  
   *\brief Base class for the Crystal classes.
   *
   *\param DIM Dimensionality of the structure being
   *           represented. Legal values are 1,2,3.
   */ 
  template<typename Field, 
	   size_t   DIM, 
	   typename Occupant,
	   typename LatticeType,
	   typename Algorithms>
  class LatticeWithPattern: 
    public LatticeType
  {
  public:

    //====================================================================== typedefs

    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>    ThisType;

    typedef Pattern<Field,DIM,Occupant,Algorithms>                        PatternType;
    typedef PatternData<Field,DIM,Occupant,Algorithms>                    PatternDataType;

    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>        TestPatternType;
    typedef CellPosition<Field,DIM,Algorithms>                            CellPositionType;
    //    typedef SpaceGroup<Field, DIM, LatticeType, Algorithms>               SpaceGroupType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                       SymmetryOperationType;

    //====================================================================== members

    PatternType    pattern;

  public:
    
    //====================================================================== Constructors
    /**
     * \brief Construct a Crystal with a default lattice and an empty pattern
     *
     */
    LatticeWithPattern(): 
      LatticeType(),
      pattern()
    {
      generateCartesianPositions();
    }

    /**
     * \brief Copy Construct a Crystal.
     *
     */
    LatticeWithPattern(const LatticeWithPattern& latticeWithPattern):
      LatticeType(latticeWithPattern), 
      pattern(latticeWithPattern.pattern)
    { 
      generateCartesianPositions();
    }
    
    /**
     * \brief Construct a Crystal with the given lattice and a trivial
     *        pattern with the given occupant at the origin.
     *
     */
    LatticeWithPattern(LatticeType& lat, const Occupant& occupant):
      LatticeType(lat), 
      pattern(occupant)
    { 
      generateCartesianPositions();
    }
    
    /**
     * \brief Construct a Crystal with the given lattice and a trivial
     *        pattern with the given occupant at the origin.
     *
     */
    LatticeWithPattern(const LatticeType& lat, const PatternDataType& patternData):
      LatticeType(lat), 
      pattern(patternData)
    { 
      generateCartesianPositions();
    }
    
    /**
     * \brief Construct a Crystal with a lattice and pattern of a
     *        different lattice type.
     *
     */
    LatticeWithPattern(const LatticeType& lat, const PatternType& pat):
      LatticeType(lat), 
      pattern(pat)
    { 
      generateCartesianPositions();
    }

    /**
     * \brief Construct a Crystal with a lattice and pattern of a
     *        different lattice type.
     *
     */
    template<template<typename, size_t, typename> class OtherLatticeTemplate>
    LatticeWithPattern(const OtherLatticeTemplate<Field,DIM,Algorithms>& lat, 
		       const PatternType& pat):
      LatticeType(lat), 
      pattern(pat)
    { 
      generateCartesianPositions();
    }

    /**
     * \brief Copy Construct a Crystal from another Crystal (with
     *        possibly a lattice and pattern of a different lattice
     *        type.)
     *
     */
    template<typename OtherLatticeType>
    LatticeWithPattern(const LatticeWithPattern<Field,DIM,Occupant,OtherLatticeType,Algorithms>& other):
      LatticeType(other), 
      pattern(other.pattern)
    { 
      //need to copy the  spacegroup info!
      generateCartesianPositions();
    }

    //============================================================ Methods 

    /** 
     * \brief Use this patterns cart. positions and reference lattice
     *        to generate cell positions from this patterns cartesian
     *        positions.
     **/
    void generateCellPositions()  {
      for(size_t pos=0; pos < pattern.NUMPOS; pos++)
	Multiply(this->inverse(),
		 pattern.cartesianPositions[pos],
		 pattern.cellPositions[pos]);
    }
    
    //============================================================

    /** 
     * \brief Regenerate the cartesian positions within the cell from
     *        the normalized cell positions.
     **/
    void generateCartesianPositions() {
      for(size_t pos=0; pos < pattern.NUMPOS; pos++)
	Multiply((*this), 
		 pattern.cellPositions[pos],
		 pattern.cartesianPositions[pos]);
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
      TestPatternType testPattern((*this),symmetryOperation);
      return testPattern.distanceTo(pattern);
    }


};

  //====================================================================== 
  
  /** \ingroup XML
   * Crystall Base XML output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms>
  Tag toXML(const LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>& latticeWithPattern,
	    std::string name="LatticeWithPattern",bool doTables=true) {
    
    const LatticeType&                            lattice = latticeWithPattern;
    const Pattern<Field,DIM,Occupant,Algorithms>& pattern = latticeWithPattern.pattern;

    Tag tag = toXML( lattice, name );
    
    Tag baseTag("IsA");
    baseTag["type"] = "LatticeWithPattern";
    tag.add(baseTag);
    
    tag.add(toXML(pattern,"Pattern",doTables));

    return tag;
  }
    
  //======================================================================

  /** \ingroup ostream
   * Crystall Base output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>& latticeWithPattern) {

    os << ( (LatticeType) latticeWithPattern ) << std::endl << std::endl;
    os << " --------------------------------------------- Pattern:" << std::endl;
    os << latticeWithPattern.pattern << std::endl;
    os << "---------------------------------------------- End Pattern" << std::endl;

    return os;

  }

} /** namespace psimag */

#endif
/*@}*/
