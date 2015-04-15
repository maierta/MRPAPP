//-*-C++-*-

#ifndef TESTPATTERN_H
#define TESTPATTERN_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file TestPattern.h
 *  Contains the class definition for TestPattern objects.
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
#include "Mat.h"
#include "PSIMAGAssert.h"
#include "Lattice.h"
#include "CellPosition.h"
#include "CartesianPosition.h"
#include "CartesianTranslation.h"
#include "LatticeTransformation.h"
#include "LatticeWithPattern.h"
#include "SymmetryOperation.h"

namespace psimag {

  template<typename Field, size_t DIM, typename> class Lattice;
  template<typename Field, 
	   size_t   DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class LatticeWithPattern;
  
  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename Algorithms>
  class Pattern;

  /** \ingroup patternConcepts
   *
   * \brief Template class for TestPatterns, which are patterns that
   * are temporarily created to test to see if a symmetry operation is
   * satisfied.
   *
   * \note Occupants must be comparable
   *  
   */
  template<typename Field,
	   size_t   DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class TestPattern: public LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>   {

  public:

    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>              TestPatternType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                             SymmetryOperationType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                              PatternType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>       LatticeWithPatternType;
    typedef CartesianTranslation<Field,DIM>                                     CartesianTranslationType;

    //======================================================================
    //
    /**
     * \brief Default Construct a TestPattern 
     **/
    TestPattern():
      LatticeWithPatternType()
    {}

    //======================================================================
    //
    /**
     * \brief Copy Construct a TestPattern 
     **/
    template<typename OtherLatticeType>
    TestPattern(const TestPattern<Field,DIM,Occupant,OtherLatticeType,Algorithms>&   testPattern):
      LatticeWithPatternType(testPattern)
    {}

    //======================================================================
    //
    /**
     * \brief Construct Pattern from another Pattern and a Cartesian
     *        SymmetryOperation. 
     *
     * \param pattern: This is the origional, untransformed, pattern.
     *
     * \param cartSymmetryOperation: This is the cartesian symmetry
     *                 operation that is being used to create this
     *                 pattern from the given pattern.
     *
     * \note The reference lattice for this pattern is the same aas
     *       the reference lattice for the given pattern.
     *
     * \note The constructed object shares the same Occupant objects
     *       with the given pattern.
     **/
    TestPattern(const LatticeWithPatternType&   latticeWithPattern,  
		const SymmetryOperationType&    cartSymmetryOperation):
      LatticeWithPatternType(latticeWithPattern)
    {
      const PatternType& pattern = latticeWithPattern.pattern;

      for (size_t pos=0; pos < pattern.NUMPOS; pos++)
	Multiply(cartSymmetryOperation,
		 pattern.cartesianPositions[pos], 
		 this->pattern.cartesianPositions[pos]);
      
      // Use the new cart. positions to generate normalized cell positions.
      this->generateCellPositions();
      
      // Move the positions to their equiv. positions with in the cell.
      this->pattern.normalizeCellPositions();
      
      // Re-generate the cartesian positions within the cell from the
      // normalized cell positions
      this->generateCartesianPositions();
    }

    /** 
     * Compute the mag squared of the difference between a position in
     * this pattern, given by pos1, and a position in otherPattern
     * given by pos2.
     *
     **/
    Field mag2(size_t pos1, size_t pos2, const PatternType& otherPattern) const {
      CartesianTranslationType diff = this->pattern.cartesianPositions[pos1] - otherPattern.cartesianPositions[pos2];
      return diff.length();
    }
    
    /** 
     * Compute the minimum mag squared difference between a position
     * in this pattern given by pos1 and any of the positions in
     * otherPattern given by the indicated range of indicies.
     *
     **/
    Field minDistance(size_t pos1, 
		      size_t pos2Start, size_t pos2End, 
		      const PatternType& otherPattern) const {
      
      Field min_ = mag2(pos1,pos2Start,otherPattern);
      for (size_t pos2=1; pos2 < pos2End; pos2++) {
	Field d = mag2(pos1,pos2,otherPattern);
	if (d<min_) min_ = d;
      }
      return min_;
    }
    
    /** 
     * Compute the 'distance' to the given pattern for the specified
     * occupant.
     *
     **/
    Field occupantDistanceTo(size_t occIndx, const PatternType& otherPattern) const {

      size_t start = this->pattern.occupantStart[occIndx];
      size_t end   = start + this->pattern.occupantNumPositions[occIndx];
      Field maxmin(0);

      for(size_t o=start; o<end; o++) {
	Field min = minDistance(o,start,end,otherPattern);
	if (min > maxmin) maxmin = min;
      }
      return maxmin;
    }

    /** 
     * Compute the 'distance' to the given pattern.
     *
     **/
    Field distanceTo(const PatternType& otherPattern) const{
      
      Field maxd(0);
      
      for(size_t o=0; o<this->pattern.occupants.size(); o++) {
	Field d = this->occupantDistanceTo(o,otherPattern);
	if (d > maxd) maxd = d;
      }
      return maxd;
    }

    /** 
     * Compute the 'distance' to the given pattern.
     *
     **/
    Field distanceTo(const SymmetryOperationType& symOp) const{
      PatternType otherPattern((*this),symOp);
      return this->distanceTo(otherPattern);
    }

  };
  //====================================================================== 
  
  /** \ingroup XML
   * Crystall Base XML output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms>
  Tag toXML(const TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>& latticeWithPattern,
	    std::string name="TestPattern") {
    
    const LatticeType&                            lattice = latticeWithPattern;
    const Pattern<Field,DIM,Occupant,Algorithms>& pattern = latticeWithPattern.pattern;

    Tag tag = toXML( lattice, name );
    
    Tag baseTag("IsA");
    baseTag["type"] = "LatticeWithPattern";
    tag.add(baseTag);
    
    tag.add(toXML(pattern,"Pattern",false));

    return tag;
  }
    
}  // end psimag  

#endif

/*@}*/
