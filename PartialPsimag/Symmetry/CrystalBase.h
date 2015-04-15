//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file CrystalBase.h  
 *
 *  Contains a the CrystalBase class.
 */

#ifndef Psimag_CrystalBase_H
#define Psimag_CrystalBase_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"
#include "Matrix.h"

//#include "Lattice.h"
#include "LatticeWithPattern.h"
#include "TestPattern.h"
#include "CellPosition.h"
#include "SpaceGroup.h"
#include "SpaceGroup2D.h"
#include "Symmetry.h"

namespace psimag {


  /** \ingroup crystallography
   *  
   *\brief Base class for the Crystal classes.
   *
   *\param DIM Dimensionality of the structure being
   *           represented. Legal values are 1,2,3.
   */ 
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class CrystalBase: 
    public LatticeWithPattern<Field,
			      DIM,
			      Occupant,
			      LatticeTemplate<Field,DIM,Algorithms>,
			      Algorithms>
  {
  public:

    //====================================================================== typedefs

    typedef CrystalBase<Field,DIM,Occupant,LatticeTemplate,Algorithms>    ThisType;

    typedef LatticeTemplate<Field,DIM,Algorithms>                         LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;
    typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>       SymmetryType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                        PatternType;
    typedef PatternData<Field,DIM,Occupant,Algorithms>                    PatternDataType;

    SymmetryType symmetry;

    // typedef SpaceGroup<Field,DIM,LatticeType,Algorithms>                  SpaceGroupType;
    // typedef std::vector<SpaceGroupType>                                   SpaceGroupVectorType;
    // typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>        TestPatternType;
    // typedef CellPosition<Field,DIM,Algorithms>                            CellPositionType;
    // typedef SymmetryOperation<Field,DIM,Algorithms>                       SymmetryOperationType;

//     enum { NUMGROUPS= SpaceGroupType::NUMGROUPS};

//     //====================================================================== members

//     SpaceGroupVectorType      spaceGroups;

//     size_t                    selectedSpaceGroupIndex;

//     std::vector<Field>        spaceGroupRatings;

//     bool                      spaceGroupsSet;

  public:
    
    //====================================================================== Constructors
    /**
     * \brief Construct a Crystal with a default lattice and an empty pattern
     *
     */
    CrystalBase(): 
      LatticeWithPatternType(),
      symmetry()
    {}

    /**
     * \brief Construct a Crystal with the given lattice and a trivial
     *        pattern with the given occupant at the origin.
     *
     */
    CrystalBase(const LatticeWithPatternType& latticeWithPattern):
      LatticeWithPatternType(latticeWithPattern), 
      symmetry(latticeWithPattern)
    { 
    }
    
//     /**
//      * \brief Construct a Crystal with the given lattice and an empty pattern
//      *
//      */
//     CrystalBase(LatticeType& lat):
//       LatticeWithPatternType(lat), 
//       spaceGroups(NUMGROUPS+1),
//       selectedSpaceGroupIndex(1),
//       spaceGroupRatings(NUMGROUPS+1),
//       spaceGroupsSet(false)
//     { 

//     }
    
    /**
     * \brief Construct a Crystal with the given lattice and a trivial
     *        pattern with the given occupant at the origin.
     *
     */
    CrystalBase(LatticeType& lat, const Occupant& occupant):
      LatticeWithPatternType(lat, occupant), 
      symmetry(lat)
    { 

    }
    
    /**
     * \brief Construct a Crystal with the given lattice and a trivial
     *        pattern with the given occupant at the origin.
     *
     */
    CrystalBase(const LatticeType& lat, const PatternDataType& patternData):
      LatticeWithPatternType(lat,patternData), 
      symmetry(*this)
    { 

    }
    
    /**
     * \brief Construct a Crystal with a lattice and pattern of a
     *        different lattice type.
     *
     */
    CrystalBase(const LatticeType& lat, const PatternType& pat):
      LatticeWithPatternType(lat,pat), 
      symmetry(*this)
    { 
    }

    /**
     * \brief Construct a Crystal with a lattice and pattern of a
     *        different lattice type.
     *
     */
    template<template<typename, size_t, typename> class OtherLatticeTemplate>
    CrystalBase(const OtherLatticeTemplate<Field,DIM,Algorithms>& lat, 
		const PatternType& pat):
      LatticeWithPatternType(lat,pat), 
      symmetry(*this)
    { 
    }

    /**
     * \brief Copy Construct a Crystal from another Crystal (with
     *        possibly a lattice and pattern of a different lattice
     *        type.)
     *
     */
    template<template<typename, size_t, typename> class OtherLatticeTemplate>
    CrystalBase(const CrystalBase<Field,DIM,Occupant,OtherLatticeTemplate,Algorithms>& other):
      LatticeWithPatternType(other), 
      symmetry(other)            // Don't copy the Symmetry chances are this does not make sense
    {}

    //======================================================================

//     void analyzeSpaceGroups() {
//       this->analyzeSpaceGroups(*this);
//     }

//     /** 
//      * Return the best spaceGroup for this Crystal
//      */ 
//     const SpaceGroupType& spaceGroup() const {
//       if (!spaceGroupsSet) {
// 	std::ostringstream buff;
// 	buff << "CrystalBase::spaceGroup(): spaceGroups have not been set in:" << std::endl;
// 	buff << (*this) << std::endl;
// 	ASSERT(spaceGroupsSet,
// 	       std::logic_error(buff.str()));
//       }

//       const SpaceGroupType& result = spaceGroups[selectedSpaceGroupIndex];

//       return result;
//     }

//     /** 
//      * \brief Set the selectedSpaceGroupIndex to the space group which
//      *        has a close match to this crystal and which has the
//      *        largest number of symmetry operations.
//      */ 
//     void selectSpaceGroup() {
//       ASSERT(spaceGroupsSet,
// 	     std::logic_error("CrystalBase::selectSpaceGroup: spaceGroups have not been set."));
//       static Field zero(0);
//       size_t maxNumOps = 1;
//       selectedSpaceGroupIndex = 1;

//       for (size_t i=1; i< NUMGROUPS + 1; i++) {
// 	if (!Algorithms::close(spaceGroupRatings[i], zero)) continue;
// 	size_t numOps = spaceGroups[i].numOperations;

// 	if (numOps > maxNumOps) {
// 	  maxNumOps = numOps;
// 	  selectedSpaceGroupIndex = i;
// 	}
//       }

//     }
  
};

//   //======================================================================
//   /** \ingroup XML
//    *
//    * \brief the Pattern XML outputer
//    *
//    */
//   template<typename Field, size_t DIM, typename Occupant, 
// 	   template<typename, size_t, typename> class LatticeTemplate,
// 	   typename Algorithms>
//   Tag movementTableXML(const CrystalBase<Field,DIM,Occupant,LatticeTemplate,Algorithms>& crystalBase) 
//   {
//     typedef LatticeTemplate<Field,DIM,Algorithms>                         LatticeType;
//     typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;

//     Matrix<int>                   movementIndex;
//     const LatticeWithPatternType& latpat = crystalBase;
//     crystalBase.symmetry.buildMovementIndex(latpat,movementIndex);

//     Tag movementTable("MovementTable");
//     movementTable["type"] = "Movement";

//     std::ostringstream buff;
//     movementIndex.print(buff);
//     movementTable.content << std::endl << buff.str();

//     return movementTable;
//   }
    
  //======================================================================

  /** \ingroup xml
   * Crystall Base xml constructor
   **/
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag toXML(const CrystalBase<Field,DIM,Occupant,LatticeTemplate,Algorithms>& crystalBase,
	    std::string name="CrystalBase") {
    
    typedef LatticeTemplate<Field,DIM,Algorithms>                            LatticeType;
    typedef LatticeWithPattern<Field, DIM, Occupant,LatticeType, Algorithms> BaseType;
    typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>          SymmetryType;

    const LatticeType&  lattice            = crystalBase;
    const BaseType&     latticeWithPattern = crystalBase;
    const SymmetryType& symmtery           = crystalBase.symmetry;

    Tag tag(name);
    tag.add(toXML( lattice ));
    tag.add(toXML( latticeWithPattern.pattern ));
    tag.add(toXML( symmtery ));
    //   tag.add(movementTableXML(crystalBase));

    return tag;
  }

  //======================================================================

  /** \ingroup ostream
   * Crystall Base output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const CrystalBase<Field,DIM,Occupant,LatticeTemplate,Algorithms>& crystalBase) {

    typedef LatticeTemplate<Field,DIM,Algorithms>                         LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;

    os << ( (LatticeWithPatternType) crystalBase ) << std::endl << std::endl;

    if (crystalBase.spaceGroupsSet) {
      os << " SpaceGroup["<< crystalBase.selectedSpaceGroupIndex << "]: " << crystalBase.spaceGroup().nam << std::endl;
      os << "   number of ops: " << crystalBase.spaceGroup().numOperations << std::endl;
      os << " spaceGroupRatings: " << std::endl;
      for (size_t i=0; i< crystalBase.NUMGROUPS; i++ )
	os << crystalBase.spaceGroups[i].nam << "  => " << crystalBase.spaceGroupRatings[i] << std::endl;

    } else {
      os << "SpaceGroups not set!" << std::endl;
    }

    return os;

  }

} /** namespace psimag */

#endif
/*@}*/
