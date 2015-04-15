//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file Crystal.h  
 *
 *  Contains a the Crystal class.
 */

#ifndef PSIMAG_CRYSTAL_H
#define PSIMAG_CRYSTAL_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"

#include "Lattice.h"
#include "ReducedLattice.h"
#include "Pattern.h"
#include "PatternData.h"
#include "LatticeWithPattern.h"
#include "SpaceGroup.h"
#include "CrystalBase.h"
#include "ReducedCrystal.h"
#include "ConventionalCrystal.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Algorithms> class ReducedLattice;  

  //======================================================================
  
  /** \ingroup crystallography
   *  
   * \brief Crystal class specialization. (This is the one the user calls.)
   *
   * \param DIM Dimensionality of the structure being
   *            represented. Legal values are 1,2,3.
   */ 
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  class Crystal: 
    public CrystalBase<Field, DIM,Occupant,Lattice,Algorithms>  
  {
  public:
    
    //====================================================================== typedefs

    typedef Crystal            <Field,DIM,Occupant,            Algorithms> ThisType;
    typedef CrystalBase        <Field,DIM,Occupant,Lattice,    Algorithms> BaseType;
    typedef Lattice            <Field,DIM,                     Algorithms> LatticeType;
    typedef ReducedLattice     <Field,DIM,                     Algorithms> ReducedLatticeType;
    typedef Pattern            <Field,DIM,Occupant,            Algorithms> PatternType;
    typedef PatternData        <Field,DIM,Occupant,            Algorithms> PatternDataType;
    typedef LatticeWithPattern <Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;
    typedef LatticeWithPattern <Field,DIM,Occupant,ReducedLatticeType,Algorithms> ReducedLatticeWithPatternType;
    typedef ReducedCrystal     <Field,DIM,Occupant,            Algorithms> ReducedCrystalType;
    typedef ConventionalCrystal<Field,DIM,Occupant,            Algorithms> ConventionalCrystalType;
    //typedef SpaceGroup         <Field,DIM,LatticeType,         Algorithms> SpaceGroupType;
    //typedef SpaceGroup         <Field,DIM,ReducedLatticeType,  Algorithms> ReducedSpaceGroupType;
    typedef SymmetryOperation  <Field,DIM,Algorithms>                      SymmetryOperationType;
    typedef Symmetry<Field,DIM,Occupant,Lattice,Algorithms>                SymmetryType;  
    typedef Symmetry<Field,DIM,Occupant,ReducedLattice,Algorithms>         ReducedSymmetryType;  
  
    //====================================================================== members

    ReducedCrystalType      reducedCrystal; 

    ConventionalCrystalType conventionalCrystal; 

    //====================================================================== constructors

    /* 
     * Construct a Crystal from a Lattice. 
     * The constructed pattern is empty.
     */
    Crystal(LatticeType& lat):
      BaseType(lat),
      reducedCrystal(*this), 
      conventionalCrystal(this->reducedCrystal)
    {
      setAppliedSymmetryElements(reducedCrystal.transform.inverse(), reducedCrystal);
      //conventionalCrystal.setAppliedSymmetryElements(reducedCrystal);
    }
    
    /* 
     * Construct a Crystal from a Lattice and a single Occupant. 
     * The pattern produced is a trivial pattern with the occupant at the origin.
     */
    Crystal(LatticeType& lat, const Occupant& occupant):
      BaseType(lat, occupant),
      reducedCrystal(*this), 
      conventionalCrystal(this->reducedCrystal)
    {
      setAppliedSymmetryElements(reducedCrystal.transform.inverse(), reducedCrystal);
      //conventionalCrystal.setAppliedSymmetryElements(reducedCrystal);
    }
    
    /* 
     * Construct a Crystal from a Lattice and a single Occupant. 
     * The pattern produced is a trivial pattern with the occupant at the origin.
     */
    Crystal(LatticeType lat, const Occupant occupant):
      BaseType(lat, occupant),
      reducedCrystal(*this), 
      conventionalCrystal(this->reducedCrystal)
    {
      setAppliedSymmetryElements(reducedCrystal.transform.inverse(), reducedCrystal);
      //conventionalCrystal.setAppliedSymmetryElements(reducedCrystal);
    }
    
    /* 
     * Construct a Crystal from a Lattice and a single Occupant. 
     * The pattern produced is a trivial pattern with the occupant at the origin.
     */
    Crystal(const LatticeType& lat, const PatternDataType& patternData):
      BaseType(lat, patternData),
      reducedCrystal(*this), 
      conventionalCrystal(this->reducedCrystal)
    {

      ReducedSymmetryType& reducedCrystalSymmetry = reducedCrystal.symmetry;
      this->symmetry.setAppliedSymmetryElements(reducedCrystal.transform.inverse(), reducedCrystalSymmetry);

    }

    /* 
     * Construct a Crystal from a Lattice and a Pattern.
     */
    Crystal(LatticeType& lat, PatternType& pat):
      BaseType(lat, pat),
      reducedCrystal(*this), 
      conventionalCrystal(this->reducedCrystal)
    {
      this->symmetry.setAppliedSymmetryElements(reducedCrystal.transform.inverse(), reducedCrystal);
      //conventionalCrystal.setSpaceGroups(reducedCrystal);
    }

    /* 
     * Construct a Crystal from a Lattice and a Pattern.
     */
    Crystal(LatticeWithPatternType& pat):
      BaseType(pat),
      reducedCrystal(*this), 
      conventionalCrystal(this->reducedCrystal)
    {

      this->symmetry.setAppliedSymmetryElements(reducedCrystal.transform.inverse(), reducedCrystal.symmetry);
      //conventionalCrystal.setSpaceGroups(reducedCrystal);
    }

    /* 
     * Copy Construct a Crystal from an other Crystal.
     *
     * \note The Lattice reduction and space group generation proceedures are not repeated.
     */
    Crystal(ThisType& otherStructure):
      BaseType(otherStructure),
      reducedCrystal(otherStructure.reducedCrystal), 
      conventionalCrystal(otherStructure.conventionalCrystal)
    { 
      // Don't set SpaceGroups here.
    }



  };

  //======================================================================

  /** \ingroup XML
   * Cell XML output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  Tag toXML(const Crystal<Field,DIM,Occupant,Algorithms>& crystal,
	    std::string name="Crystal") {

    typedef CrystalBase<Field, DIM,Occupant,Lattice,Algorithms> BaseType;

    Tag tag(name);

//     Tag baseTag("IsA");
//     baseTag.content << "Crystal";
//     tag.add(baseTag);

    tag.add(toXML( (BaseType)crystal, "InputCrystal"));
    tag.add(toXML(crystal.reducedCrystal,      "ReducedCrystal"));
    //tag.add(toXML(crystal.conventionalCrystal, "ConventionalCrystal"));
    return tag;
  }

  //======================================================================

  /** \ingroup ostream
   * Cell output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const Crystal<Field,DIM,Occupant,Algorithms>& crystalStruture) {

    typedef CrystalBase<Field, DIM,Occupant,Lattice,Algorithms> CrystalBaseType;

    os << std::endl;
    os << "========================================================== Crystal Structure" << std::endl;
    os << "==================================================== Crystal Structure" << std::endl;
    os << ( (CrystalBaseType) crystalStruture ) << std::endl << std::endl;
    os << crystalStruture.reducedCrystal << std::endl<< std::endl;
    os << crystalStruture.conventionalCrystal << std::endl << std::endl;
    os << "============================================================" << std::endl;

    
    return os;
  }


} /** namespace psimag */

#endif
/*@}*/
