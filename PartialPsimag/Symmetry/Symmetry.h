//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file Symmetry.h  
 *
 *  Contains a the Symmetry class.
 */

#ifndef Psimag_Symmetry_H
#define Psimag_Symmetry_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"
#include "Matrix.h"

//#include "Lattice.h"
#include "SymmetryGroup.h"
#include "ReducedLattice.h"
#include "LatticeWithPattern.h"
#include "InverseLatticeTransformation.h"
#include "TestPattern.h"
#include "CellPosition.h"
#include "SpaceGroup.h"
#include "SpaceGroup2D.h"
#include "SymmetryElements2D.h"
#include "AppliedSymmetryElement.h"

namespace psimag {

  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class AppliedSymmetryElement;

  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class Star;

  template<typename Field, size_t DIM,typename Algorithms>
  class ReducedLattice;

  /** \ingroup crystallography
   *  
   *\brief Base mixin class for associating spacegroups with a LatticeWithPattern.
   *
   */ 
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class Symmetry
  {
  public:

    //====================================================================== typedefs

    typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>       SymmetryType;

    typedef LatticeTemplate<Field,DIM,Algorithms>                         LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;

    typedef Pattern<Field,DIM,Occupant,Algorithms>                        PatternType;

    typedef SpaceGroup<Field,DIM,Occupant,Algorithms>                     SpaceGroupType;
    typedef std::vector<SpaceGroupType>                                   SpaceGroupVectorType;
    typedef CellPosition<Field,DIM,Algorithms>                            CellPositionType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                       SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms>                        SymmetryElementsType;
    typedef SymmetryGroup<Field,DIM,Occupant,LatticeTemplate,Algorithms>  SymmetryGroupType;

    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>            TestPatternType;

    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
    typedef std::vector<AppliedSymmetryElementType>                           AppliedSymmetryElementsType;


    enum { NUMGROUPS= SpaceGroupType::NUMGROUPS};

    //====================================================================== members

    /** \brief The lattice and pattern that the spacegroups symmetry operations refer to.
     *
     * \note A refernce to LatticeWithPatternType does not work.
     */
    //    LatticeWithPatternType      latpat;   // Can't be a reference?

    SpaceGroupVectorType        spaceGroups;

    size_t                      selectedSpaceGroupIndex;

    bool                        spaceGroupsSet;

    //The use of this is optional to the SymmetryElementsType specialization
    AppliedSymmetryElementsType allAppliedSymmetryElements;  

    AppliedSymmetryElementsType appliedSymmetryElements;

    const LatticeWithPatternType latpat;

    SymmetryGroupType            symmetryGroup;
    SymmetryGroupType            pointGroup;

  public:
    
    //====================================================================== Constructors
    /**
     * \brief Construct an un-set Symmetry
     *
     */
    Symmetry(const LatticeWithPatternType& latPatt):
      spaceGroups                (NUMGROUPS+1),
      selectedSpaceGroupIndex    (1),
      spaceGroupsSet             (false),
      allAppliedSymmetryElements (),
      appliedSymmetryElements    (),
      latpat                     (latPatt),  // copied, may be reduced
      symmetryGroup              (latpat),   // refer to the copy
      pointGroup                 (latpat)    // refer to the copy
    {}

    //======================================================================

    /**
     *\brief Use the SymmetryElementsType to set the appliedSymmetryElements field.
     *
     *\note This implies that the Symmetry is the Symmetry of a ReducedCrystal. Otherwise
     *      calling this method may not be appropriate.
     */
    void setAppliedSymmetryElements() {
      
      allAppliedSymmetryElements  = SymmetryElementsType::allAppliedElements(latpat); 
      appliedSymmetryElements     = SymmetryElementsType::validElements(latpat, allAppliedSymmetryElements);
      
      symmetryGroup.init(appliedSymmetryElements);

      pointGroup.init(symmetryGroup.pointGroupAppliedElements());

      analyzeSpaceGroups();
    }

    //======================================================================

    /**
     *\brief Use the appliedSymmetryElements of another Symmetry to
     *       set this Symmetries appliedSymmetryElements. It is
     *       assumed that the lattices of the two symmetries are
     *       related by the given LatticeTransform.
     *
     *\note This implies that the Symmetry is the NOT Symmetry of a ReducedCrystal. Otherwise
     *      calling this method may not be appropriate.
     */
    void setAppliedSymmetryElements(const InverseLatticeTransformation<Field,DIM>&          transform, 
				    Symmetry<Field,DIM,Occupant,ReducedLattice,Algorithms>& reducedSymmetry) {
      
      allAppliedSymmetryElements.resize(reducedSymmetry.allAppliedSymmetryElements.size());
      for (size_t i=0; i< reducedSymmetry.allAppliedSymmetryElements.size(); i++) 
 	transform(reducedSymmetry.allAppliedSymmetryElements[i], allAppliedSymmetryElements[i]);
      
      appliedSymmetryElements.resize(reducedSymmetry.appliedSymmetryElements.size());
      for (size_t i=0; i< reducedSymmetry.appliedSymmetryElements.size(); i++) 
	transform(reducedSymmetry.appliedSymmetryElements[i], appliedSymmetryElements[i]);

      symmetryGroup.init(appliedSymmetryElements);
      pointGroup.init(symmetryGroup.pointGroupAppliedElements());
      analyzeSpaceGroups();
    }

    void checkAppliedElements() {}

    int indexFor(int id) {
      for(size_t i=0; i<appliedSymmetryElements.size(); i++)
	if (appliedSymmetryElements[i].element.id == id)
	  return i;
      return -1;
    }
    
    AppliedSymmetryElementType& getAppliedElement(int id) {

      int i = indexFor(id);
      if (i==-1) 
	throw std::range_error("Symmetry.getAppliedElement found no element for the given id.");

      return appliedSymmetryElements[i];
    }

    //======================================================================
    
    void analyzeSpaceGroups() {

      for(size_t i=1; i<=NUMGROUPS; i++) 
	spaceGroups[i].initialize(i, (*this)); 

      spaceGroupsSet = true;
      selectSpaceGroup();
    }

    //======================================================================

    /** 
     * Return the best spaceGroup for this Crystal
     */ 
    const SpaceGroupType& spaceGroup() const {
      if (!spaceGroupsSet) {
	std::ostringstream buff;
	buff << "Symmetry::spaceGroup(): spaceGroups have not been set in:" << std::endl;
	buff << (*this) << std::endl;
	ASSERT(spaceGroupsSet,
	       std::logic_error(buff.str()));
      }

      const SpaceGroupType& result = spaceGroups[selectedSpaceGroupIndex];

      return result;
    }

    /** 
     * \brief Set the selectedSpaceGroupIndex to the space group which
     *        has a close match to this crystal and which has the
     *        largest number of symmetry operations.
     */ 
    void selectSpaceGroup() {

      ASSERT(spaceGroupsSet,
	     std::logic_error("Symmetry::selectSpaceGroup: spaceGroups have not been set."));

      static Field zero(0);
      size_t maxNumOps = 1;
      selectedSpaceGroupIndex = 1;

      for (size_t i=1; i< NUMGROUPS + 1; i++) {

	SpaceGroupType& spaceGroup = spaceGroups[i];

	if (!Algorithms::close(spaceGroup.rating, zero))
	  continue;

	size_t numOps = spaceGroup.numOperations();

	if (numOps > maxNumOps) {
	  maxNumOps = numOps;
	  selectedSpaceGroupIndex = i;
	}
      }
      spaceGroups[selectedSpaceGroupIndex].selected = true;

    }
  
  };


  //======================================================================
  
  /** \ingroup xml
   * Crystall Base xml constructor
   **/
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag opsToXML(const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry) {
    
    Tag opsTag("SymmetryOperations");

    opsTag["number"] = symmetry.appliedSymmetryElements.size();
    
    std::ostringstream buff;
    buff << "[";
    for (size_t i=0; i< symmetry.appliedSymmetryElements.size(); i++) {
      if (i!=0) buff << ",";
      buff << symmetry.appliedSymmetryElements[i].element.id;
    }
    buff << "]";
    opsTag["operationIds"] = buff.str();

    opsTag.content << std::endl << std::endl;

    for (size_t i=0; i< symmetry.appliedSymmetryElements.size(); i++) {
      opsTag.content << symmetry.appliedSymmetryElements[i].id 
		     << "  '" 
		     << symmetry.appliedSymmetryElements[i].element.name << "'" << std::endl;
    }
    return opsTag;
  }
  //======================================================================
  
  /** \ingroup xml
   * Crystall Base xml constructor
   **/
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag toXML(const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry,
	    std::string name="Symmetry") {
    
    typedef LatticeTemplate<Field,DIM,Algorithms>                         LatticeType;
    typedef SpaceGroup<Field,DIM,Occupant,Algorithms>                     SpaceGroupType;
   
    //    const LatticeType& lattice = symmetry.latpat;
    
    Tag tag(name);

    tag.add(toXML(symmetry.symmetryGroup));
    tag.add(toXML(symmetry.pointGroup,"PointGroup"));
    
//     tag.add(opsToXML(symmetry));

//     Tag allRatingsTag("AllOperations");

//     for (size_t i=0; i< symmetry.allAppliedSymmetryElements.size(); i++) {

//       allRatingsTag.add(toXML(symmetry.allAppliedSymmetryElements[i]));
//     }
    
//     tag.add(allRatingsTag);

    Tag spaceGroupsTag("SpaceGroups");
    
    spaceGroupsTag["selected"] = symmetry.selectedSpaceGroupIndex;

    if (symmetry.spaceGroupsSet) {
      for (size_t i=1; i<symmetry.spaceGroups.size(); i++) {
	const SpaceGroupType& spaceGroup = symmetry.spaceGroups[i];
	Tag sgtag = toXML(spaceGroup);
	sgtag.add(asymmetricUnitToXML(spaceGroup, symmetry.latpat));
	spaceGroupsTag.add(sgtag);
      }
    }

    tag.add(spaceGroupsTag);

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
			     const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry) {

    if (symmetry.spaceGroupsSet) {
      os << " SpaceGroup["<< symmetry.selectedSpaceGroupIndex << "]: " << symmetry.spaceGroup().nam << std::endl;
      os << "   number of ops: " << symmetry.spaceGroup().numOperations() << std::endl;
    } else {
      os << "SpaceGroups not set!" << std::endl;
    }

    return os;

  }

} /** namespace psimag */

#endif
/*@}*/
