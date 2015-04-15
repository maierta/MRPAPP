//-*-C++-*-

#ifndef  Psimag_SymmetryGroup
#define  Psimag_SymmetryGroup

/** \ingroup symmetryConcepts */
/*@{*/

/** \file Group.h
 *
 */  

#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>

#include "SymmetryOperation.h"
#include "AppliedSymmetryElement.h"
#include "Pattern.h"
#include "LatticeWithPattern.h"
#include "TestPattern.h"
#include "CellPosition.h"
#include "CartesianPosition.h"
#include "SymmetryElement.h"
#include "AppliedSymmetryElement.h"
#include "Matrix.h"
#include "Orbits.h"
#include "GroupAction.h"
#include "GroupMultiplicationTable.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class Symmetry;

  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class AppliedSymmetryElement;

  /** \ingroup symmetryConcepts 
   *
   * \brief The Group class
   *
   * \ note http://en.wikipedia.org/wiki/Group_action
   */
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class SymmetryGroup
  {
  public:

    typedef LatticeTemplate<Field,DIM,Algorithms>                             LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>     LatticeWithPatternType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                            PatternType;
    typedef CartesianPosition<Field,DIM>                                      CartesianPositionType;
    typedef CellPosition<Field,DIM, Algorithms>                               CellPositionType;

    typedef SymmetryGroup<Field,DIM,Occupant,LatticeTemplate,Algorithms>      ThisType;
    typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>           SymmetryType;
    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
    typedef std::vector<AppliedSymmetryElementType>                           AppliedSymmetryElementsType;
    typedef std::vector<AppliedSymmetryElementType*>                          AppliedSymmetryElementsPtrType;

    typedef std::pair<int,AppliedSymmetryElementType>                         AppliedSymmetryElementPairType;
    typedef std::map<std::string,AppliedSymmetryElementPairType>              GroupType;
    typedef SymmetryElement<Field,DIM,Algorithms>                             SymmetryElementType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                           SymmetryOperationType;
    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>            TestPatternType;
    typedef const std::pair<int,int>                                          IntPairType;
    typedef std::map<IntPairType,int>                                         MultiplicationTableMapType;
    typedef std::vector<int>                                                  PointPermutationType;

    typedef Orbits<Field>                                                     OrbitsType;
    typedef GroupAction<Field,DIM,Occupant,LatticeTemplate,Algorithms>        GroupActionType;

    typedef GroupMultiplicationTable<Field,DIM,Occupant,LatticeTemplate,Algorithms> GroupMultiplicationTableType;

    size_t                                   numberOfOperations;
    size_t                                   numPatternPos;

    const LatticeWithPatternType&            latpat;

    GroupMultiplicationTableType             multiplicationTable;
    const GroupType&                         group;           // kept in multiplication table
    const AppliedSymmetryElementsPtrType&    appliedElements; // kept in multiplication table, an index into group
    GroupActionType                          groupAction;
    OrbitsType                               orbits;


    /** \ingroup symmetryConcepts 
     *
     * \brief The SymmetryGroup class
     *
     * \note We cannot use latpat right away since it may be reduced
     *       after construction?  This also means that the
     *       initialization must be postponed.
     */
    SymmetryGroup(const LatticeWithPatternType&  lpat):
      numberOfOperations(0),
      numPatternPos(0),
      latpat(lpat),
      multiplicationTable(latpat),
      group(multiplicationTable.group),
      appliedElements(multiplicationTable.appliedElements),
      groupAction(appliedElements,latpat.pattern),
      orbits()
    {}
    
    /** \ingroup symmetryConcepts 
     *
     * \brief Initialize the SymmetryGroup with an itinial set of Applied Elements
     */
    void init(const AppliedSymmetryElementsType& initalAppliedElements)
    {
      multiplicationTable.init(initalAppliedElements);
      groupAction.init();
      orbits.init(appliedElements,groupAction);
    }

    /** \ingroup symmetryConcepts 
     *
     * \brief Initialize the SymmetryGroup with an itinial set of Applied Element Ptrs
     */
    void init(const AppliedSymmetryElementsPtrType& initalAppliedElements)
    {
      multiplicationTable.init(initalAppliedElements);
      groupAction.init();
      orbits.init(appliedElements,groupAction);
    }

    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief 
     */
    void printGroup(std::ostream& os) {

      for (size_t i =0; i < appliedElements.size(); i++) {
	AppliedSymmetryElementType& appliedElement = *appliedElements[i];
	os << appliedElement.element.name << std::endl;
      }
    }
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief 
     */
    AppliedSymmetryElementsPtrType pointGroupAppliedElements() const {
      
      AppliedSymmetryElementsPtrType result;
      CellPositionType               zeroPos;
      int numOps = appliedElements.size();
      
      for (int i =0; i < numOps; i++) {
	AppliedSymmetryElementType& appliedElement = *appliedElements[i];
	if (appliedElement.element.type == "translation") continue;
	if (appliedElement.element.type == "glide") continue;
	if (appliedElement.element.cellPosition.closeTo(zeroPos)) {
	  result.push_back(appliedElements[i]);
	}
      }
      return result;
    }
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief 
     */
    Tag makeAppliedSymmetryElementsTag(const std::string& name="AppliedSymmetryElements") const {
      Tag result(name);
      for (size_t i=0; i < appliedElements.size(); i++) 
	result.add(toXML(*appliedElements[i]));
      return result;
    }

  };

  //============================================================ 
  /** \ingroup xml
   *
   * SymmetryGroup XML Tag factory
   *
   **/
  template<typename Field, size_t DIM, typename Occupant,
 	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag toXML(const SymmetryGroup<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetryGroup,
	    std::string name="SymmetryGroup") {
    
    Tag result(name);

    result.add(toXML(symmetryGroup.multiplicationTable));
    result.add(toXML(symmetryGroup.groupAction));
    result.add(toXML(symmetryGroup.orbits));

    //std::vector<int> unique = symmetryGroup.uniquePermutationOps();
    //result.add(symmetryGroup.makeGroupActionTag("UniqueGroupAction",unique));

    result.add(symmetryGroup.makeAppliedSymmetryElementsTag());
    //result.add(symmetryGroup.makePointGroupTag());

    return result;
  }
  
} /** namespace spimag **/

#include "Symmetry.h"

#endif

/*@}*/
