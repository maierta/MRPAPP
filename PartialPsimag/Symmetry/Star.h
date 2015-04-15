//-*-C++-*-

#ifndef  Psimag_Star
#define  Psimag_Star

/** \ingroup symmetryConcepts */
/*@{*/

/** \file Star.h
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
   * \brief The Star class
   */
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class Star
  {
  public:

    typedef Star<Field,DIM,Occupant,LatticeTemplate,Algorithms>           StarType;
    typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>       SymmetryType;
    typedef LatticeTemplate<Field,DIM,Algorithms>                         LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                        PatternType;

    typedef CartesianPosition<Field,DIM>                                  CartesianPositionType;
    typedef CellPosition<Field,DIM, Algorithms>                           CellPositionType;

    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
    typedef std::pair<int,AppliedSymmetryElementType>                         AppliedSymmetryElementPairType
    typedef std::map<std::string,AppliedSymmetryElementPairType>              GroupType;
    typedef SymmetryElement<Field,DIM,Algorithms>                             SymmetryElementType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                           SymmetryOperationType;
    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>            TestPatternType;
    
    // Orbits are sets of pattern elements given by position indicies.
    typedef std::set<int>                                                 OrbitSetType;
    typedef std::vector<int>                                              OrbitType;

    // The Quotient of the action also called the OrbitSpace is the set of all of the Orbits
    typedef std::vector<OrbitType>                                        OrbitSpaceType;

    Matrix<int>              groupAction;
    Matrix<int>              multiplicationTable;

    GroupType                group;

    size_t                   numPatternPos;
    OrbitSpaceType           orbitSpace;
    std::map<size_t, size_t> posIndexToOrbit;

    
    /** \ingroup symmetryConcepts 
     *
     * \brief The Star class
     */
    Star():
      numPatternPos(0),
      orbitSpace(),
      posIndexToOrbit()
    {}
    
    /** \ingroup symmetryConcepts 
     *
     * \brief The Star class
     */
    void init(const SymmetryType& symmetry)
    {
      buildGroup(symmetry);
      buildGroupAction(symmetry);
      buildMultiplicationTable(symmetry);
      
      // build the OrbitSpace and pos to orbit map
      numPatternPos = symmetry.latpat.pattern.NUMPOS;

      for (size_t p = 0; p < numPatternPos; p++) {
	int orbitIndex = -1; // orbitFor(p);
	//if (orbitIndex == -1) 
	  orbitIndex = makeOrbitFor(p, symmetry);
	  posIndexToOrbit[p] = orbitIndex;
      }

    }
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Add the given AppliedSymmetryElement to the given pair,
     *        using appliedElements to determine the proper index and
     *        saving a reference in appliedElements.
     */
    void addElementToPair(AppliedSymmetryElementPairType&          pair,
			  const AppliedSymmetryElementType&        appliedElement,
			  std::vector<AppliedSymmetryElementType&> appliedElements)
    {
      const SymmetryElementType&        element        = appliedElement.element;
      int                               newIndex       = appliedElements.size();
      pair.first      = newIndex;             // save the index of the new appliedElement
      pair.second     = appliedElement;       // copy the new applied element into the pair 
      appliedElements.puh_back(pair.second);  // save a reference in appliedElements
    }
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Create the closed group and multiplication table.
     */
    void buildGroup(const SymmetryType& symmetry) {
      
      std::vector<AppliedSymmetryElementType&> appliedElements;
      
      size_t numOperations = symmetry.appliedSymmetryElements.size();
      
      // Initial Load of the group
      for (size_t i = 0; i < numOperations; i++) 
	addElementToPair(group[symmetry.appliedSymmetryElements[i].element.name],
			 symmetry.appliedSymmetryElements[i],
			 appliedElements);

      // Allocate a map-type multiplication table
      std::map<pair<int,int>,int> multiplicationTable;

      size_t numberAdded = 1;   // value is to ensure one iteration
      size_t lastSize    = 0;
      while(numberAdded > 0) 
	numberAdded = extendMultiplicationTable(lastSize, multiplicationTable, appliedElements,symmetry);
    }
  

    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Extend the multiplication table. Return the number of operations added.
     */
    size_t extendMultiplicationTable(size_t                                    lastSize,
				     std::map<pair<int,int>,int>&              multiplicationTable,
				     std::vector<AppliedSymmetryElementType&>& appliedElements,
				     const SymmetryType&                       symmetry)
    {
      size_t numberAdded   = 0;
      size_t numOperations = appliedElements.size();

      // Do the bottom extension
      for (size_t i = lastSize; i < numOperations; i++) 
	for (size_t j = 0; j < numOperations; j++) 
	  numberAdded += extendMultiplicationTable(i,j,multiplicationTable, appliedElements);
	  
      // Do the right-side extension
      for (size_t i = 0; i < lastSize; i++) 
	for (size_t j = lastSize; j < numOperations; j++) 
	  numberAdded += extendMultiplicationTable(i,j,multiplicationTable, appliedElements);
	  
      return numberAdded;
    }

    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Extend the multiplication table for the product cspd to
     *        i and j. Return the number of operations added, zero or 1;
     */
    size_t extendMultiplicationTable(size_t i, size_t j, 
				     std::map<pair<int,int>,int>&              multiplicationTable,
				     std::vector<AppliedSymmetryElementType&>& appliedElements,
				     const SymmetryType&                       symmetry)
    {
      const SymmetryOperationType& oi           = appliedElements[i].operation;
      const SymmetryOperationType& oj           = appliedElements[j].operation;
      SymmetryOperationType        product;
      Multiply(oi,oj,product);
      
      SymmetryElementType& productElement = symmetryElement(product, symmetry.latpat);
      std::string&         name           = productElement.name; 
      
      GroupType::Iterator itr = group.find(name); 
	  
      if(itr == group.end()) {   // new operation
	AppliedSymmetryElementPairType&   pair              = group[name];  // create new group entry
	AppliedSymmetryElementType&       newAppliedElement = pair.second;  // take a reference to its AppliedSymmetryElement
	int                               newIndex          = appliedElements.size();
	newAppliedElement = AppliedSymmetryElement(newIndex, symmetry.latpat, productElement);
	pair.first = newIndex;
	appliedElements.push_back(newAppliedElement);
	return 1;
      }
      else { // existing operation
	AppliedSymmetryElementPairType& pair = itr*; 
	multiplicationTable[std::pair<int,int>(i,j)] = pair.first;	    
	return 0;
      }
    }
    //============================================================ 
    // 

    void buildMultiplicationTable(const SymmetryType& symmetry) {

      size_t numOperations = symmetry.appliedSymmetryElements.size();

      multiplicationTable.resize(numOperations, numOperations);

      for (size_t i = 0; i < numOperations; i++) {
	for (size_t j = 0; j < numOperations; j++) {
 
	  const SymmetryOperationType& oi           = symmetry.appliedSymmetryElements[i].operation;
	  const SymmetryOperationType& oj           = symmetry.appliedSymmetryElements[j].operation;

	  SymmetryOperationType  product;
	  Multiply(oi,oj,product);

	  const SymmetryElementType& productElement = symmetryElement(product,symmetry.latpat);
	  //product.normalize();

	  int productIndex = findIndex(symmetry, productElement);

	  if ( productIndex == -1) {
	    std::ostringstream msg;
	    msg << "\nCould not find a match for product:" << "\n";
	    msg << toXML(product) << "\n";
	    msg << "producedFrom:" <<  symmetry.appliedSymmetryElements[i].element.name 
		      << " times " << symmetry.appliedSymmetryElements[j].element.name 
		      << " => " << productElement.name << "\n";
	    throw std::logic_error(msg.str());
	  }
      	
	  multiplicationTable(i,j) = productIndex;
	}
      }      
    }

    //============================================================ 
    // 
    int findIndex(const SymmetryType&          symmetry, 
		  const SymmetryElementType&   productElement) {

      std::vector<AppliedSymmetryElementType> elements = symmetry.allAppliedSymmetryElements;  // For testing
      //std::vector<AppliedSymmetryElementType> elements = symmetry.appliedSymmetryElements;  // This is the real one

      size_t numOperations = elements.size();

      for (size_t i = 0; i < numOperations; i++) {
	const SymmetryElementType& candidateEl   = elements[i].element;
	if (candidateEl.name == productElement.name)
	  return i;
      }
      
      return -1;
    }

    //============================================================ 
    //
    // http://en.wikipedia.org/wiki/Group_action

    void buildGroupAction(const SymmetryType& symmetry) {

      size_t                        numOperations = symmetry.appliedSymmetryElements.size();
      const LatticeWithPatternType& latpat        = symmetry.latpat;
      const PatternType&            pattern       = latpat.pattern;
      
      groupAction.resize(numOperations, pattern.NUMPOS);

      for (size_t o = 0; o < numOperations; o++) 
 	for (size_t p = 0; p < pattern.NUMPOS; p++) 
	  groupAction(o,p) = -1;    // This means that the index is undefined, i.e. it would map an occupant into a different occupant

      for (size_t o = 0; o < numOperations; o++) {
	
 	const AppliedSymmetryElementType& op           = symmetry.appliedSymmetryElements[o];
	const TestPatternType&            testPattern  = op.testPattern;

	for (size_t p = 0; p < pattern.NUMPOS; p++) {
	  try {
	    groupAction(o,p) = testPattern.pattern.indexOf(pattern.cellPositions[p]);
	  }
	  catch (std::out_of_range&) {	  
	    groupAction(o,p) = -2;
	  }
	}
      }
    }       


    //============================================================ 
    //
    // Note that orbits partition the pattern
    //
    int orbitFor(size_t p) {
      for (size_t i=0; i< orbitSpace.size();i++)
	if (std::find(orbitSpace[i].begin(), orbitSpace[i].end(), int(p)) != orbitSpace[i].end())
	  return i;
      return -1;
    }

    size_t makeOrbitFor(size_t posIndex, const SymmetryType& symmetry) {

      OrbitSetType  orbitSet;
      
      size_t             numOperations = symmetry.appliedSymmetryElements.size();

      for (size_t o = 0; o < numOperations; o++) {

	int otherPosIndex = groupAction(o,posIndex);
	if (otherPosIndex >= 0) 
	  orbitSet.insert(otherPosIndex);
	else {  // complain
	  const AppliedSymmetryElementType& op = symmetry.appliedSymmetryElements[o];
	  std::ostringstream                buff;
	  buff << "Error! Star.makeOrbitFor(" << posIndex << ") op[" << op.element.id << "," << op.element.name << "]" << std::endl;
	  buff << "    groupAction undefined:  error " << otherPosIndex << std::endl;
	  throw std::out_of_range(buff.str());
	}
      }

      size_t newIndex = orbitSpace.size();
      orbitSpace.resize(newIndex+1);
      OrbitType& newOrbit = orbitSpace[newIndex];
      newOrbit.resize(orbitSet.size());
      std::copy(orbitSet.begin(), orbitSet.end(), newOrbit.begin());
      return newIndex;
    }

    const Occupant& getOccupant(size_t posIndex, const SymmetryType& symmetry) const {
      return symmetry.latpat.pattern.getOccupant(posIndex);
    }

    const CartesianPositionType& getCartesianPosition(size_t posIndex, const SymmetryType& symmetry) const {
      return symmetry.latpat.pattern.cartesianPositions[posIndex];
    }

    const CellPositionType&   getCellPosition(size_t posIndex, const SymmetryType& symmetry) const {
      return symmetry.latpat.pattern.cellPositions[posIndex];
    }

    std::string orbitString(size_t orbitIndex) const {

      const OrbitType& orbit = orbitSpace[orbitIndex];

      std::ostringstream buff;
      buff << "[";
      for (size_t j =0; j<orbit.size();j++) {
	if (j!=0)
	  buff << ", ";
	buff << orbit[j];
      }
      buff << "]";
      return buff.str();
    
    }
    
  };

  template<typename Field, size_t DIM, typename Occupant,
 	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag toXML(const Star<Field,DIM,Occupant,LatticeTemplate,Algorithms>& star,
	    std::string name="Star") {
    
    Tag result(name);

    Tag groupActionTag("GroupAction");
    std::ostringstream buff;
    star.groupAction.print(buff);
    groupActionTag.content << std::endl << std::endl <<buff.str();

    result.add(groupActionTag);

    Tag multiplicationTableTag("MultiplicationTable");
    std::ostringstream buff2;
    star.multiplicationTable.print(buff2);
    multiplicationTableTag.content << std::endl << std::endl <<buff2.str();

    result.add(multiplicationTableTag);

    Tag orbitSpaceTag("OrbitSpace");
    size_t numOrbits = star.orbitSpace.size();
    orbitSpaceTag["numOrbits"] = numOrbits;

    for (size_t i=0; i < numOrbits; i++) {
      Tag orbitTag("orbit");
      orbitTag["orbitNumber"]  = i;
      orbitTag["positions"] = star.orbitString(i);
      orbitSpaceTag.add(orbitTag);
    }

    orbitSpaceTag.content << "pos     ->  orbit  [orbitPositions]" << std::endl;
    orbitSpaceTag.content << "------------------------------------------" << std::endl;
    for (size_t posIndex=0; posIndex < star.numPatternPos; posIndex++) {

      if (star.posIndexToOrbit.count(posIndex) == 0)
	orbitSpaceTag.content << posIndex << "->  -1 [] !" << std::endl;
      else {
	std::map<size_t,size_t>& indexToOrbit = const_cast<std::map<size_t,size_t>&>(star.posIndexToOrbit);
	const size_t orbitIndex = indexToOrbit[posIndex];
	orbitSpaceTag.content << posIndex << "->  " << orbitIndex << " " << star.orbitString(orbitIndex) << std::endl;
      }
    }
    result.add(orbitSpaceTag);
    return result;
  }
  
} /** namespace spimag **/

#include "Symmetry.h"

#endif

/*@}*/
