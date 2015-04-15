//-*-C++-*-

#ifndef  Psimag_GroupMultiplicationTable
#define  Psimag_GroupMultiplicationTable

/** \ingroup symmetryConcepts */
/*@{*/

/** \file GroupMultiplicationTable.h
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
  class SymmetryGroup;

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
  class GroupMultiplicationTable:
    public Matrix<int>
  {
  public:

    typedef CartesianPosition<Field,DIM>                                           CartesianPositionType;
    typedef CellPosition<Field,DIM, Algorithms>                                    CellPositionType;

    typedef Matrix<int>                                                             BaseType;
    typedef GroupMultiplicationTable<Field,DIM,Occupant,LatticeTemplate,Algorithms> ThisType;
    typedef LatticeTemplate<Field,DIM,Algorithms>                                   LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>           LatticeWithPatternType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                                  PatternType;
    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>                  TestPatternType;

    typedef SymmetryOperation<Field,DIM,Algorithms>                                 SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>                                   SymmetryElementType;
    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms>       AppliedSymmetryElementType;
    typedef std::vector<AppliedSymmetryElementType>                                 AppliedSymmetryElementsType;
    typedef std::vector<AppliedSymmetryElementType*>                                AppliedSymmetryElementsPtrType;
    typedef std::pair<int,AppliedSymmetryElementType>                               AppliedSymmetryElementPairType;

    typedef std::map<std::string,AppliedSymmetryElementPairType>                    GroupType;
    typedef const std::pair<int,int>                                                IntPairType;
    typedef std::map<IntPairType,int>                                               MultiplicationTableMapType;

    size_t                             numberOfOperations;
    const LatticeWithPatternType&      latpat;

    GroupType                          group;
    AppliedSymmetryElementsPtrType     appliedElements; // an index into group 
   
    /** \ingroup symmetryConcepts 
     *
     * \brief The SymmetryGroup class
     */
    GroupMultiplicationTable(const LatticeWithPatternType&      lpat):
      numberOfOperations(0),
      latpat(lpat)
    {}
      
    void init(const AppliedSymmetryElementsPtrType& initialElements)
    {
      // Initial Load of the group
      for (size_t i = 0; i < initialElements.size(); i++) 
	addAppliedElement(*initialElements[i]);

      build();
    }
      
    void init(const AppliedSymmetryElementsType& initialElements)
    {
      // Initial Load of the group
      for (size_t i = 0; i < initialElements.size(); i++) 
	addAppliedElement(initialElements[i]);

      build();
    }
      
    void build()
    {
      // Allocate a map-type multiplication table
      MultiplicationTableMapType multiplicationTableMap;
      
      size_t numberAdded = 1;   // value is to ensure one iteration
      size_t lastSize    = 0;
      while(numberAdded > 0) 
	numberAdded = extendMultiplicationTable(lastSize, multiplicationTableMap);
      
      numberOfOperations = appliedElements.size();
      
      (*this).resize(numberOfOperations,numberOfOperations);
      
      // Unload the multiplicationTableMap into this table
      for (MultiplicationTableMapType::iterator itr = multiplicationTableMap.begin(); 
	   itr != multiplicationTableMap.end();
	   itr++) {
	std::pair<int,int> foo;
	
	std::pair<IntPairType,int>& tablepair = *itr;
	(*this)(tablepair.first.first,tablepair.first.second) = tablepair.second;
      }
    }
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Add the given AppliedSymmetryElement to the given pair,
     *        using appliedElements to determine the proper index and
     *        saving a reference in appliedElements.
     */
    int addAppliedElement(const AppliedSymmetryElementType&        appliedElement)
    {
      const SymmetryElementType&        element        = appliedElement.element;

      // creates entry if needed
      AppliedSymmetryElementPairType&   pair           = group[element.name]; 

      int                               newIndex       = appliedElements.size();
      pair.first      = newIndex;             // save the index of the new appliedElement
      pair.second     = appliedElement;       // copy the new applied element into the pair 
      pair.second.id  = newIndex;
      appliedElements.push_back(&pair.second);  // save a reference in appliedElements
      return newIndex;
    }
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Extend the multiplication table. Return the number of operations added.
     */
    size_t extendMultiplicationTable(size_t                                    lastSize,
				     MultiplicationTableMapType&               multiplicationTableMap)
    {
      size_t numberAdded   = 0;
      size_t numOperations = appliedElements.size();

      if (numOperations > 200) {
	//printGroup();
	throw std::range_error("extendMultiplicationTable failed to converge");
      }

      // Do the bottom extension
      for (size_t i = lastSize; i < numOperations; i++) 
	for (size_t j = 0; j < numOperations; j++) 
	  numberAdded += extendMultiplicationTable(i,j,multiplicationTableMap);
	  
      // Do the right-side extension
      for (size_t i = 0; i < lastSize; i++) 
	for (size_t j = lastSize; j < numOperations; j++) 
	  numberAdded += extendMultiplicationTable(i,j,multiplicationTableMap);

      return numberAdded;
    }
  
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Extend the multiplication table for the product cspd to
     *        i and j. Return the number of operations added, zero or 1;
     */
    size_t extendMultiplicationTable(size_t i, size_t j, 
				     MultiplicationTableMapType&  multiplicationTableMap)
    {
      //      static bool trace(false);

      // If one of the operations is the identify don't bother with the multiplication
      if (appliedElements[i]->element.type == "identity") {
	IntPairType intPair(i,j);
	multiplicationTableMap[intPair] = j;
	return 0;
      }

      if (appliedElements[j]->element.type == "identity") {
	IntPairType intPair(i,j);
	multiplicationTableMap[intPair] = i;
	return 0;
      }

      const SymmetryOperationType& oi           = appliedElements[i]->operation;
      const SymmetryOperationType& oj           = appliedElements[j]->operation;
      SymmetryOperationType        product;
      Multiply(oi,oj,product);
      
      // Note that productElement is normaized by the symmetryElement factory function
      const SymmetryElementType& productElement = symmetryElement(product, latpat);  
      const std::string&         name           = productElement.name; 
      
      typename GroupType::iterator itr = group.find(name); 
	  
      if(itr == group.end()) {   // new operation

	int newIndex = addAppliedElement(AppliedSymmetryElementType(-1, latpat, productElement));
	AppliedSymmetryElementType& newElement = *appliedElements[newIndex];

	if(! Algorithms::close(Field(0), newElement.distance)) {
	  std::ostringstream buff;
	  buff << " Product operation produced during closure did not have a zero distance!" << std::endl;
	  buff <<  appliedElements[i]->element.name << " * " << appliedElements[j]->element.name 
	       << " = " << newElement.element.unNormalizedName << " -> " << newElement.element.name <<std::endl;
	  buff << " latpat.cellParameters = " << std::endl;
	  buff << toXML(latpat.parameters) << std::endl;
	  buff << " latpat.pattern = " << std::endl;
	  buff << latpat.pattern << std::endl;
	  buff << " newElement.testPattern..pattern = " << std::endl;
	  buff << newElement.testPattern.pattern << std::endl;
	  buff << " Bad Applied element is: " << std::endl;
	  buff << toXML(newElement) << std::endl;
	  throw std::logic_error(buff.str());
	}
	return 1;
      }
      else { // existing operation
	IntPairType intPair(i,j);
	multiplicationTableMap[intPair] = (*itr).second.first;	    
	return 0;
      }
    }

  };

  template<typename Field, size_t DIM, typename Occupant,
 	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag toXML(const GroupMultiplicationTable<Field,DIM,Occupant,LatticeTemplate,Algorithms>& multiplicationTable,
	    std::string name="MultiplicationTable") {
    
    Tag result(name);

    result["rows"] = multiplicationTable.appliedElements.size();
    result["cols"] = multiplicationTable.appliedElements.size();
    std::ostringstream buff2;
    multiplicationTable.print(buff2);
    result.content << std::endl << std::endl <<buff2.str();
    return result;
  }
  
} /** namespace spimag **/

#endif

/*@}*/
