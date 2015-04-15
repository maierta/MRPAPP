//-*-C++-*-

#ifndef  Psimag_Orbits
#define  Psimag_Orbits

/** \ingroup symmetryConcepts */
/*@{*/

/** \file Orbits.h
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
#include "SymmetryGroup.h"

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
  template<typename Field>
  class Orbits:
    public std::vector<std::vector<int> > 
  {
  public:
    
    // Orbits are sets of pattern elements given by position indicies.
    typedef std::set<int>                                                     OrbitSetType;
    typedef std::vector<int>                                                  OrbitType;
    
    // The Quotient of the action, also called the OrbitSpace, is the set of all of the Orbits
    typedef std::vector<OrbitType>                                            OrbitSpaceType;
    typedef OrbitSpaceType                                                    BaseType;

    //typedef LatticeTemplate<Field,DIM,Algorithms>                             LatticeType;
    //typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>           SymmetryType;
    //typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
    //typedef std::vector<AppliedSymmetryElementType*>                          AppliedSymmetryElementsType;
    //typedef SymmetryGroup<Field,DIM,Occupant,LatticeTemplate,Algorithms>      SymmetryGroupType;
    
        size_t                                   numPatternPos;
    std::map<size_t, size_t>                 posIndexToOrbit;
    
    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Built the groups orbits
     */
    Orbits()
    //:
    //appliedElements(applied),
    // groupAction(ga)
    {}

    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Built the groups orbits
     */
    Orbits(const Orbits& other)
    {
      BaseType& mat            = *this;
      const BaseType& otherMat = other;
      mat = otherMat;
      numPatternPos   = other.numPatternPos;
      posIndexToOrbit = other.posIndexToOrbit;
    }

    template<typename Occupant,size_t DIM, typename LatticeType, typename Algorithms>
    void init(const std::vector<AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> *>& appliedElements,
	      const Matrix<int>&  groupAction) {

      numPatternPos = groupAction.n_col();
      
      for (size_t p = 0; p < numPatternPos; p++) {
	
	int orbitIndex = orbitFor(p);
	if (orbitIndex == -1) 
	  orbitIndex         = makeOrbitFor(p, appliedElements, groupAction);
	posIndexToOrbit[p] = orbitIndex;
      }
    }
    
    OrbitType& orbit(size_t i) {
      return (*this)[i];
    }

    const OrbitType& orbit(size_t i) const {
      return (*this)[i];
    }

    int operator () (size_t i,size_t j) const {
      return orbit(i)[j];
    }

    /** \ingroup symmetryConcepts 
     *
     * \brief Find an index for the orbit containing the site at index
     *        p.
     *
     * \note that orbits partition the pattern
     */
    int orbitFor(size_t p) {
      for (size_t i=0; i< (*this).size();i++)
	if (std::find(orbit(i).begin(), orbit(i).end(), int(p)) != orbit(i).end())
	  return i;
      return -1;
    }

    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Make an orbit containing for the site at index
     *        posIndex.
     *
     */
    template<typename Occupant, size_t DIM, typename LatticeType, typename Algorithms>
    size_t makeOrbitFor(size_t posIndex,
			const std::vector<AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> *>& appliedElements,
			const Matrix<int>&  groupAction) {
      
      typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;

      OrbitSetType  orbitSet;
      
      size_t        numOperations = appliedElements.size();
      
      for (size_t o = 0; o < numOperations; o++) {

	int otherPosIndex = groupAction(o,posIndex);
	if (otherPosIndex >= 0) 
	  orbitSet.insert(otherPosIndex);

	else {  // complain
	  const AppliedSymmetryElementType& op = *appliedElements[o];
	  std::ostringstream                buff;
	  buff << "Error! SymmetryGroup.makeOrbitFor(" << posIndex << ") op[" << op.element.id << "," << op.element.name << "]" << std::endl;
	  buff << "    groupAction undefined:  error " << otherPosIndex << std::endl;
	  throw std::out_of_range(buff.str());
	}
      }

      size_t newIndex = this->size();
      this->resize(newIndex+1);
      OrbitType& newOrbit = orbit(newIndex);
      newOrbit.resize(orbitSet.size());
      std::copy(orbitSet.begin(), orbitSet.end(), newOrbit.begin());
      return newIndex;
    }

    std::string orbitString(size_t orbitIndex) const {
      
      std::ostringstream buff;
      buff << "[";
      for (size_t j =0; j<(*this)[orbitIndex].size();j++) {
	if (j!=0)
	  buff << ", ";
	buff << (*this)[orbitIndex][j];
      }
      buff << "]";
      return buff.str();
    }
  };

  template<typename Field>
  Tag toXML(const Orbits<Field>& orbits,
	    std::string name="Orbits") {
    
    Tag result(name);

    size_t numOrbits = orbits.size();
    result["number"] = numOrbits;
    
    for (size_t i=0; i < numOrbits; i++) {
      Tag orbitTag("orbit");
      orbitTag["orbitNumber"]  = i;
      orbitTag["positions"] = orbits.orbitString(i);
      result.add(orbitTag);
    }
      
    result.content << "pos     ->  orbit  [orbitPositions]" << std::endl;
    result.content << "------------------------------------------" << std::endl;
    for (size_t posIndex=0; posIndex < orbits.numPatternPos; posIndex++) {
	
      if (orbits.posIndexToOrbit.count(posIndex) == 0)
	result.content << posIndex << "->  -1 [] !" << std::endl;
      else {
	std::map<size_t,size_t>& indexToOrbit = const_cast<std::map<size_t,size_t>&>(orbits.posIndexToOrbit);
	const size_t orbitIndex = indexToOrbit[posIndex];
	result.content << posIndex << "->  " << orbitIndex << " " << orbits.orbitString(orbitIndex) << std::endl;
      }
    }
    return result;
  }
  
} /** namespace spimag **/

#endif

/*@}*/
