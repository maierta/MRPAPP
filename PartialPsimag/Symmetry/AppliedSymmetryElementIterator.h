//-*-C++-*-

#ifndef  Psimag_SpaceGroup2D_Iterator
#define  Psimag_SpaceGroup2D_Iterator

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SpaceGroup2D.h
 *
 *  1999, Grosse-Kunstleve, "Algorithms for deriving crystallographic space-group information"
 *
 *  some of this will be implemented in the cctbx where I may be able to reuse it.
 */  

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <bitset>
#include <iterator>
#include <stdexcept>

#include "SymmetryElements2D.h"
#include "AppliedSymmetryElement.h"
#include "SpaceGroup.h"
#include "SpaceGroup2D.h"

namespace psimag {
  
  //====================================================================== 
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroup class specialization.
   */
  template<typename Field, 
	   size_t DIM,
	   typename Occupant, 
	   template<typename, size_t, typename> class LatticeTemplate, 
	   typename Algorithms>
  class AppliedSymmetryElementIterator:
    public std::iterator< std::forward_iterator_tag, 
			  AppliedSymmetryElement<Field,DIM,Occupant,LatticeTemplate<Field,DIM,Algorithms>,Algorithms> >
  {
  public:
    typedef AppliedSymmetryElementIterator<Field,DIM,Occupant,LatticeTemplate,Algorithms> AppliedSymmetryElementIteratorType;

    typedef LatticeTemplate<Field,DIM,Algorithms>                             LatticeType;
    typedef SpaceGroup<Field,DIM,Occupant, Algorithms>                        SpaceGroupType;
    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
    typedef std::vector< AppliedSymmetryElementType >                         AppliedSymmetryElementsType;
    typedef Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>           SymmetryType;

  protected:

    bool                               atEnd;
    size_t                             position;
    const SpaceGroupType&              spaceGroup;
    const AppliedSymmetryElementsType& allElements;
    size_t                             numElements;

    bool validPosition() { 
      if (position < numElements) 
	return (spaceGroup.opFlags[position] == 1); 
      if (position == numElements)
	return true;   //end position
      return false;
    }
    
  public:


    explicit AppliedSymmetryElementIterator(const SpaceGroupType& sg, const AppliedSymmetryElementsType& allEl):
      atEnd(false),
      position(0),
      spaceGroup(sg),
      allElements(allEl),
      numElements(allEl.size())
    {
      if (!validPosition())
	(*this)++;
    }

     explicit AppliedSymmetryElementIterator(const AppliedSymmetryElementIterator& other):
      atEnd(other.atEnd),
      position(other.position),
      spaceGroup(other.spaceGroup),
      allElements(other.allElements),
      numElements(allElements.size())
    {
      if (!validPosition())
	(*this)++;
    }

   static AppliedSymmetryElementIterator end(const SpaceGroupType& sg, const AppliedSymmetryElementsType& allEl) {
      AppliedSymmetryElementIterator result(sg,allEl);
      result.position = result.numElements;
      return result;
    }

    const AppliedSymmetryElementType& operator* () {
      return allElements[position];
    }

    AppliedSymmetryElementIteratorType& operator++ () {
      position++;
      if (position > numElements) 
	throw std::out_of_range("AppliedSymmetryElementIterator iterated to far!");
      while (!validPosition())
	position++;
      return *this;
    }

    AppliedSymmetryElementIteratorType operator++ (int dummy) {
      AppliedSymmetryElementIteratorType save = *this;
      ++*this;
      return save;
    }

    bool operator== (const AppliedSymmetryElementIteratorType& other) { return (position == other.position); }
    bool operator!= (const AppliedSymmetryElementIteratorType& other) { return (position != other.position); }


  };

} /** namespace spimag **/

#endif // Psimag_SpaceGroup2D_Iterator

/*@}*/
