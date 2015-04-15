//-*-C++-*-

#ifndef PSIMAG_FiniteGroup_H
#define PSIMAG_FiniteGroup_H

/** \ingroup symmetryConcepts */
/*@{*/

/** \file FiniteGroup.h
 *  Contains the  class definition for FiniteGroup objects.
 */
 
#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"

#include "PSIMAGAssert.h"

namespace psimag {

  typedef std::map<pair<SymmetryOperation*, SymmetryOperation*>, SymmetryOperation*> MultiplicationTable;

  /** \ingroup symmetryConcepts
   *
   *  \brief A finite group in the mathematical sense. 
   *
   * \note This is the base class.
   *
   *  For finite groups we provide a closure operation, group
   *  multiplication, etc.
   *
   *  The FinitGroup is templated on a GroupElement which must:
   * 
   *  Be comparable, support a multiplication operation.
   *  
   */
  template <GroupElement>
  class FiniteGroup {
  public:
    
    
  };
} /* namespace psimag */
#endif

/*@}*/
