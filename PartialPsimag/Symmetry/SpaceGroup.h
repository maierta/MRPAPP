//-*-C++-*-

#ifndef  Psimag_SpaceGroup
#define  Psimag_SpaceGroup

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SpaceGroup.h
 *
 *  1999, Grosse-Kunstleve, "Algorithms for deriving crystallographic space-group information"
 *
 *  some of this will be implemented in the cctbx where I may be able to reuse it.
 */  

#include <vector>
#include <map>
#include <iostream>

#include "SymmetryOperation.h"

namespace psimag {

  /** \ingroup symmetryConcepts 
   *
   * \brief The SpaceGroup class
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms> 
  class SpaceGroup
  {};

  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms> 
  std::ostream& operator << (std::ostream& os, SpaceGroup<Field,DIM,Occupant,Algorithms> spaceGroup) {
    
    os << "Generic SpaceGroup: This should not be used, only specializations of this class." << std::endl;
    return os;
  }
  
} /** namespace spimag **/
#endif

/*@}*/
