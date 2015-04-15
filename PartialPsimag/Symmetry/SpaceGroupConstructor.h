//-*-C++-*-

#ifndef  Psimag_SpaceGroupConstructor
#define  Psimag_SpaceGroupConstructor

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SpaceGroupConstructor.h
 *
 */  

#include <vector>
#include <map>
#include <iostream>
#include <stdexcept>

#include "Crystal.h"
#include "ConventionalCrystal.h"
#include "ConventionalLattice.h"

namespace psimag {

  // Forward declaration of CrystalStructure
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms> class Crystal;
  
  /** \ingroup symmetryConcepts 
   *
   * \brief 
   *
   * Following "Algorithms for deriving crysallographic space-group
   * information" by R.W. Grosse-Kunstleve, Section 3: Efficient
   * generation of the symmetry operations of a spage group.
   *
   * \warning The GroupElement type must provide a * operator.  This
   *          operator will be deemed to implement the group
   *          operation.  SpaceGroupConstructor
   */
  template<typename Field, size_t DIM> 
  class SpaceGroupConstructor  
  {
  public: 

    template<typename Occupant, typename Algorithms>
    static void classify(ConventionalCrystal<Field, DIM, Occupant, Algorithms>& structure) {

      std::ostringstream msg;
      msg <<  "SpaceGrpConstructor<DIM> should not have been called, call a specialization!\n";

      throw std::domain_error(msg.str());
    }
    
  };
  
} /** namespace psimag **/
#endif

/*@}*/
