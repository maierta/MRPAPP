//-*-C++-*-

#ifndef  Psimag_SpaceGroupConstructor_2D
#define  Psimag_SpaceGroupConstructor_2D

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SpaceGroupConstructor.h
 *
 */  

#include <vector>
#include <map>
#include <iostream>

#include "Crystal.h"
#include "SymmetryOperation.h"
#include "SymmetryOperations2D.h"
#include "LatticeTransformation.h"
#include "Lattice.h"
#include "Simple2DReducer.h"
#include "SpaceGroupConstructor.h"
#include "ConventionalLattice.h"
#include "ForEachCentering.h"

namespace psimag {
  
  template<typename Field, size_t DIM, typename Algorithms> class ConventionalLattice;

  template<typename Field, size_t DIM> class SpaceGroupConstructor;

  template<typename CenteringType, typename ArgType>
  class Explorer {
  public:
    static bool ForCentering(ArgType& arg) {
      return false;
    }
  };

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
   *          operation.
   */
  template<typename Field> 
  class SpaceGroupConstructor<Field,2>  {
  public: 
    
    template<typename Occupant, typename Algorithms>
    static void classify(Crystal<Field, 2, Occupant, ConventionalLattice, Algorithms>& structure) {
      

      double v = 3;
      ForEachCentering<double, 2,Explorer,double>::EXEC(v);

    }
  };

} /** namespace spimag **/
#endif

/*@}*/
