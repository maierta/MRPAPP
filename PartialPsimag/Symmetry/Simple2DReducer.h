//-*-C++-*-

#ifndef PSIMAG_Simple2DReducer_H
#define PSIMAG_Simple2DReducer_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file  Simple2DReducer.h
 *  Contains a class providing a static method for the reduction of 2D lattices.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"
#include "Lattice.h"
#include "Crystal.h"
#include "ReducedCrystal.h"

#include "PSIMAGAssert.h"

namespace psimag {

  template<typename, size_t,typename,typename> class Crystal;
  template<typename, size_t,typename,typename> class ReducedCrystal;
  template<typename, size_t,typename>          class Lattice;
  template<typename, size_t,typename>          class ReducedLattice;

  /** \ingroup latticeConcepts
   *
   * \brief This class provides a static method for the reduction of 2D cells.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   */
  template<typename Field, typename Occupant, typename Algorithms>
  class Simple2DReducer  
  {
  public:

    enum {DIM=2};

    typedef Crystal       <Field,DIM,Occupant,Algorithms>  CrystalType;
    typedef ReducedCrystal<Field,DIM,Occupant,Algorithms>  ReducedCrystalType;
    
    static void reduce(ReducedCrystalType& reduced) {

      size_t count =0;
      
      while(true) {

	if (reduced.parameters.alpha > 90)
	  reduced.makeAcute(0,1);
	
	if (!Algorithms::close(reduced.parameters.a, reduced.parameters.b) &&
	    reduced.parameters.a > reduced.parameters.b)
	  reduced.swapBasisVectors(0,1);
	
	if (!reduced.reduceDiagonal(0,1)) {

	  reduced.transform.finalize();   
	  
	  const ReducedLattice<Field,DIM,Algorithms> latticeCopyOfReduced(reduced);
	  ReducedLattice<Field,DIM,Algorithms>&      reducedLattice = reduced;
	  reduced.transform(latticeCopyOfReduced, reducedLattice); 

	  return;
	}

	if (count++ > 100) throw std::range_error("SimpleReducer2D failed to converge!");
      }
    }
    
  };
  
} /* namespace psimag */

#endif /* PSIMAG_Simple2DReducer_H */

/*@}*/
