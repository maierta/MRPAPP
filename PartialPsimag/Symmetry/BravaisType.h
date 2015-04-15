//-*-C++-*-

#ifndef PSIMAG_BravaisType_H
#define PSIMAG_BravaisType_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file BravaisType.h
 *
 *  Contains the class definition for BravaisType objects and
 *  enumerations of the possible bravaisType objects by DIM
 *  specialization.
 */
 

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"
#include "SeitzVectors.h"

#include "PSIMAGAssert.h"

namespace psimag {

/** \ingroup latticeConcepts
 *
 * \brief Lattices can be classified into one of 14 BravaisTypes.
 *
 */

  template <Field, size_t DIM> class BravaisType {
    
  public:

    std::string name;
    
  };

}
#endif

/*@}*/
