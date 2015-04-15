//-*-C++-*-

#ifndef PSIMAG_Cell_H
#define PSIMAG_Cell_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file  PrimitiveCell.h Contains the PrimitiveCell class definition.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"

#include "PSIMAGAssert.h"

namespace psimag {

  /** \ingroup latticeConcepts
   *
   * \brief The PrimitiveCell marker class.
   *
   * A primitive cell is a Cell whose basis vectors can be used to
   * generate the translation group of the lattice.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM Dimensionality of the lattice being represented. Legal values 
   *            are 1,2,3.
   */
  template <typename Field, size_t DIM>
  class PrimitiveCell: public Cell<Field,DIM> {
    
  public:
    
    
  }
}
#endif

/*@}*/
