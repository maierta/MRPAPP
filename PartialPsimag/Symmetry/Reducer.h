//-*-C++-*-

#ifndef PSIMAG_Reducer_H
#define PSIMAG_Reducer_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file  Reducer.h
 *
 *  Contains a prototype reducer algorithm illustrating the Reducer interface.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"
#include "ReducedCell.h"

#include "PSIMAGAssert.h"

namespace psimag {

  /** \ingroup latticeConcepts
   *
   *  \brief An example class illustrating the reducer interface.
   *
   *  A Reducer is a class which provides a static method for reducing
   *  the basis of a given lattice. Examples of classes that implement
   *  this interface are \ref LLL_Reducer and \ref GruberReducer.
   *
   *  The crystallography/lattice/reduction directory is where the
   *  various reduction codes are kept.
   *
   * <pre>
   *    template<class Reducer>
   *    ReducedCell<Field, DIM>& getReducedCell() {
   *      if (reducedCellSet) 
   *	     return reducedCell;
   *      reducedCell = Reducer::reduce(*this);
   *      return reducedCell;
   *    }
   * </pre>
   *
   * \param DIM Dimensionality of the lattice being represented. Legal values
   *            are 1,2,3.
   */
  template <typename Field, size_t DIM>
  class Reducer {

  public:

    static ReducedCell reduce(Lattice& lattice) {

      // Note that a Lattice is a cell (i.e. inherits from public
      // Cell<Field, DIM>). It is understood to be the origional cell
      // that the user defined the lattice with.
      
      ReducedCell result;

      // Compute the reduced cell . . . . .

      return result;
    }

} /* namespace psimag */

#endif /* PSIMAG_Reducer_H */

/*@}*/
