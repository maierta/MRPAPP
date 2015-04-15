//-*-C++-*-

#ifndef PSIMAG_Centering_H
#define PSIMAG_Centering_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file Centering.h
 *
 *  Contains the class definition for Centering objects and
 *  enumerations of the possible centering objects by DIM
 *  specialization.
 */
 

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"
#include "CellDirection.h"

#include "PSIMAGAssert.h"

namespace psimag {

/** \ingroup latticeConcepts
 *
 * \brief A Centering is a set of translations which are used with
 *        non-primitive cells to define a the translations of lattice.
 *
 *  The conventional centerings are (ITC page. 4):
 *
 *  DIM == 1:
 *
 *      Primitive: One lattice point per cell (i.e.) no additional
 *                 centering translations.
 *  
 *  DIM == 2:
 *
 *      Primtive: One lattice point per cell (i.e.) no additional
 *                centering translations.
 *
 *      Centered: Two lattice points per cell, additional point at (1/2, 1/2) 
 *
 *      Hexagonally Centered: Three lattice points per cell,
 *                            additional point at (2/3, 1/3), and (1/3, 2/3).
 *
 *  DIM == 3:
 *    
 *      Primtive: One lattice point per cell (i.e.) no additional
 *                centering translations.
 *
 *      A Face Centered: 
 *      B Face Centered: 
 *      C Face Centered: 
 *      Body Centered:
 *      All Face Centered: 
 *      Rhombohedrally Centered (Obverse Setting):
 *      Rhombohedrally Centered (Reverse Setting):
 *      Hexagonally Centered:
 *
 */

  template<size_t DIM>
  class Centerings {
  public:
    enum{ NumCenterings=0 };
  };

  template<>
  class Centerings<2> {
  public:
    enum{ NumCenterings=2 };
  };

  template<>
  class Centerings<3> {
  public:
    enum{ NumCenterings=6 }; // Probably more check this.
  };

  //======================================================================
  /**
   * \brief The Generic Centering2D template. Not used directly.
   */
  template<typename Field, size_t DIM>
  class CenteringBase {
  public:
  };

  //======================================================================
  /**
   * \brief The Generic Centering2D template. Not used directly.
   */
  template<typename Field, size_t DIM, size_t NUM>
  class Centering: public CenteringBase<Field,DIM> {
  public:
  };


}
#endif

/*@}*/
