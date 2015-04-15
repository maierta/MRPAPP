//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file ReducedCrystal.h  
 *
 *  Contains a the ReducedCrystal class.
 */

#ifndef Psimag_ReducedCrystal_H
#define Psimag_ReducedCrystal_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"

#include "ReducedLattice.h"
#include "CrystalBase.h"
#include "Crystal.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Occupant,typename Algorithms> class Crystal;
  template<typename Field, size_t DIM, typename Algorithms> class ReducedLattice;

  //============================================================ Reduced Crystal Structure
  /** \ingroup crystallography
   *  
   *\brief Reduced Crystal Class.
   *
   *\param DIM Dimensionality of the structure being
   *           represented. Legal values are 1,2,3.
   */ 
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  class ReducedCrystal: 
    public CrystalBase<Field, DIM,Occupant,ReducedLattice,Algorithms>  
  {  
  public:

    typedef CrystalBase<Field, DIM, Occupant,ReducedLattice,Algorithms> BaseType;
    typedef Crystal<Field,DIM,Occupant,Algorithms>                      CrystalType;
    /* 
     * \brief Construct a Reduced Crystal from a given Crystal.
     */
    ReducedCrystal(const CrystalType& crystal):
      BaseType(crystal)
    {
      Algorithms::Reducer::reduce(*this);

      this->symmetry.setAppliedSymmetryElements();
      
    }
    
  };

  //======================================================================

  /** \ingroup ostream
   * Cell output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const ReducedCrystal<Field,DIM,Occupant,Algorithms>& reducedCrystal) {

    typedef CrystalBase<Field, DIM,Occupant,ReducedLattice,Algorithms>  CrystalBaseType;

    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    
    os << "================================================== Reduced Crystal Structure" << std::endl;
    os << ( (CrystalBaseType) reducedCrystal ) << std::endl;
    os << "============================================================" << std::endl;

    
    return os;
  }

} /** namespace psimag */

#endif
/*@}*/
