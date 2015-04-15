//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file ConventionalCrystal.h  
 *
 *  Contains a the ConventionalCrystal class.
 */

#ifndef PSIMAG_CONVENTIONALCRYSTAL_H
#define PSIMAG_CONVENTIONALCRYSTAL_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"

#include "ConventionalLattice.h"
#include "CrystalBase.h"
#include "ReducedCrystal.h"

namespace psimag {

  //============================================================ Conventional Crystal Structure
  /** \ingroup crystallography
   *  
   *\brief Conventional  Crystal Class.
   *
   *\param DIM Dimensionality of the structure being
   *           represented. Legal values are 1,2,3.
   */ 
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  class ConventionalCrystal:
    public CrystalBase<Field, DIM, Occupant,ConventionalLattice,Algorithms>  
  {  
  public:

    typedef CrystalBase<Field,DIM,Occupant,ConventionalLattice,Algorithms> BaseType;
    typedef ConventionalCrystal<Field,DIM,Occupant,Algorithms>             ConventionalCrystalType;
    typedef ReducedCrystal<Field,DIM,Occupant,Algorithms>                  ReducedCrystalType;
    /* 
     * \brief Construct a Conventional Crystal from a given Reduced Crystal.
     */
    ConventionalCrystal(ReducedCrystalType& structure):
      BaseType(structure)
    {
      this->transform = structure.transform; // Initialize to ReducedLattice's transform
      //Algorithms::Classifier::classify(*this);
      //updatePattern(structure.pattern, this->pattern);

    }
  };

  //======================================================================

  /** \ingroup ostream
   * Cell output stream operator 
   **/
  template<typename Field, size_t DIM, typename Occupant,typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const ConventionalCrystal<Field,DIM,Occupant,Algorithms>& conventionalCrystal) {

    typedef CrystalBase<Field, DIM,Occupant,ConventionalLattice,Algorithms>  CrystalBaseType;

    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    
    os << "================================================== Conventional Crystal Structure" << std::endl;
    os << ( (CrystalBaseType) conventionalCrystal ) << std::endl;
    os << "============================================================" << std::endl;

    
    return os;
  }


} /** namespace psimag */

#endif
/*@}*/
