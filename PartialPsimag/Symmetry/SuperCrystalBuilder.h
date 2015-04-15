//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file SuperCrystalBuilder.h  
 *
 *  Contains the SuperCrystalBuiler classes
 */

#ifndef Psimag_SuperCrystalBuilder_H
#define Psimag_SuperCrystalBuilder_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"
#include "Matrix.h"

#include "CellPosition.h"
#include "LatticeCoordinates.h"
#include "Pattern.h"
#include "CrystalBase.h"


namespace psimag {


  //======================================================================

  /** 
   * \brief This is a Helper Class for SuperCrystal construction.
   *
   * Objects of this class are constructed and passed to a
   * SuperCrystal Constructor which uses the getSubLatticeVecCoords()
   * and getReciprocalSuperCrystalPattern() methods as it initializes
   * its components.
   *
   * Only specializations of this template are ever used.
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms>
  class SuperCrystalBuilder {};

  //======================================================================
  /** 
   * \brief This is the 2D version of SuperCrystalBuilder, a Helper
   *        class for SuperCrystal. 
   *
   *  It is constructed given an 2D super-lattice specification which
   *  is given by four integers. It assembles the the sub-lattice
   *  vector coordinates from the given integers, builds a generic
   *  occupant called k-point and a simple pattern from the
   *  k-point. (k-point at the origin of the cell).
   */
  template<typename Field, typename Occupant, typename Algorithms>
  class SuperCrystalBuilder<Field,2,Occupant,Algorithms> {
  public:
    enum {DIM=2};
    typedef LatticeCoordinates<DIM>                    SubLatticeVecCoordType;
    typedef std::vector<SubLatticeVecCoordType>        SubLatticeVecCoordsType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>     PatternType;
    typedef PatternData<Field,DIM,Occupant,Algorithms> PatternDataType;

    SubLatticeVecCoordsType subLatticeVecCoords;
    Occupant                occupant;
    PatternType             reciprocalSuperCrystalPattern; 

    SuperCrystalBuilder(int n0, int m0, int n1, int m1):
      subLatticeVecCoords(2),
      occupant("kpoint","red"),
      reciprocalSuperCrystalPattern(occupant)
    {
      subLatticeVecCoords[0][0] = n0;
      subLatticeVecCoords[0][1] = m0;
      subLatticeVecCoords[1][0] = n1;
      subLatticeVecCoords[1][1] = m1;
    }

    SuperCrystalBuilder(int n0, int m0, int n1, int m1, PatternDataType patternData):
      subLatticeVecCoords(2),
      reciprocalSuperCrystalPattern(patternData)
    {
      subLatticeVecCoords[0][0] = n0;
      subLatticeVecCoords[0][1] = m0;
      subLatticeVecCoords[1][0] = n1;
      subLatticeVecCoords[1][1] = m1;
    }

    const PatternType getReciprocalSuperCrystalPattern() const {

      return reciprocalSuperCrystalPattern;
    }

    const SubLatticeVecCoordsType& getSubLatticeVecCoords() const {
      return subLatticeVecCoords;
    }
  };

  /** \ingroup ostream
   *
   * Output stream operator for SuperCrystal
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms>
  inline  std::ostream& operator << (std::ostream& os, 
				     const SuperCrystalBuilder<Field,DIM,Occupant,Algorithms>& builder) {

    typedef CrystalBase<Field,DIM,Occupant,Lattice,Algorithms>            BaseType;
    
    os << "====================================================================== SuperCrystalBuilder" << std::endl;
    os << "subLatticeVecCoords:" << std::endl;
    for (size_t i=0; i< DIM; i++)
      os << "coord[" << i << "] = " << builder.subLatticeVecCoords[i] << std::endl;
    os << "occupant:" << std::endl;
    os << builder.occupant << std::endl;
    os << "reciprocalSuperCrystalPattern:" <<std::endl;
    os << builder.reciprocalSuperCrystalPattern << std::endl;
    return os;
  }

} /** namespace psimag */

#endif
/*@}*/
