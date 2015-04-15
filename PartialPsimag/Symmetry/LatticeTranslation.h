//-*-C++-*-
#ifndef  Psimag_LatticeTranslation
#define  Psimag_LatticeTranslation

/** \ingroup extendedMatrices */
/*@{*/

/** \file LatticeTranslation.h
 *
 * Contains a class for implementing LatticeTranslation (SeitzTranslation subclass) objects.
 */

#include "SeitzVector.h"
#include "Lattice.h"

template <typename Field,size_t DIM,typename Algorithms> class Lattice;

namespace psimag {
  
  /** \ingroup extendedMatrices
   *
   * \brief A marker class indicating that the Translation indicates a direction from within a Cell.
   *
   * These can be converted to  and from MillerIndices for a given Cell.
   */
  template<size_t DIM> 
  class LatticeTranslation: public SeitzVector< int, DIM, 0 > {
    
  public:

    /** The Default Constructor produces the zero translation. */
    LatticeTranslation(): SeitzVector<int, DIM, 0>() {}

    /** Construct a translation whose components are set to the given value. */
    LatticeTranslation(int val): SeitzVector<int, DIM, 0>(val) {}

    /** Copy Construct a translation. */
    LatticeTranslation(const LatticeTranslation<DIM>& v): SeitzVector<int, DIM, 0>(v) {}

    template<typename Field, typename Algorithms>
    CartesianTranslation<Field,DIM> cartesianTranslation(Lattice<Field,DIM,Algorithms> lat) {
      CartesianTranslation<Field,DIM> result;
      for (size_t i=0; i<DIM; i++) 
	result += lat[i] * (*this)[i];
      return result;
    }

  };

  LatticeTranslation<2> latticeTranslation(int t0, int t1) {
    LatticeTranslation<2> result;
    result[0] = t0;
    result[1] = t1;
    return result;
  }

  LatticeTranslation<2> latticeTranslation(int t0, int t1, int t2) {
    LatticeTranslation<2> result;
    result[0] = t0;
    result[1] = t1;
    result[2] = t2;
    return result;
  }

  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<size_t DIM> 
  Tag toXML(const LatticeTranslation<DIM> lTrans,
	    std::string name="LatticeTranslation") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << lTrans[i] ;
    return result;
  }


}  /* namspace psimag */

#endif 
/*@}*/
