//-*-C++-*-

#ifndef PSIMAG_Reciprocal_H
#define PSIMAG_Reciprocal_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file Reciprocal.h Contains the template class definitions for
 *        producing Reciprocals of Lattices and structures.
 **/

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Mat.h"
#include "Real.h"

#include "Lattice.h"
#include "PSIMAGAssert.h"
#include "Pattern.h"
#include "Crystal.h"
#include "ReducedCrystal.h"
#include "CartesianTranslation.h"

namespace psimag {


  template<typename Field, size_t DIM, typename Algorithms> class Lattice;

  //================================================== 1D case
  //
  /** \ingroup latticeConcepts
   *
   * \brief Function for computing 1D Reciprocal Lattices
   *
   */
  template <typename Field, typename Algorithms> 
  void Reciprocal(const Lattice<Field,1,Algorithms>& givenLattice, 
		  Lattice<Field,1,Algorithms>&       result) {
	
    typedef Lattice<Field,1,Algorithms>   LatticeType;
    typedef CartesianTranslation<Field,1> CartesianTranslationType;

    Field vol = givenLattice->volume();
    CartesianTranslationType b0 = (*givenLattice)[0];
    std::vector<CartesianTranslationType> newBasis(1);
	
    newBasis[0] = b0 / vol;
	
    LatticeType tmp(newBasis);
    result = tmp;
  }

  //================================================== 2D case

  /** \ingroup latticeConcepts
   *
   * \brief Function for computing 2D Reciprocal Lattices
   *
   */
  template <typename Field, 
	    template<typename,size_t,typename> class LatticeTemplate1,
	    template<typename,size_t,typename> class LatticeTemplate2,
	    typename Algorithms> 
  LatticeTemplate2<Field,2,Algorithms> reciprocate(const LatticeTemplate1<Field,2,Algorithms>& givenLattice) {
    
    typedef LatticeTemplate2<Field,2,Algorithms>    LatticeType;
    typedef CartesianTranslation<Field,3>           CartesianTranslationType;
    
    Field vol = givenLattice.volume();
    
    CartesianTranslationType b0 = cartesianTranslation<Field>( givenLattice(0,0), givenLattice(1,0), Field(0));
    CartesianTranslationType b1 = cartesianTranslation<Field>( givenLattice(0,1), givenLattice(1,1), Field(0));
    CartesianTranslationType b2 = cartesianTranslation<Field>( Field(0), Field(0), Field(1));
    
    CartesianTranslationType newB0 = ((b1 % b2) / vol) * Field(2) * convert<Field>(M_PI);
    CartesianTranslationType newB1 = ((b0 % b2) / vol) * Field(2) * convert<Field>(M_PI);

    LatticeType result;

    result(0,0) = newB0[0];
    result(1,0) = newB0[1];
    result(0,1) = newB1[0];
    result(1,1) = newB1[1];

    result.update();

    return result;
  }


  /** \ingroup latticeConcepts
   *
   * \brief Function for computing 2D Reciprocal Lattices
   *
   */
  template <typename Field, 
	    template<typename,size_t,typename> class LatticeTemplate1,
	    template<typename,size_t,typename> class LatticeTemplate2,
	    typename Algorithms> 
  void Reciprocal(const LatticeTemplate1<Field,2,Algorithms>& givenLattice, 
		  LatticeTemplate2<Field,2,Algorithms>&       result) {
	
    typedef LatticeTemplate2<Field,2,Algorithms>   LatticeType;
    typedef CartesianTranslation<Field,3>          CartesianTranslation3Type;
    typedef CartesianTranslation<Field,2>          CartesianTranslation2Type;

    Field vol = givenLattice->volume();
      
    CartesianTranslation3Type b0 = (*givenLattice)[0];
    CartesianTranslation3Type b1 = (*givenLattice)[1];
    CartesianTranslation3Type b2(0,0,1);

    std::vector<CartesianTranslation2Type> newBasis(2);

    newBasis[0] = (b1 % b2) / vol;
    newBasis[1] = (b0 % b2) / vol;

    LatticeType tmp(newBasis);
    result = tmp;
  }

  /** \ingroup latticeConcepts
   *
   * \brief Function for computing 2D Reciprocal Lattices
   *
   */
  template <typename Field, 
	    template<typename,size_t,typename> class LatticeTemplate,
	    typename Algorithms> 
  LatticeTemplate<Field,2,Algorithms> Reciprocal(const LatticeTemplate<Field,2,Algorithms>& givenLattice)
  { 
	
    typedef LatticeTemplate<Field,2,Algorithms>   LatticeType;
    typedef CartesianTranslation<Field,3>         CartesianTranslationType;

    Field vol = givenLattice->volume();
      
    CartesianTranslationType b0 = (*givenLattice)[0];
    CartesianTranslationType b1 = (*givenLattice)[1];
    CartesianTranslationType b2 = cartesianTranslation<Field>("0","0","1");

    std::vector<CartesianTranslationType> newBasis(2);

    newBasis[0] = (b1 % b2) / vol;
    newBasis[1] = (b0 % b2) / vol;

    return LatticeType(newBasis);
  }

  //================================================== 3D case
  //
  /** \ingroup latticeConcepts
   *
   * \brief Function for computing 2D Reciprocal Lattices
   *
   */
  template <typename Field, typename Algorithms> 
  void Reciprocal(const Lattice<Field,3,Algorithms>& givenLattice, 
		  Lattice<Field,3,Algorithms>&       result) {
	
    typedef Lattice<Field,3,Algorithms>   LatticeType;
    typedef CartesianTranslation<Field,3> CartesianTranslationType;

    Field vol = givenLattice->volume();
    
    CartesianTranslationType b0 = (*givenLattice)[0];
    CartesianTranslationType b1 = (*givenLattice)[1];
    CartesianTranslationType b2 = (*givenLattice)[2];
    
    std::vector<CartesianTranslationType> newBasis(3);
    
    newBasis[0] = (b1 % b2) / vol;
    newBasis[1] = (b0 % b2) / vol;
    newBasis[2] = (b0 % b1) / vol;
    
    LatticeType tmp(newBasis);
    result = tmp;
  }


} // namespace psimag


#endif //PSIMAG_Reciprocal_H


/*@}*/

