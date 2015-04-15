//-*-C++-*-

#ifndef PSIMAG_ReciprocalLattice_H
#define PSIMAG_ReciprocalLattice_H

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
#include "CartesianTranslation.h"

namespace psimag {


  template<typename Field, size_t DIM, typename Algorithms> class Lattice;

  //================================================== generic case
  template <typename Field, 
	    size_t   DIM, 
	    template<typename, 
		     size_t  , 
		     typename > class LatticeTemplate,
	    typename Algorithms> 
  class ReciprocalLattice {};

  //========================================================================================== 1D case
  //========================================================================================== 

  template <typename Field, 
	    template<typename , 
		     size_t   , 
		     typename > class LatticeTemplate, 
	    typename Algorithms> 
  class ReciprocalLattice<Field,1,LatticeTemplate,Algorithms>
  {
    //====================================================================== typedefs
    enum {DIM=1};
    typedef LatticeTemplate<Field,DIM,Algorithms>  LatticeType;
    typedef CartesianTranslation<Field,DIM>        CartesianTranslationType;
    typedef std::vector<CartesianTranslationType>  BasisVectorsType;
   //====================================================================== Member
    const LatticeType& lattice;
    //====================================================================== Constructor
    ReciprocalLattice(LatticeType& lat): lattice(lat) {}
    //====================================================================== Conversion Operator
//     template< template<typename, size_t, typename> class OtherLatticeTemplate >
//     operator OtherLatticeTemplate<Field,DIM,Algorithms>() {

//       typedef OtherLatticeTemplate<Field,DIM,Algorithms> OtherLatticeType;
     
//       Field vol = lattice.volume();
//       CartesianTranslationType b0 = lattice[0];
//       BasisVectorsType newBasis(1);
      
//       newBasis[0] = b0 / vol;
      
//       return OtherLatticeType(newBasis);
//     }
  };

  //========================================================================================== 2D case
  //========================================================================================== 

  template <typename Field, 
	    template<typename , size_t, typename > class LatticeTemplate, 
	    typename Algorithms> 
  class ReciprocalLattice<Field,2,LatticeTemplate,Algorithms>
  {
  public:
    //====================================================================== typedefs
    enum {DIM=2};
    typedef LatticeTemplate<Field,DIM,Algorithms>  LatticeType;
    typedef CartesianTranslation<Field,DIM+1>      CartesianTranslationType;
    typedef std::vector<CartesianTranslation<Field,DIM> >  BasisVectorsType;
    //====================================================================== Member
    const LatticeType& lattice;
    //====================================================================== Constructor
    ReciprocalLattice(const LatticeType& lat): lattice(lat) {}
    //====================================================================== Converter

    LatticeType getLattice() {
    
      Field vol = lattice.volume();
      
      CartesianTranslationType b0 = cartesianTranslation<Field>( lattice(0,0), lattice(1,0), Field(0) );
      CartesianTranslationType b1 = cartesianTranslation<Field>( lattice(0,1), lattice(1,1), Field(0) );
      CartesianTranslationType b2 = cartesianTranslation<Field>( Field(0),     Field(0),     Field(1) );
      
      CartesianTranslationType newB0 = ((b1 % b2) / vol) * Field(2) * convert<Field>(M_PI);
      CartesianTranslationType newB1 = ((b0 % b2) / vol) * Field(2) * convert<Field>(M_PI);

      BasisVectorsType basis(DIM);

      basis[0][0] = newB0[0];
      basis[1][0] = newB0[1];
      basis[0][1] = newB1[0];
      basis[1][1] = newB1[1];

      return LatticeType(basis);
    }

    //====================================================================== Conversion Operator

    operator LatticeType() {
      
      Field vol = lattice.volume();
      
      CartesianTranslationType b0 = cartesianTranslation<Field>( lattice(0,0), lattice(1,0), Field(0) );
      CartesianTranslationType b1 = cartesianTranslation<Field>( lattice(0,1), lattice(1,1), Field(0) );
      CartesianTranslationType b2 = cartesianTranslation<Field>( Field(0), Field(0), Field(1) );
      
      CartesianTranslationType newB0 = ((b1 % b2) / vol) * Field(2) * convert<Field>(M_PI);
      CartesianTranslationType newB1 = ((b0 % b2) / vol) * Field(2) * convert<Field>(M_PI);

      BasisVectorsType basis;

      basis(0,0) = newB0[0];
      basis(1,0) = newB0[1];
      basis(0,1) = newB1[0];
      basis(1,1) = newB1[1];
      
      return LatticeType(basis);
    }
  };

  //========================================================================================== 3D case
  //========================================================================================== 

  template <typename Field, 
	    template<typename, size_t, typename> class LatticeTemplate, 
	    typename Algorithms> 
  class ReciprocalLattice<Field,3,LatticeTemplate,Algorithms>
  {
  public:
    //====================================================================== typedefs
    enum {DIM=3};
    typedef LatticeTemplate<Field,DIM,Algorithms>  LatticeType;
    typedef CartesianTranslation<Field,DIM>        CartesianTranslationType;
    typedef std::vector<CartesianTranslationType>  BasisVectorsType;
    //====================================================================== Member
    const LatticeType& lattice;
    //====================================================================== Constructor
    ReciprocalLattice(LatticeType& lat): lattice(lat) {}
    //====================================================================== Conversion Operator
 //    template< template<typename, 
// 		       size_t, 
// 		       typename> class OtherLatticeTemplate >
//     operator OtherLatticeTemplate<Field,DIM,Algorithms>() {
      
//       typedef OtherLatticeTemplate<Field,DIM,Algorithms> OtherLatticeType;
//       Field vol = lattice.volume();

//       CartesianTranslationType b0 = lattice[0];
//       CartesianTranslationType b1 = lattice[1];
//       CartesianTranslationType b2 = lattice[2];
      
//       BasisVectorsType basis;

//       basis[0] = ((b1 % b2) / vol) * Field(2) * convert<Field>(M_PI);
//       basis[1] = ((b0 % b2) / vol) * Field(2) * convert<Field>(M_PI);
//       basis[2] = ((b0 % b1) / vol) * Field(2) * convert<Field>(M_PI);

//       return OtherLatticeType(basis);
//     }
  };
    
} // namespace psimag


#endif //PSIMAG_Reciprocal_H


/*@}*/

