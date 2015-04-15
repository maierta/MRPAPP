//-*-C++-*-

#ifndef PSIMAG_ReducedLattice_H
#define PSIMAG_ReducedLattice_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file  ReducedLattice.h Contains the ReducedLattice class definition.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"

#include "PSIMAGAssert.h"
#include "Lattice.h"
#include "LatticeTransformation.h"
#include "Pattern.h"
#include "Simple2DReducer.h"

namespace psimag {

  template <typename Field, size_t DIM, typename Algorithms> class Lattice;

  template <typename Field, size_t DIM, typename Occupant, typename Algorithms>  class Pattern;
  
  template <typename Field, size_t DIM, typename Occupant, 
	    typename LatticeType,typename Algorithms>  class PatternWithLattice;
  
  /** \ingroup latticeConcepts
   *
   * \brief The ReducedLattice class.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM Dimensionality of the lattice being represented. Legal values 
   *            are 1,2,3.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  class ReducedLattice: 
    public Lattice<Field,DIM,Algorithms> 
  {
  public:

    typedef Lattice<Field,DIM,Algorithms>     LatticeType;
    typedef LatticeTransformation<Field,DIM>  LatticeTransformationType;
    typedef CartesianTranslation<Field,DIM>   CartesianTranslationType;

  protected:
    
//     /**   
//      * \brief The given pattern is transformed to reflect the
//      *        reduction transformation to this Lattice.
//      */
//     template<typename Occupant, size_t NUMPOS>
//     void updatePattern(Pattern<Field,DIM,NUMPOS,Occupant,Algorithms>& origPattern,
// 		       Pattern<Field,DIM,NUMPOS,Occupant,Algorithms>& newPattern) {

//       typedef Pattern<Field,DIM,NUMPOS,Occupant,Algorithms> PatternType;
//       PatternType result(*this);
      
//       result[occupant] = transform.transform(cell)
//       pattern = transform(pattern);

//     }
    
  public:

    static std::string typeName() { return "ReducedLattice";   }

    static bool close(const Field v1, const Field v2) {
      return Algorithms::close(v1,v2);
    }

    /**
     * \brief The transformation which takes the origional cell into this reduced cell.
     */
    //    LatticeTransformationType transform; now inherited from lattice
    
    /**
     * \brief Construct a Reduced Lattice given an input/origional lattice.
     *
     * \note Reduced Lattice is initialized to the given origional cell
     *       and the identity transformation.
     */
    ReducedLattice():
      LatticeType()
    {
      //      throw std::logic_error("ReducedLattice() should not be called");
    }

    /**
     * \brief Construct a Reduced Lattice given an input/origional lattice.
     *
     * \note Reduced Lattice is initialized to the given origional cell
     *       and the identity transformation.
     */
    ReducedLattice(const LatticeType& origionalLattice):
      LatticeType(origionalLattice)
    {  }

    /**
     * \brief Perform a coordinated swap of the specified basis
     *        vectors and transform rows.
     *
     */
    void swapBasisVectors(size_t b1i, size_t b2i) {

      this->swapCols(b1i,b2i);
      this->transform.swapRows(b1i,b2i);
      this->update(); // update the metric and parameters
    }

    /**
     * \brief Find the multiple of the b1i vector that produces the
     *        shortest diagonal. Use this diagonal as a replacement for
     *        the b2i vector. 
     *
     * \note It is assumed that this cell is normalized at the thime
     *       of this call. That is to say: 1) the length of b1i is
     *       less that that of b2i and the angle between b1 and b2 is
     *       acute.
     */
    size_t shortestDiagonalFactor(size_t b1i, size_t b2i) {

      CartesianTranslationType b1    = (*this)[b1i];
      CartesianTranslationType b(b1);
      CartesianTranslationType b2    = (*this)[b2i];
      CartesianTranslationType lastDiag  = b - b2;

      Field b2Len = b2.length();

      // A private member
      Field lastDiagLen  = lastDiag.length();
      size_t factor = 1;
      
      while(true) {
	b = b + b1;  // b = factor*b1
	CartesianTranslationType diag = b - b2;
	Field dl = diag.length();
	if (dl > lastDiagLen) {
	  if(close(lastDiagLen, b2Len))
	    return 0;
	  if (lastDiagLen > b2Len)
	    return 0;
	  return factor;
	}
	factor++;
	lastDiag = diag;
	lastDiagLen = dl;
	if (factor > 100) 
	  throw std::range_error("reduceDiagonal failed!");
      }
    }

    /**
     * \brief Use the shortest diagonal between a multiple of the
     *        first vector and the second vector to relace the second
     *        vector with. Change the transform to reflect the
     *        change.  
     *
     * \note It is assumed that this cell is normalized at the time of
     *       this call and i < j. That is to say: 1) the length of bi
     *       is less that that of bj and the angle between bi and bj
     *       is acute.
     */
    bool reduceDiagonal(size_t bi, size_t bj) {
      
      typedef  typename LatticeTransformation<Field,DIM>::MatrixType::RotationType RMType;

      size_t factor = shortestDiagonalFactor(bi, bj);
      
      if (factor == Field(0)) return false;

      // Change the bj basis vector to f*bi - bj
      for (size_t r=0; r < DIM; r++)
	(*this)(r,bj) = ((*this)(r,bi) * Field(factor)) - (*this)(r,bj);

      // Change the transform to reflect this change
      RMType reduceDiagIncr;
      MakeIdentity(reduceDiagIncr);
      
      for (size_t c=0; c<DIM; c++) {
	if (c == bj) {
	  reduceDiagIncr(bj,c) = Field(-1);
	  continue;
	}
	if (c == bi) {
	  reduceDiagIncr(bj,c) = Field(factor);
	  continue;
	}
	reduceDiagIncr(bj,c) = Field(0);
      }
      RMType temp(this->transform.rotation); // Remember the part that we are about to overwrite.
      Multiply(reduceDiagIncr, temp, this->transform.rotation); //left multiply in the incremental transform

      this->update();

      return true;
    }
    
    /**
     * \brief Changing the direction of the second basis vector
     *        changes the angle from obtuse to acute.
     *
     * \note This could just as easily have been done with the first
     *       basis vector.
     */
    void makeAcute(size_t b1i, size_t b2i) { negate(b2i); }
    
    /**
     * \brief Changing the direction of the indicated
     *        changes the angle from obtuse to acute.
     */
    void negate(size_t r) {
      
      this->negateCol(r);
      this->transform.rotation.negateRow(r);
      this->update();
    }

 //    bool diagLessThanOrEqBasisVectors(size_t i, size_t j) {
      
//       CartesianTranslationType bi = (*this)[i];
//       CartesianTranslationType bj = (*this)[j];
//       CartesianTranslationType diagonal = bi-bj;
//       Field dl  = diagonal.length();
      
//       return (dl <= bi.length() && dl <= bj.length());
//     }
    
  };

  //======================================================================

  /** \ingroup ostream
   * Lattice output stream operator 
   **/
  template<typename Field, size_t DIM, typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const ReducedLattice<Field,DIM,Algorithms>& cell) {

    typedef typename Lattice<Field,DIM,Algorithms>::MatType LatticeMatType;
    
    //    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    //    os.precision(6);
    
    os << "ReducedLattice { det: "<< cell.det << ", "
       << "order: [" << cell.order << "], "
       << "params : "<< cell.parameters << ", " 
       << "type: " << cell.parameters.typeString() << "}" << std::endl;
    os << " ---------------------------------- Basis:"                << std::endl;
    os << ( (LatticeMatType) cell ) << std::endl;
    os << " ---------------------------------- End Basis:"  << std::endl;
    os << cell.metric << std::endl;
    os << cell.transform << std::endl;
    os << "----------------------------------- End Reduced Lattice " << std::endl;
    
    return os;
  }

  
}
#endif //PSIMAG_ReducedLattice_H

/*@}*/
