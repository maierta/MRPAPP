//-*-C++-*-

#ifndef PSIMAG_ConventionalLattice_H
#define PSIMAG_ConventionalLattice_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file  ConventionalLattice.h
 *  Contains class definitions for Conventional Lattice objects.
 */
 
#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"

#include "PSIMAGAssert.h"
#include "Pattern.h"
#include "Lattice.h"
#include "ReducedLattice.h"
#include "Centering.h"

namespace psimag {

  template<typename, size_t, typename> class ReducedLattice;

  template <typename Field, size_t DIM> class LatticeTransformation;

  /** \ingroup latticeConcepts
   *
   *  \brief A Cell of the type used in 2002 ITC to represent lattices
   *         and space groups.
   *
   *  Since many Cells can be chosen to represent the same lattice,
   *  crystallography defines conventional cells which are chosen in ways that make the
   *  specification of the symmetry operation's (i.e SymmetryOperation objects)
   *  and elements of the Lattice (or SpaceGroup) as simple as
   *  possible. 
   *
   *  The information in the International Tables for
   *  Crystallography (ITC) is provided in the context of such
   *  conventional cells. 
   *
   *  If a lattice's \b BravaisType has a Centering other than
   *  primitive it means that the the lattice basis of the
   *  corresponding ConventionalLattice is not primitive. In this case
   *  the basis vectors must be supplemented with additional Centering
   *  vectors before they can be used to generate the translational
   *  symmetry of the lattice.
   *
   *  For each dimension we will enumerate the of ConventionalLattices
   *  and assocaite them with lattice types.
   *
   *  \param Field: The scalar type used in the representation of the
   *                cell. Legal values include double, rational and
   *                sqrtExtendedRational.
   *
   *  \param DIM Dimensionality of the lattice being
   *             represented. Legal values are 1,2,3.
   */
  template <typename Field, size_t DIM, typename Algorithms>
  class ConventionalLattice: 
    public Lattice<Field,DIM,Algorithms>
  {
  public:
    
    typedef Lattice<Field,DIM,Algorithms>              LatticeType;
    typedef ReducedLattice<Field,DIM,Algorithms>       ReducedLatticeType;
    typedef LatticeTransformation<Field,DIM>           LatticeTransformationType;
    typedef CartesianTranslation<Field,DIM>            CartesianTranslationType;
    typedef ConventionalLattice<Field,DIM,Algorithms>  ConventionalLatticeType;

  protected:   
//     /**   
//      * \brief The given pattern was copied from the Reduced
//      *        Crystal.  It is transformed to reflect the
//      *        conventional transformation to this Lattice.
//      */
//     template<typename Occupant, size_t NUMPOS>
//     void updatePattern(Pattern<Field,DIM,NUMPOS,Occupant,Algorithms>& pattern) {

//     }
    
  public:
    
    static std::string typeName() { return "ConventionalLattice";   }

    /**
     * \brief The transformation which takes the reduced lattice into this lattice.
     */
    // LatticeTransformationType transform; now inherited from Lattice
    
    /**
     * \brief Construct the appropriate Conventional Lattice from the
     *        given reduced lattice. 
     *
     * The appropriate lattice (which may be the primitive reduced
     * lattice) will be the one which leads to a space group match.
     *
     * For each dimension there will be a set of standard,
     * primitive->conventional transformations. These standard
     * transformations will be evaluated (in a dimension and
     * Conventional Lattice Type specific manner?) to find those that
     * lead to a space group match (within a given delta). This
     * lattice is made to be the best one found and maintains a list
     * of alternates.
     *
     * \note This lattice is initialized with reduced lattice.C
     */
    ConventionalLattice(const ReducedLatticeType& reducedLattice): 
      LatticeType(reducedLattice)
    { 
      // "When Coded a ConventionalLattice is constructed from the given reducedLattice!" << std::endl;
    }

    /** The centering that is asocciated with this ConventionalLattice. */
    CenteringBase<Field,DIM> centering;
    
  };

  /** \ingroup ostream
   * Lattice output stream operator 
   **/
  template<typename Field, size_t DIM,typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const ConventionalLattice<Field,DIM,Algorithms>& cell) {

    typedef typename ConventionalLattice<Field,DIM,Algorithms>::MatType LatticeMatType;
    
    //    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    //    os.precision(6);
    
    os << "ConventionalLattice { det: "<< cell.det << ", "
       << "order: [" << cell.order << "], "
       << "params : "<< cell.parameters << ", " 
       << "type: " << cell.parameters.typeString() << "}" << std::endl;
    os << " ---------------------------------- Basis:"                << std::endl;
    os << ( (LatticeMatType) cell ) << std::endl;
    os << " ---------------------------------- End Basis:"  << std::endl;
    os << cell.metric << std::endl;
    os << cell.transform << std::endl;
    os << "----------------------------------- End ConventionalLattice " << std::endl;
    
    return os;
  }

} /* namespace psimag */

#endif

/*@}*/

