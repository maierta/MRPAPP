//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file SuperCrystal.h  
 *
 *  Contains the SuperCrystal class.
 */

#ifndef Psimag_SuperCrystal_H
#define Psimag_SuperCrystal_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"
#include "Matrix.h"
#include "XMLHeading.h"

#include "Lattice.h"
#include "ReciprocalLattice.h"
#include "Pattern.h"
#include "PatternData.h"
#include "TestPattern.h"
#include "CrystalBase.h"
#include "Reciprocal.h"
#include "FloodTiler.h"
#include "SuperCrystalBuilder.h"

namespace psimag {

  //======================================================================

  /** \ingroup crystallography
   *  
   *\brief SuperCrystal class which is built from a crystal structure and a sublattice specification.
   *
   *  The super crystal has 4 main parts:
   *
   *  - A given origional crystal
   *  - the reciprocal Lattice of the origional
   *  - A 'Super Crystal'
   */ 
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms, 
	   template<typename, size_t, typename, typename> class BuilderHelperTemplate = SuperCrystalBuilder>
  class SuperCrystal: 
    public CrystalBase<Field,DIM,Occupant,Lattice,Algorithms>
  // not Crystal<Field,DIM,Occupant,Algorithms> since we don't want to consider reductions etc.
  {
  public:

    typedef SuperCrystal<Field,DIM,Occupant,Algorithms,BuilderHelperTemplate> ThisType;      

    typedef BuilderHelperTemplate<Field,DIM,Occupant,Algorithms>          BuilderHelperType;

    typedef CrystalBase<Field,DIM,Occupant,Lattice,Algorithms>            BaseType;
    typedef Crystal<Field,DIM,Occupant,Algorithms>                        CrystalType;
    typedef Lattice<Field,DIM,Algorithms>                                 LatticeType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                        PatternType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;

    typedef ReciprocalLattice<Field,DIM,Lattice,Algorithms>               ReciprocalLatticeType;
    typedef FloodTiler<Field,DIM,Occupant,Lattice,Algorithms>             FloodTilerType;
    typedef LatticeCoordinates<DIM>                                       SubLatticeVecCoordType;
    typedef std::vector<SubLatticeVecCoordType>                           SubLatticeVecCoordsType;
    SubLatticeVecCoordsType        subLatVecCoords;
    const CrystalType&             crystal;                       /**  The Crystal that this is a SuperCrystal of. */
    LatticeType                    reciprocalCrystalLattice;      /**  The reciprocal lattice of the crystal that this is a SuperCrystal of. */
    LatticeType                    reciprocalSuperCrystalLattice; /**  The reciprocal lattice of this SuperCrystal's lattice. */
    BaseType                       reciprocalSuperCrystal;        /**  The reciprocal of this SuperCrystal. */
    LatticeWithPatternType         reciprocalCrystalPattern;      /**  The superCrystal generated pattern for the 'reciprocal' crystal. */
    CrystalType                    reciprocalCrystal;             /**  The 'reciprocal' crystal. */

    /** \ingroup crystallography
     *  
     *\brief SuperCrystal constructor.
     *
     */ 
    SuperCrystal(SubLatticeVecCoordsType& subLatticeVecCoords, 
		 const CrystalType&       cryst,
		 PatternType&             reciprocalSuperCrystalPattern): 
      BaseType( FloodTilerType(subLatticeVecCoords, cryst).getTiledLatticeWithPattern() ),
      subLatVecCoords(subLatticeVecCoords),
      crystal(cryst),                   // Save a reference to the input crystal structure
      reciprocalCrystalLattice(ReciprocalLatticeType(crystal)),
      reciprocalSuperCrystalLattice(ReciprocalLatticeType(*this)),
      reciprocalSuperCrystal( reciprocalSuperCrystalLattice, reciprocalSuperCrystalPattern),
      reciprocalCrystalPattern( FloodTilerType(reciprocalCrystalLattice, reciprocalSuperCrystal).getTiledLatticeWithPattern()),
      reciprocalCrystal(reciprocalCrystalPattern)
    {
      this->analyzeSpaceGroups();
      reciprocalSuperCrystal.analyzeSpaceGroups();
    }

    /** \ingroup crystallography
     *  
     *\brief SuperCrystal constructor.
     *
     */ 
    SuperCrystal(const CrystalType&       cryst,
		 const BuilderHelperType& builder): 
      BaseType( FloodTilerType(builder.getSubLatticeVecCoords(), cryst).getTiledLatticeWithPattern() ),
      subLatVecCoords(builder.getSubLatticeVecCoords()),
      crystal(cryst),                   // Save a reference to the input crystal structure
      reciprocalCrystalLattice      (ReciprocalLatticeType(crystal).getLattice()),
      reciprocalSuperCrystalLattice (ReciprocalLatticeType(*this  ).getLattice()),
      reciprocalSuperCrystal( reciprocalSuperCrystalLattice, builder.getReciprocalSuperCrystalPattern()),
      reciprocalCrystalPattern(FloodTilerType(reciprocalCrystalLattice, reciprocalSuperCrystal).getTiledLatticeWithPattern()),
      reciprocalCrystal(reciprocalCrystalPattern)
    {
      this->symmetry.setAppliedSymmetryElements();
      this->symmetry.analyzeSpaceGroups();
      reciprocalSuperCrystal.symmetry.setAppliedSymmetryElements();
      reciprocalSuperCrystal.symmetry.analyzeSpaceGroups();
    }
    
    template<typename ObjectType>
    void dumpXML(const std::string outputDir, const std::string name, const ObjectType& obj) const {
      
      std::ofstream xmlFile;
      xmlFile.open((outputDir+name+".xml").c_str());
      xmlFile << XMLHeading(XMLHeading::XHTML) << std::endl;
      xmlFile << toXML(obj)  << std::endl;
      xmlFile.close();

    }

    void dumpXML(std::string directoryName) const {
      dumpXML<ThisType>   (directoryName, "SuperCrystal",         *this);
      dumpXML<CrystalType>(directoryName, "GivenCrystal",          this->crystal);
      dumpXML<BaseType>   (directoryName, "ReciprocalCrystal",     this->reciprocalCrystal);
      dumpXML<BaseType>   (directoryName, "ReciprocalSuperCrystal",this->reciprocalSuperCrystal);
    }

  };
    
  /** \ingroup ostream
   *
   * Output stream operator for SuperCrystal
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms, 
	   template<typename, size_t, typename, typename> class BuilderHelperTemplate>
  Tag toXML(const SuperCrystal<Field,DIM,Occupant,Algorithms,BuilderHelperTemplate>& superCrystal,
	    std::string name="SuperCrystal") {

    typedef CrystalBase<Field,DIM,Occupant,Lattice,Algorithms>  BaseType;

    Tag tag = toXML( (BaseType)superCrystal, "SuperCrystal" );

    Tag subLatticeCoordsTag("SubLatticeCoords");
    for(size_t i=0; i< superCrystal.subLatVecCoords.size(); i++) 
      subLatticeCoordsTag.add(toXML(superCrystal.subLatVecCoords[i]));

    tag.add(subLatticeCoordsTag);
    //tag.add(toXML(superCrystal.crystal,               "GivenCrystal"));
    //tag.add(toXML(superCrystal.reciprocalCrystal,     "ReciprocalCrystal"));
    //tag.add(toXML(superCrystal.reciprocalSuperCrystal,"ReciprocalSuperCrystal"));

    return tag;
  }

  /** \ingroup ostream
   *
   * Output stream operator for SuperCrystal
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms, 
	   template<typename, size_t, typename, typename> class BuilderHelperTemplate>
  inline  std::ostream& operator << (std::ostream& os, 
				     const SuperCrystal<Field,DIM,Occupant,Algorithms,BuilderHelperTemplate>& superCrystal) {

    typedef CrystalBase<Field,DIM,Occupant,Lattice,Algorithms>            BaseType;
    
    os << "====================================================================== SuperCrystal" << std::endl;
    os << "====================================================================== ^^^^^^^^^^^^" << std::endl;
    os << "==============================================****************************************** Origional Crystal:" << std::endl;
    os << superCrystal.crystal << std::endl;
    os << "==============================================****************************************** This SuperCrystal:" << std::endl;
    os << (BaseType) superCrystal    << std::endl;
    os << "==============================================****************************************** Origional Crystal's recipricol Lattice:" << std::endl;
    os <<  superCrystal.reciprocalCrystalLattice    << std::endl;
    os << "==============================================****************************************** Super Crystal's recipricol Lattice:" << std::endl;
    os <<  superCrystal.reciprocalSuperCrystalLattice    << std::endl;
    os << "==============================================****************************************** Reciprocal SuperCrystal:" << std::endl;
    os << superCrystal.reciprocalSuperCrystal   << std::endl;
    os << "==============================================****************************************** Generated Pattern for the reciprocal crystal:" << std::endl;
    os << superCrystal.reciprocalCrystalPattern << std::endl;

    os << "==============================================********************************************* Generated reciprical crystal:" << std::endl;
    os << superCrystal.reciprocalCrystal << std::endl;
    return os;
  }
} /** namespace psimag */

#endif
/*@}*/
