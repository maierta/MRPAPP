//-*-C++-*-

#ifndef  Psimag_APPLIED_SYMMETRY_ELEMENT
#define  Psimag_APPLIED_SYMMETRY_ELEMENT

/** \ingroup symmetryConcepts */
/*@{*/

/** \file AppliedSymmetryElement.h
 *
 *  Contains a class for implementing symmetry elements.
 *
 */
 
#include <iostream>
#include <sstream>
#include <string>

#include "Vec.h"
#include "Mat.h"
#include "PSIMAGAssert.h"

#include "FieldConvert.h"

#include "CartesianTranslation.h"
#include "CartesianPosition.h"
#include "SymmetryOperation.h"
#include "TestPattern.h"
#include "LatticeWithPattern.h"
#include "SymmetryElement.h"
#include "SymmetryElements2D.h"

namespace psimag {
  
  template<typename Field, size_t DIM, typename Algorithms>
  class SymmetryElements;

  /** \ingroup symmetryConcepts
   *
   * \brief A class for implementing SymmetryElements which record
   *        their application to a pattern.
   *
   * When a SymmetryElement is applied to a LatticeWithPattern we can construct:
   * 
   * - A SymmetryOperation which represents the SymmetryElement applied
   *   to a given Lattice,
   *
   * - A TestPattern which is the image of the given Pattern under this
   *   objects SymmetryOperation.
   */
  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class AppliedSymmetryElement
  {
    
  public:

    typedef SymmetryElement<Field,DIM,Algorithms>                  SymmetryElementType;
    typedef SymmetryElements<Field,DIM,Algorithms>                 SymmetryElementsType;

    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;

    typedef SymmetryOperation<Field,DIM,Algorithms>                SymmetryOperationType;
    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms> TestPatternType;

    typedef CartesianTranslation<Field,DIM>                        CartesianTranslationType;
    typedef CartesianPosition<Field,DIM>                           CartesianPositionType;


    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms> LatticeWithPatternType;

    SymmetryElement< Field, DIM, Algorithms>  element;

    CartesianTranslationType                  cartesianTranslation;
    CartesianPositionType                     cartesianPosition;

    SymmetryOperationType                     operation;

    TestPatternType                           testPattern;
    Field                                     latticeDistance;
    Field                                     patternDistance;
    Field                                     distance;
    int                                       id;

    //====================================================================== Constructors
    /**
     * \brief  The default constructor sets to identity
     */
    AppliedSymmetryElement(): 
      element              (),
      cartesianTranslation (),
      cartesianPosition    (),
      operation            (),
      testPattern          (),
      latticeDistance      (Field(0)),
      patternDistance      (Field(0)),
      distance             (Field(0)),
      id                   (-1)
    {}

    /**
     * \brief  The copy constructor 
     */
    AppliedSymmetryElement(const AppliedSymmetryElement& other): 
      element              (other.element), 
      cartesianTranslation (other.cartesianTranslation),
      cartesianPosition    (other.cartesianPosition),
      operation            (other.operation),
      testPattern          (other.testPattern),
      latticeDistance      (other.latticeDistance),
      patternDistance      (other.patternDistance),
      distance             (other.distance),
      id                   (other.id)
    {   
//       if(element.id == 0 && element.trace == convert<Field>(-911)) {
// 	std::ostringstream buff;
// 	buff << "Error in AppliedSymmetryElement Copy Constructor " << std::endl;
// 	buff << toXML(*this) << std::endl;
// 	buff << toXML(other) << std::endl;
// 	throw std::logic_error(buff.str());
//       }
    }

    template< typename OtherLatticeType >
    /**
     * \brief  The cross lattice type copy constructor 
     */
    AppliedSymmetryElement(const AppliedSymmetryElement<Field,DIM,Occupant,OtherLatticeType,Algorithms>& other): 
      element              (other.element), 
      cartesianTranslation (other.cartesianTranslation),
      cartesianPosition    (other.cartesianPosition),
      operation            (other.operation),
      testPattern          (other.testPattern),
      latticeDistance      (other.latticeDistance),
      patternDistance      (other.patternDistance),
      distance             (other.distance),
      id                   (other.id)
    {
//       if(element.trace == convert<Field>(-911)) {
// 	std::ostringstream buff;
// 	buff << "Error in AppliedSymmetryElement Other lattice Copy Constructor " << std::endl << toXML(*this) << std::endl;
// 	throw std::logic_error(buff.str());
//       }
    }
    
    /**
     * \brief The constructor which applies a given SymmetryElement
     *        (and its associated SymmetryOperation) to a
     *        latticeWithPattern.
     */
    AppliedSymmetryElement(int index,
			   const LatticeWithPatternType& latpat,
			   const SymmetryElementType&    symEl): 
      element             (symEl), 
      cartesianTranslation(latpat.cartesianTranslation(element.netDirection)),
      cartesianPosition   (latpat.cartesianPosition   (element.cellPosition)  ),
      operation           (symEl(static_cast< Lattice<Field,DIM,Algorithms> >(latpat))),
      testPattern         (latpat,this->operation),
      latticeDistance     (latpat.symmetryDistance(this->operation)),
      patternDistance     (this->testPattern.distanceTo(latpat.pattern)),
      distance            (std::max(this->latticeDistance, this->patternDistance)),
      id                  (index)
    {	
//       if(element.trace == convert<Field>(-911)) {
// 	std::ostringstream buff;
// 	buff << "Error in AppliedSymmetryElement latpat symEl Constructor " << std::endl << toXML(*this) << std::endl;
// 	throw std::logic_error(buff.str());
//       }
      circleCheck();
    }

    //======================================================================
    //
    bool circleCheck() {
      return element.circleCheck(operation, testPattern);
    }

    //======================================================================
    //
    /**
     * \brief AppliedSymmetryElelment Assignment Operator 
     * 
     */
    AppliedSymmetryElement& operator = (const AppliedSymmetryElement& op) {

      cartesianTranslation = op.cartesianTranslation;
      cartesianPosition    = op.cartesianPosition;
      operation            = op.operation;
      testPattern          = op.testPattern;
      latticeDistance      = op.latticeDistance;
      patternDistance      = op.patternDistance;
      distance             = op.distance;
      id                   = op.id;
      element              = op.element;

//       if(element.trace == convert<Field>(-911)) {
// 	std::ostringstream buff;
// 	buff << "Error in AppliedSymmetryElement operator = " << std::endl;
// 	buff << toXML(*this) << std::endl;
// 	buff << toXML(op) << std::endl;
// 	throw std::logic_error(buff.str());
//       }
      return (*this);
    }
  };

  //====================================================================== 
  /** \ingroup XML
   *
   * SpaceGroup  XML output
   **/
  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms> 
  Tag toXML(const AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms>& symEl,
	    std::string name="AppliedSymmetryElement") {

    Tag result(name);

    result["latticeDistance"]  = symEl.latticeDistance;
    result["patternDistance"]  = symEl.patternDistance;
    result["distance"]         = symEl.distance;
    result["id"]               = symEl.id;

    result.add(toXML(symEl.element));
    result.add(toXML(symEl.cartesianTranslation,"CartesianDirection"));
    result.add(toXML(symEl.cartesianPosition   ,"CartesianPosition"));
    
    const LatticeType& lat = symEl.testPattern;
    result.add(toXML(symEl.operation,lat));
    //result.add(toXML(symEl.testPattern,"TestPattern"));

    return result;

  }


  //====================================================================== 
  /** \ingroup ostream
   *
   * \brief  SymmetryOperation output stream operator.
   *
   */
  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms> 
  std::ostream& operator << (std::ostream& os, 
			     AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms>& S) {
      
    typedef SymmetryElement<Field,DIM,Algorithms>   SymmetryElementType;
    SymmetryElementType&       asElement = S;

    os << "Applied Symmetry Element:  " << S.name << "------------------------------" << std::endl;
    os << asElement  << std::endl;
    os << "End Applied Symmetry Element ------------------------------" << std::endl;
    return os;
  }

} /* namespace psimag */

#endif   //Psimag_SYMMETRY_OPERATION

/*@}*/
