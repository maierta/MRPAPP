//-*-C++-*-

#ifndef  Psimag_SYMMETRY_ELEMENT
#define  Psimag_SYMMETRY_ELEMENT

/** \ingroup symmetryConcepts */
/*@{*/

/** \file Symmetryelement.h
 *
 *  Contains a class for implementing symmetry elements.
 *
 *  \author  Mike Summers
 *
 */
 
#include <iostream>
#include <sstream>
#include <string>

#include "Vec.h"
#include "Mat.h"
#include "PSIMAGAssert.h"

#include "FieldConvert.h"

#include "SeitzMatrix.h"
#include "SeitzVectors.h"
#include "CellPosition.h"
#include "LatticeCoordinates.h"
#include "Lattice.h"
#include "SymmetryOperation.h"
#include "SymmetryElementName.h"
#include "CellTranslation.h"
#include "CellRotation.h"

namespace psimag {
  
  template<typename Field, size_t DIM, typename Algorithms> class SymmetryOperation;

  // These will be specialized elsewhere
  template<typename Field, size_t DIM, typename Algorithms> class IdentityElement {};
  template<typename Field, size_t DIM, typename Algorithms> class Mirror    {};
  template<typename Field, size_t DIM, typename Algorithms> class Glide     {};
  template<typename Field, size_t DIM, typename Algorithms> class TwoFold   {};
  template<typename Field, size_t DIM, typename Algorithms> class ThreeFold {};
  template<typename Field, size_t DIM, typename Algorithms> class FourFold  {};
  template<typename Field, size_t DIM, typename Algorithms> class SixFold   {};
  template<typename Field, size_t DIM, typename Algorithms> class TwoFoldN  {};
  template<typename Field, size_t DIM, typename Algorithms> class ThreeFoldN{};
  template<typename Field, size_t DIM, typename Algorithms> class FourFoldN {};
  template<typename Field, size_t DIM, typename Algorithms> class SixFoldN  {};


  //======================================================================
  /** 
   * \brief Return the value of the trace that is expected for a
   *        given type of SymmetryElement.  
   */
  template<size_t DIM>
  int expectedTrace(const std::string& type) {
    
    if (type=="identity")   return  2;
    if (type=="translation")return  2;
    if (type=="twoFold")    return -2;
    if (type=="threeFold")  return -1;
    if (type=="fourFold")   return  0;
    if (type=="sixFold")    return  1;
    if (type=="threeFoldN") return -1;
    if (type=="fourFoldN")  return  0;
    if (type=="sixFoldN")   return  1;
    if (type=="glide")      return  0;
    if (type=="mirror")     return  0;
    
    std::ostringstream buff;
    buff << "SymmetryElement.h DIM="<<DIM<<" expectedTrace("<<type<<") unknown type!";
    throw std::logic_error(buff.str());
  }

  //======================================================================
  /** 
   * \brief Return the value of the trace that is expected for a
   *        given type of SymmetryElement.  
   */
  template<typename Field, size_t DIM> 
  Field expectedRotationAngle(const std::string& type) {

    if (type=="identity")   return  Field(-911);
    if (type=="translation")return  Field(-911);
    if (type=="twoFold")    return  Field(180);
    if (type=="threeFold")  return  Field(120);
    if (type=="fourFold")   return  Field(90);
    if (type=="sixFold")    return  Field(60);
    if (type=="threeFoldN") return  Field(-120);
    if (type=="fourFoldN")  return  Field(-90);
    if (type=="sixFoldN")   return  Field(-60);
    if (type=="glide")      return  Field(-911);
    if (type=="mirror")     return  Field(-911);
    
    std::ostringstream buff;
    buff << "SymmetryElement.h DIM="<<DIM<<" expectedRotationAngle("<<type<<") unknown type!";
    throw std::logic_error(buff.str());
  }

  //======================================================================
  /** 
   * \brief Return the value of the trace that is expected for a
   *        given type of SymmetryElement.  
   */
  template<size_t DIM> 
  int expectedSense(const std::string& type) {

    if (type=="identity")   return   0;
    if (type=="translation")return   0;
    if (type=="twoFold")    return   0;
    if (type=="threeFold")  return   1;
    if (type=="fourFold")   return   1;
    if (type=="sixFold")    return   1;
    if (type=="threeFoldN") return  -1;
    if (type=="fourFoldN")  return  -1;
    if (type=="sixFoldN")   return  -1;
    if (type=="glide")      return   0;
    if (type=="mirror")     return   0;
    
    std::ostringstream buff;
    buff << "SymmetryElement.h DIM="<<DIM<<" expectedSense("<<type<<") unknown type!";
    throw std::logic_error(buff.str());
  }

  /** \ingroup symmetryConcepts
   *
   * \brief A class for implementing symmetry elements.
   *
   * Symmetry elements store the geometric characteristics of a
   * SymmetryOperation.
   *
   */

  template<typename Field, size_t DIM, typename Algorithms>
  class SymmetryElement {
    
  public:

    typedef CellTranslation<Field,DIM>              CellTranslationType;
    typedef CellPosition<Field,DIM,Algorithms>      CellPositionType;
    typedef LatticeCoordinates<DIM>                 LatticeCoordinatesType;
    typedef SymmetryElement<Field,DIM,Algorithms>   SymmetryElementType;
    typedef SymmetryOperation<Field,DIM,Algorithms> SymmetryOperationType;

    std::string             type;
    std::string             name;
    std::string             unNormalizedName;

    /* \brief The direction of linear symmetry elements given in
     *        lattice coordinates.
     *
     *  For mirrors, this is the direction of a mirror line in lattice
     *  coordinates, the coordinates are always chosen so that the
     *  first coordinate is positive.
     *
     *  For glides and translations this is the direction of the glide
     *  or translation. The translationPart is parallel to the
     *  netDirection but is usually not the same length and may not be
     *  in the same direction.
     *
     *  In 3D this is used to give the axis of rotation.
     */
    LatticeCoordinatesType  netDirection;

    /* \brief The translation part of a screw, glide or translation.
     */
    CellTranslationType     translationPart;

    /* The glide fraction is given by netDirection*glideFraction = translationPart
     */
    Field                   glideFraction;

    CellPositionType        cellPosition;

    Field                   trace;
    Field                   rotationAngle;
    int                     sense;
    int                     id;

    virtual ~SymmetryElement() {}

    //====================================================================== Constructors
    /**
     * \brief  The default constructor sets to no-name
     */
    SymmetryElement(): 
      type("no-type"),
      name("no-name"),
      unNormalizedName("no-name"),
      glideFraction(convert<Field>(0)),
      trace        (convert<Field>(-911)),
      rotationAngle(convert<Field>(-911)),
      sense        (-911),
      id(-1)
    { }

    /**
     * \brief  The generic constructor 
     */
    SymmetryElement(std::string                   elType,
		    const LatticeCoordinatesType& orientation,
		    const CellPositionType&       pos,
		    Field                         glide=convert<Field>("1/2")): 
      type            (elType),
      name            ("no-name?"),
      unNormalizedName("no-name?"),
      netDirection    (orientation),
      translationPart (orientation.cellTranslation<Field>() * glide),
      glideFraction   (glide),
      cellPosition    (pos),
      trace           (expectedTrace<DIM>(elType)),
      rotationAngle   (expectedRotationAngle<Field,DIM>(elType)),
      sense           (expectedSense<DIM>(elType)),
      id              (-1)
    { 
      unNormalizedName = SymmetryElementName<DIM>(*this);
      normalize(*this);
      name = SymmetryElementName<DIM>(*this);
    }

    /**
     * \brief  The copy constructor 
     */
    SymmetryElement(const SymmetryElementType& other): 
      type            (other.type),
      name            (other.name),
      unNormalizedName(other.unNormalizedName),
      netDirection    (other.netDirection),
      translationPart (other.translationPart),
      glideFraction   (other.glideFraction),
      cellPosition    (other.cellPosition),
      trace           (other.trace),
      rotationAngle   (other.rotationAngle),
      sense           (other.sense),
      id              (other.id)
    {}

    //======================================================================

    SymmetryElement& operator = (const SymmetryElement& el) {
      type             = el.type;
      name             = el.name;
      unNormalizedName = el.unNormalizedName;
      netDirection     = el.netDirection;
      cellPosition     = el.cellPosition;
      translationPart  = el.translationPart;
      trace            = el.trace;
      rotationAngle    = el.rotationAngle;
      sense            = el.sense;
      glideFraction    = el.glideFraction;
      id               = el.id;
      return (*this);
    }
    
    //======================================================================
    /** 
     * \brief Throw a logic error if the netDirection, cellTranslation,
     *        glideFraction combo is not right.
     */
    void checkDirection() const {

      static Field zero(0);

      if (type == "translation" || type == "glide") {

	for (size_t i=0; i<DIM; i++) {
	  Field comparator = Field(netDirection[i])*glideFraction;
	  if (! (Algorithms::close(comparator,
				   translationPart[i])) ) {
	    std::ostringstream buff;
	    buff << "In SymmetryElement.checkDirection, error! " << std::endl
		 << " i          = " << i << std::endl
		 << " comparator = " << comparator << std::endl
		 << " element    = " << std::endl << toXML(*this) ;
	    throw std::logic_error(buff.str());
	  }
	}
	return;
      }
      
      if (! (Algorithms::close(translationPart[0], zero) &&
	     Algorithms::close(translationPart[0], zero) ) ) {
	std::ostringstream buff;
	buff << "In SymmetryElement.checkDirection error! " << std::endl 
	     << "   glideFraction error! non-zero translationPart " 
	     << " element = " << std::endl << toXML(*this);
	throw std::logic_error(buff.str());
      }  

      if (! Algorithms::close(glideFraction, zero) ) {
	std::ostringstream buff;
	buff << "In SymmetryElement.checkDirection error! " << std::endl 
	     << "   glideFraction should be zero for " 
	     << " element = " << std::endl << toXML(*this);
	throw std::logic_error(buff.str());
      }
    }
    
    //======================================================================
    /**
     * \brief Sets the glide fraction from the netDirection and the translationPart.
     *
     */
    void setGlideFraction() {

      static Field zero(0);
      glideFraction = zero;

      if (! (type == "glide"       ||
	     type == "translation" ||
	     type == "screw" ) ) return;

      if (netDirection[0] != 0) 
	glideFraction = translationPart[0] / Field(netDirection[0]);
      else if (netDirection[1]  != 0)
	glideFraction = translationPart[1] / Field(netDirection[1]);
    }

    //======================================================================
    /** 
     * \brief Return the symmetry operation corresponding to this
     *        element in a given lattice.
     */
    virtual SymmetryOperationType operator () (const Lattice<Field,DIM,Algorithms>& lattice) const
    {
      return operationFor(*this, lattice);
    }

    bool circleCheck(const Lattice<Field,DIM,Algorithms>& lattice) const {
      return circleCheck(operationFor(*this, lattice), lattice);
    }

    bool circleCheck(const SymmetryOperationType&         operation, 
		     const Lattice<Field,DIM,Algorithms>& lattice) const {
      
      SymmetryElementType   otherEl = symmetryElement(operation, lattice);

      if( otherEl == *this) return true;
      
      std::ostringstream buff;
      buff << "---------------------------------  CircleCheck: Symmetry Elements not == " << std::endl;
      diff(buff,*this,otherEl);
      buff << "---------------------------------  original:" << std::endl;
      buff << toXML(*this)      << std::endl;
      buff << toXML(operation)      << std::endl;
      buff << "---------------------------------  from generated operation:" << std::endl;
      buff << toXML(otherEl) << std::endl;
      throw std::logic_error(buff.str());
      
      return false;
    }

    //======================================================================
    /** 
     * \brief Return the symmetry operation corresponding to this
     *        element in a given lattice.
     */
    static SymmetryOperationType operationFor (const SymmetryElementType& element,
					       const Lattice<Field,DIM,Algorithms>& lattice)
    {
      if (element.type == "identity") {
	return IdentityElement<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "translation") {
	SymmetryOperationType result = IdentityElement<Field,DIM,Algorithms>::operationFor(element,lattice);
	result.setTranslation(lattice.cartesianTranslation(element.translationPart));
	result.name = element.name;
	return result;
      }
      if (element.type == "twoFold") {
	return TwoFold<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "threeFold") {
	return ThreeFold<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "threeFoldN") {
	return ThreeFoldN<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "fourFold") {
	return FourFold<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "fourFoldN") {
	return FourFoldN<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "sixFold") {
	return SixFold<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "sixFoldN") {
	return SixFoldN<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "glide") {
	return Glide<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      if (element.type == "mirror") {
	return Mirror<Field,DIM,Algorithms>::operationFor(element,lattice);
      }
      std::ostringstream buff;
      buff << "In SymmetryElement.oparationFor, SymmetryElement type (=" << element.type << ") not found!";
      buff << "   " << toXML(element) << std::endl;
      throw std::logic_error(buff.str());
    }
  };

  //======================================================================
  /** 
   * \Brief Normalize the components of this symmetry element.
   *
   * \note This is a 2D function only!
   */
  template<typename Field, typename Algorithms> 
  void normalize(SymmetryElement<Field,2,Algorithms>& element) {   

    enum {DIM=2};
    typedef CellTranslation<Field,DIM> CellTranslationType;
    
    static Field zero(0);
    //static Field one(1);
    CellTranslationType zeroTranslation;
    
    //Normalize the position
    element.cellPosition.normalize();
    
    //Normalize the translation part 
    if (element.type == "translation") 
      normalizeTranslation<Field,DIM,Algorithms>(element.translationPart);
    
    if (element.type == "glide")
      normalizeGlide<Field,DIM,Algorithms>(element.translationPart);

    // If the normalized translationPart is zero then a translation is
    // has become the identity and a glide has become a mirror
    
    bool isZeroTranslation = closeTo<Field,DIM,0,Algorithms>(element.translationPart, zeroTranslation);
    
    if (isZeroTranslation) {

      if (element.type == "translation") 
	element.type    = "identity";
      
      if (element.type == "glide") {
	element.type    = "mirror";
	element.glideFraction = zero;
      }
      if (element.type == "mirror")
	element.netDirection = element.netDirection.normalize<Algorithms>();

      if (element.type == "identity") 
	element.netDirection = latticeCoordinate(0,0);
    }
    else {  // not isZeroTranslation
      // In 2D only glides and translations have translation parts
      if (element.type != "glide" && element.type != "translation" ) {
	std::ostringstream buff;
	buff << "SymmetryElement.normalize(" << element.unNormalizedName << ")" << std::endl;
	buff << " type is now " << element.type << " but the tranlationPart is:" << element.translationPart << std::endl;
	throw std::logic_error(buff.str());
      }
      // Set the netDirection for the glide from the transltionPart
      element.netDirection = latticeCoordinate<Field,Algorithms>(element.translationPart);
      // Set the glideFraction
      element.setGlideFraction();
    }

    if (element.type == "mirror" || element.type == "glide")
      // The location point was given on the b axis 
      if (element.cellPosition[0] == zero &&
	  element.cellPosition[0] != zero ) {  
	
	Field slope = getSlope<Field>(element.netDirection);
	element.cellPosition = cellPosition<Field,Field>((slope / element.cellPosition[1]), 0);
	// now it should be on the a axis!
      }
  }

  //====================================================================== 
  /** \ingroup XML
   *
   * SpaceGroup  XML Symmetry Element Output
   **/
  template<typename Field, size_t DIM, typename Algorithms> 
  Tag toXML(const SymmetryElement<Field,DIM,Algorithms>& el,
	    std::string name="SymmetryElement") {

    Tag result(name);
    result["type"]             = el.type;
    result["name"]             = el.name;
    result["unNormalizedName"] = el.unNormalizedName;
    result["DIM"]              = DIM;
    result["glideFraction"]    = el.glideFraction;
    result["id"]               = el.id;
    result["trace"]            = el.trace;

    if(el.type != "glide" && 
       el.type != "mirror" && 
       el.type != "translation" &&
       el.type != "identity") {
      result["rotationAngle"]  = el.rotationAngle;
      result["sense"]          = el.sense;
    }

    result.add(toXML(el.netDirection,    "NetDirection"));
    result.add(toXML(el.cellPosition,    "CellPosition"));
    result.add(toXML(el.translationPart, "TranslationPart"));

    return result;
  }

  //====================================================================== 
  /** \ingroup ostream
   *
   * \brief  Symmetryelement output stream operator.
   *
   */
  template<typename Field,size_t DIM, typename Algorithms>
  bool operator == (const SymmetryElement<Field,DIM,Algorithms>& el,
		    const SymmetryElement<Field,DIM,Algorithms>& otherEl) {
      
    if (el.type            == otherEl.type &&
	el.name            == otherEl.name &&
	Algorithms::close(el.trace, otherEl.trace) &&
	el.sense           == otherEl.sense &&
	el.netDirection    == otherEl.netDirection &&
	el.cellPosition.closeTo(otherEl.cellPosition) &&
	closeTo<Field,DIM,0,Algorithms>(el.translationPart,otherEl.translationPart))
      return true;

    return false;
  }

  //====================================================================== 
  /** \ingroup ostream
   *
   * \brief  Symmetryelement output stream operator.
   *
   */
  template<typename Field,size_t DIM, typename Algorithms>
  void diff (std::ostream& os, 
	     const SymmetryElement<Field,DIM,Algorithms>& el,
	     const SymmetryElement<Field,DIM,Algorithms>& otherEl) {
    
  } 
  
  //====================================================================== 
  /** \ingroup ostream
   *
   * \brief  Symmetryelement output stream operator.
   *
   */
  template<typename Field,size_t DIM, typename Algorithms>
  std::ostream& operator << (std::ostream& os, const SymmetryElement<Field,DIM,Algorithms>& el) {
      
    os << "SymmetryElement[DIM=" << DIM << ", glideFraction=" << el.glideFraction << "] " << el.type << "---------" << std::endl;
    os << "net direction: " << el.direction << std::endl;
    os << "cell cellPosition:    " << el.cellPosition    << std::endl;
    os << "End Symmetry Element ------------------------------" << std::endl;
    return os;
  }

} /* namespace psimag */

#endif   //Psimag_SYMMETRY_OPERATION

/*@}*/
