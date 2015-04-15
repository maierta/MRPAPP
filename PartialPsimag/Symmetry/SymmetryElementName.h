//-*-C++-*-

#ifndef  Psimag_SYMMETRY_ELEMENT_NAME
#define  Psimag_SYMMETRY_ELEMENT_NAME

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SymmetryElementName.h
 *
 *  Contains a class for implementing the naming of symmetry elements and the parsing of names to elements.
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

#include "SymmetryElement.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Algorithms>
  class SymmetryElement;
  

  template<typename Field, typename Algorithms>
  class Fraction {
  public:
    static std::string GET(Field value) {
  
      static std::string negStr                 ("-");
      static std::string noStr                  ("");

      static std::string zeroStr                ("0");
      static std::string oneHalfStr             ("1/2");
      static std::string oneThirdStr            ("1/3");
      static std::string oneQuarterStr          ("1/4");
      static std::string oneStr                 ("1");
      static std::string twoThirdsStr           ("2/3");
      static std::string threeQuartersStr       ("3/4");
      static std::string fiveQuartersStr        ("5/4");

      static std::string minusOneHalfStr        ("-1/2");
      static std::string minusOneThirdStr       ("-1/3");
      static std::string minusOneQuarterStr     ("-1/4");
      static std::string minusOneStr            ("-1");
      static std::string minusTwoThirdsStr      ("-2/3");
      static std::string minusThreeQuartersStr  ("-3/4");
      static std::string minusFiveQuartersStr   ("-5/4");
  
      static Field zero               = convert<Field>(zeroStr);
      static Field oneHalf            = convert<Field>(oneHalfStr);
      static Field oneThird           = convert<Field>(oneThirdStr);
      static Field oneQuarter         = convert<Field>(oneQuarterStr);
      static Field one                = convert<Field>(oneStr);
      static Field twoThirds          = convert<Field>(twoThirdsStr);
      static Field threeQuarters      = convert<Field>(threeQuartersStr);
      static Field fiveQuarters       = convert<Field>(fiveQuartersStr);

      static Field minusOneHalf       = convert<Field>(minusOneHalfStr);
      static Field minusOneThird      = convert<Field>(minusOneThirdStr);
      static Field minusOneQuarter    = convert<Field>(minusOneQuarterStr);
      static Field minusOne           = convert<Field>(minusOneStr);
      static Field minusTwoThirds     = convert<Field>(minusTwoThirdsStr);
      static Field minusThreeQuarters = convert<Field>(minusThreeQuartersStr);
      static Field minusFiveQuarters  = convert<Field>(minusFiveQuartersStr);

      if (Algorithms::close(value , one))           return oneStr;
      if (Algorithms::close(value , zero))          return zeroStr;
      if (Algorithms::close(value , oneHalf))       return oneHalfStr;
      if (Algorithms::close(value , oneThird))      return oneThirdStr;
      if (Algorithms::close(value , oneQuarter))    return oneQuarterStr;
      if (Algorithms::close(value , twoThirds))     return twoThirdsStr;
      if (Algorithms::close(value , threeQuarters)) return threeQuartersStr;
      if (Algorithms::close(value , fiveQuarters))  return fiveQuartersStr;

      if (Algorithms::close(value , minusOneHalf))        return minusOneHalfStr;
      if (Algorithms::close(value , minusOneThird))       return minusOneThirdStr;
      if (Algorithms::close(value , minusOneQuarter))     return minusOneQuarterStr;
      if (Algorithms::close(value , minusOne))            return minusOneStr;
      if (Algorithms::close(value , minusTwoThirds))      return minusTwoThirdsStr;
      if (Algorithms::close(value , minusThreeQuarters))  return minusThreeQuartersStr;
      if (Algorithms::close(value , minusFiveQuarters))   return minusFiveQuartersStr;

      std::ostringstream buff;
      buff << value;
      return buff.str();
    }
  };
  /** \ingroup symmetryConcepts
   *
   * \brief A to-be-specialized class for implementing symmetry element names.
   *
   *
   */
  template<size_t DIM>
  class SymmetryElementName: 
    public std::string 
  {};
    
  /** \ingroup symmetryConcepts
   *
   * \brief A class for implementing 2D symmetry element names.
   *
   * The names uniquely identify a symmetry element class.
   */
  template<>
  class SymmetryElementName<2>: 
    public std::string 
  {
  public:
    typedef SymmetryElementName                     ThisType;

    //====================================================================== Constructors
    /**
     * \brief  Sets the name from a std::string
     */
    SymmetryElementName(const std::string& name):
      std::string(name)
    { }

    /**
     * \brief  Sets the name from a char*
     */
    SymmetryElementName(const char*& name):
      std::string(name)
    { }

    /**
     * \brief  The copy constructor 
     */
    SymmetryElementName(const ThisType& name):
      std::string(name)
    { }

    /**
     * \brief  Sets the name from a given Symmetry Element
     */
    template<typename Field, typename Algorithms>
    SymmetryElementName(const SymmetryElement<Field,2,Algorithms>& element):
      std::string("?")
    { 
      std::ostringstream buff;
      setDirectionPart(buff, element);
      buff << element.type;
      setPositionPart(buff, element);
      (*this) = buff.str();
    }

    //======================================================================
    //======================================================================

    /** 
     * \brief Check to make sure the glide ratio is right and return it.
     *
     * \note This is a 2D function, should be seperated out! 
     */
    void setDirComponent(std::ostream& buff, 
			 const int&  netDirComponent,
			 bool isFirst,
			 std::string postfix) const {
      if (netDirComponent == 0) 
	return;
      if (!isFirst)
	buff << " ";
      switch (netDirComponent) {
      case -1: {
	buff << "-";
	break;
      }
      case 1: {
	break;
      }
      default:
	buff << netDirComponent;
      }
      buff << postfix;
    }

    /**
     * \brief Add to the given ostream the direction part of the name
     *        corresponding to the given Symmetry Element.
     *
     */
    template<typename Field, typename Algorithms>
    void setDirectionPart(std::ostream& buff, 
			  const SymmetryElement<Field,2,Algorithms>& element) {
      
      static Field zero(0);

      // Don't print anything for a zero direction (e.g. rotations)
      if (element.netDirection[0] == 0 && 
	  element.netDirection[1] == 0)
	return;
      
      element.checkDirection();
      Field gRatio         = element.glideFraction;
      bool  zeroGlideRatio = Algorithms::close(gRatio,zero);

      if (!zeroGlideRatio)
	buff << Fraction<Field,Algorithms>::GET(gRatio) << " ";

      // The a direction
      setDirComponent(buff,element.netDirection[0],true,"a");

     // The b diirection
      setDirComponent(buff,element.netDirection[1],(element.netDirection[0] == 0),"b");

      if (!zeroGlideRatio)
	buff << " ";
      else
	buff << " ";
    }
      
    template<typename Field, typename Algorithms>
    void setPositionPart(std::ostream& buff, 
			 const SymmetryElement<Field,2,Algorithms>& element) {
      
      Field zero(0);
      // Don't write anything for a zero position
      if (Algorithms::close(element.cellPosition[0], zero) &&
	  Algorithms::close(element.cellPosition[1], zero) )
	return;

      buff << "(";
      buff << Fraction<Field,Algorithms>::GET(element.cellPosition[0]);
      buff << ",";
      buff << Fraction<Field,Algorithms>::GET(element.cellPosition[1]);
      buff << ")";
    }



  };
  
} /* namespace psimag */

#endif   //Psimag_SYMMETRY_ELEMENT_NAME

/*@}*/

