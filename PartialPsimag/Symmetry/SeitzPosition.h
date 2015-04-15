//-*-C++-*-
#ifndef  Psimag_Seitz_Position
#define  Psimag_Seitz_Position

/** \ingroup extendedMatrices */
/*@{*/

/** \file SeitzPosition.h
 *
 * Contains a class for implementing seitz position vectors.
 */

#include "SeitzVector.h"
#include "SeitzTranslation.h"

namespace psimag {
  
  //====================================================================== 
  //============================================================

  /** \ingroup extendedMatrices
   *
   * \brief A class for implementing seitz position vectors
   */
  template<typename Field, size_t DIM> 
  class SeitzPosition: public SeitzVector< Field, DIM, 1 > 
  {
  public:

    typedef          SeitzVector<Field,DIM,1> BaseType;
    typedef typename BaseType::ArrayType      ArrayType;

   /** The Default Constructor produces the zero position. */
    SeitzPosition(): BaseType() {}

    /** Construct a position whose components are set to the given value. */
    SeitzPosition(const Field& val): BaseType(val) {}

    /** Construct a position whose components are set to the given value. */
    //SeitzPosition(const Field& val): BaseType(val) {}

    /** Construct a position whose components are set to the given value array. */
    SeitzPosition(const ArrayType& vals):  BaseType(vals) {}

     /** Construct a position whose components are set to the given value. */
    SeitzPosition(const Vec<Field, DIM>& v):  BaseType(v) {}

    /** Construct a position whose components are set to the given value. */
    SeitzPosition(const BaseType& v): BaseType(v) {}

    /** Construct a position whose components are set to the given value. */
    SeitzPosition(const SeitzPosition<Field, DIM>& v): BaseType(v) {}

  };
  
  
  //====================================================================== 

  template<typename Field>
  SeitzPosition<Field,3> seitzPosition(Field x, Field y, Field z) {
    SeitzPosition<Field,3> result;
    result[0] = x;
    result[1] = y;
    result[2] = z;
    return result;
  }

  template<typename Field>
  SeitzPosition<Field,2> seitzPosition(Field x, Field y) {
    SeitzPosition<Field,2> result;
    result[0] = x;
    result[1] = y;
    return result;
  }


  template<typename Field, typename IN_TYPE> 
  SeitzPosition<Field,2> seitzPosition(IN_TYPE t0, IN_TYPE t1) {
    SeitzPosition<Field,2> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    return result;
  }

  template<typename Field, typename IN_TYPE> 
  SeitzPosition<Field,3> seitzPosition(IN_TYPE t0, IN_TYPE t1, IN_TYPE t2) {
    SeitzPosition<Field,3> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    result[2] = convert<Field>(t2);
    return result;
  }
  
  //======================================================================
  
  /**
   * Return the SeitzTranslation from the difference of two positions.
   */
  template<typename Field, size_t DIM>
  SeitzPosition<Field,DIM> bisector(const SeitzPosition<Field,DIM>& p1, 
				    const SeitzPosition<Field,DIM>& p2) {    
    
    SeitzPosition<Field,DIM>  result( (p2-p1)/Field(2) + p1 );
    
    return result;
  }
  
  /**
   * Return the SeitzTranslation from the difference of two positions.
   */
  template<typename Field, typename Algorithms>
  Field slope(const SeitzPosition<Field,2>&  p1, 
	      const SeitzPosition<Field,2>&  p2) {    
    SeitzTranslation<Field,2> d  = p2 - p1;
    return slope<Field,Algorithms>(d);
  }
  
  /**
   * Return the SeitzTranslation from the difference of two positions.
   */
  template<typename Field, typename Algorithms>
  SeitzPosition<Field,2> xIntercept(const SeitzPosition<Field,2>& position1, 
				    const SeitzPosition<Field,2>& position2) {    
    
    static Field zero(0);
    
    if (Algorithms::close(position1[1], position2[1]))
      throw std::range_error("No xIntercept exists!");
    
    if (Algorithms::close(position1[1], zero))
      return position1;
    
    if (Algorithms::close(position2[1], zero))
      return position2;
    
    if (Algorithms::close(position1[0], position2[0]))
      return seitzPosition<Field,Field>(position1[0],zero);
    
    Field m = slope<Field,Algorithms>(position1, position2);
    return seitzPosition<Field,Field>(position1[0] - position1[1]/m, zero);
  }

  /**
   * Return the SeitzTranslation from the difference of two positions.
   */
  template<typename Field, typename Algorithms>
  SeitzPosition<Field,2> yIntercept(const SeitzPosition<Field,2>& position1, 
				    const SeitzPosition<Field,2>& position2) {    
    
    static Field zero(0);
    
    if (Algorithms::close(position1[0], position2[0]))
      throw std::range_error("No yIntercept exists!");
    
    if (Algorithms::close(position1[0], zero))
      return position1;
    
    if (Algorithms::close(position2[0], zero))
      return position2;
    
    if (Algorithms::close(position1[1], position2[1]))
      return seitzPosition<Field,Field>(zero, position1[1]);

    Field m = slope<Field,Algorithms>(position1, position2);
    return seitzPosition<Field,Field>(zero, position1[1] - m* position1[0]);
  }

  /**
   * Return the SeitzTranslation from the difference of two positions.
   */
  template<typename Field, typename Algorithms>
  SeitzPosition<Field,2> intercept(const SeitzPosition<Field,2>& position1, 
				   const SeitzPosition<Field,2>& position2) {    

    if (Algorithms::close(position1[0], position2[0]) &&
	Algorithms::close(position1[1], position2[1])) {
      std::ostringstream buff;
      buff << " Logic_error in SeitzPosition.h intercept(" << position1 << "," << position2 << ")" << std::endl;
      throw std::logic_error(buff.str());
    }
    
    try { 
      return xIntercept<Field,Algorithms>(position1, position2);
    }
    catch (std::range_error&) {
      return yIntercept<Field,Algorithms>(position1, position2);
    }
  }

  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<typename Field, size_t DIM> 
  Tag toXML(const SeitzPosition<Field,DIM> pos, 
	    std::string name="SeitzPosition") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << pos[i] ;
    return result;
  }

  
}  /* namspace psimag */

#endif //  Psimag_Seitz_Position
/*@}*/
