//-*-C++-*-
#ifndef  Psimag_Seitz_Translation
#define  Psimag_Seitz_Translation

/** \ingroup extendedMatrices */
/*@{*/

/** \file SeitzTranslation.h
 *
 * Contains classes for implementing seitz translation vectors.
 */

#include "SeitzVector.h"
#include "Mat.h"    // For Multiply

namespace psimag {
  
  //============================================================

  /** \ingroup extendedMatrices
   *
   * \brief A class for implementing seitz translation vectors
   */
  template<typename Field, size_t DIM> 
  class SeitzTranslation: public SeitzVector< Field, DIM, 0 > 
  {    
  public:
    /** The Default Constructor produces the zero translation.*/
    SeitzTranslation(): SeitzVector< Field, DIM > ()  {}
    
    /** Construct a translation whose components are set to the given value. */
    template<typename IN_TYPE> 
    SeitzTranslation(IN_TYPE val): SeitzVector<Field, DIM>(val) {}

    /** Construct a translation whose components are set to the given value array. */
    template<typename IN_TYPE>
    SeitzTranslation(IN_TYPE vals[DIM]): SeitzVector< Field, DIM > (vals) {}
  
//     /** 
//      * Construct a translation whose components are set to the given value. 
//      *
//      * Convert data types as necessary.
//      */
//     template<typename IN_TYPE> 
//     SeitzTranslation(SeitzTranslation<IN_TYPE,DIM> val):       
//       SeitzVector< Field, DIM > ((SeitzVector< IN_TYPE, DIM >)val)  
//     {}

    /** 
     * Return the length of the translation given the appropriate Metric. 
     */
    Field length(const Mat< Field, DIM, DIM  >& metric) {
      SeitzTranslation temp;
      Multiply(metric, (*this), temp);
      return (*this) * temp;
    }

  };
  
  /**
   * Return the slope of a translation
   */
  template<typename Field, typename Algorithms>
  Field slope(const SeitzTranslation<Field,2>&  t) {
    static Field zero(0);
    if (Algorithms::close(t[0],zero))
      throw std::range_error("vertical slope");
    return t[1]/t[0];
  }
  

  //====================================================================== 

  template<typename Field>
  SeitzTranslation<Field,3> seitzTranslation(Field x, Field y, Field z) {
    SeitzTranslation<Field,3> result;
    result[0] = x;
    result[1] = y;
    result[2] = z;
    return result;
  }

  template<typename Field>
  SeitzTranslation<Field,2> seitzTranslation(Field x, Field y) {
    SeitzTranslation<Field,2> result;
    result[0] = x;
    result[1] = y;
    return result;
  }


  template<typename Field, typename IN_TYPE> 
  SeitzTranslation<Field,2> seitzTranslation(IN_TYPE t0, IN_TYPE t1) {
    SeitzTranslation<Field,2> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    return result;
  }

  template<typename Field, typename IN_TYPE> 
  SeitzTranslation<Field,3> seitzTranslation(IN_TYPE t0, IN_TYPE t1, IN_TYPE t2) {
    SeitzTranslation<Field,3> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    result[2] = convert<Field>(t2);
    return result;
  }
  

  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for SeitzTranslations.
   */
  template<typename Field, size_t DIM> 
  Tag toXML(const SeitzTranslation<Field,DIM> pos, 
	    std::string name="SeitzTranslation") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << pos[i] ;
    return result;
  }

}  /* namspace psimag */

#endif // Psimag_Seitz_Translation
/*@}*/
