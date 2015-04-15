//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file BasicCrystalAlgorithms.h  
 *
 *  Contains a the CrystalBase class.
 */

#ifndef Psimag_BasicCrystalAlgorithms_H
#define Psimag_BasicCrystalAlgorithms_H

#include <iostream>                                  
#include <sstream>                                  
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "SpaceGroupConstructor.h"
#include "Occupant.h"

namespace psimag {
  
//======================================================================
 
class BasicCrystalAlgorithms {
public:

  typedef Simple2DReducer<double,Occupant,BasicCrystalAlgorithms>  Reducer;       // Occupant should be a template arg? &*&*&*&*
  typedef SpaceGroupConstructor<double,2>                          Classifier;
   
  static double threshold() {
    static const double result(0.00001);
    return result;
  }

  template<typename Field>
  static bool close(const Field v1, const Field v2) { 
    Field diff = v1-v2;
    if (diff < Field(0))
      diff = diff * Field(-1);

    if (diff < BasicCrystalAlgorithms::threshold())
      return true;
    return false;
  }

  template<typename Field>
  static int sign(const Field v) { 
    static Field zero(0);
    if (close(v,zero)) return 0;
    if (v < zero) return -1;
    return 1;
  }

  template<typename Field>
  static Field normalize(const Field v1) { 
    static Field zero(0);
    static Field one(1);
    static Field minusone(-1);
    if (close(v1, zero))     return zero;
    if (close(v1, one))      return one;
    if (close(v1, minusone)) return minusone;
    return v1;
  }

  template<typename Field, size_t DIM, typename VectorType>
  static void normalizeVector(VectorType& v) { 
    static Field one(1);
    for (size_t i=0; i< DIM; i++)
      v[i] = modulus(v[i],one);
    return v;
  }

  template<typename Field>
  static Field modulus(const Field v1, const Field v2) {  

    static const Field zero(0);
    static const Field one(1);
    static const Field minusOne(-1);

    if (close(v1,zero))     return zero;
    if (close(v1,one))      return zero;
    if (close(v1,minusOne)) return zero;

    Field result = v1 - v2 * floor(v1/v2);

    if (close(result,zero))     return zero;
    if (close(result,one))      return zero;
    if (close(result,minusOne)) return zero;

    if (close(result,one)) 
     throw std::logic_error("BasicAlgorithms::modulus error!");
    return result;
  } 

}; 

} /** namespace psimag */

#endif
/*@}*/
