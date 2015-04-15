//-*-C++-*-

#ifndef PSIMAG_FieldConvert_H
#define PSIMAG_FieldConvert_H

///
///\file FieldConvert.h 
///
///\brief Contains the generic template functions and select
///       specializations for the conversion of scalar 'Field'
///       types. 
///
///       These templates includes a capability to convert char*
///       values to a given Field value. Since Field may be something
///       like a Rational or extended rational object we need to be
///       able to define it's input with strings like "3/2", or
///       perhaps "3/sqrt(2)". At the same time Field could also be
///       they type double so that these same string constants need to
///       be converted to double. 
///
///       As far as I know (MSS) there is no way in c++ to define a
///       meaning for double("1/2*sgrt(3)"). As a result we must resort to a
///       template like, convert<double>("1/2*sgrt(3)").
///        
///\author Mike Summers.           
///                                                           

#include <cstddef>
#include <limits>

#include <iostream>
#include <sstream>

#include <string>
#include <map>
#include <vector>
#include "FieldParser.h"
#include "PSIMAGAssert.h"

namespace psimag {

  template<typename Field>
  inline  Field convert(int value) {
    return Field(value);
  }
  
  template<typename Field>
  inline  Field convert(double value) {
    return Field(value);
  }
  
  template<typename Field>
  inline  Field convert(const char* text) {
    AbstractRat<long> map = FieldParser<long>(text).parse();
    return Field(map);
  }
  
  template<typename Field>
  inline  Field convert(const std::string text) {
    AbstractRat<long> map = FieldParser<long>(text).parse();
    return Field(map);
  }
 
template<>
inline float convert<float>(const char* text) {
    
    float val = FieldParser<long>(text).parseToFloat();

    return val;
  }
 
  template<>
  inline  double convert<double>(const char* text) {
    
    double val = FieldParser<long>(text).parseToDouble();

    return val;
  }

  template<>
  inline  double convert<double>(const std::string text) {
    
    double val = FieldParser<long>(text).parseToDouble();

    return val;
  }

  template<>
  inline  float convert<float>(const std::string text) {
    
    float val = FieldParser<long>(text).parseToDouble();

    return val;
  }
  
  template<>
  inline    int convert<int>(const char* text) {
    return (int) FieldParser<long>(text).parseToDouble();
  }
  
} /* namespace psimag */

#endif //PSIMAG_FieldConvert_H
