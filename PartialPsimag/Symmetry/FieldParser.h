//-*-C++-*-

#ifndef PSIMAG_FieldParser_H
#define PSIMAG_FieldParser_H

///
///\file FieldParser.h 
///
///\brief Contain a templated parser for reading field data from strings.
///
///       The input strings can be of the form:
///
///       "2/7 + 2/5*sqrt(3) + 5/11*sqrt(5) + . . ."
///
///       The template parameter, IntType, if the type used by the parser
///       to represent integers.
///
///       The parsed result of the above example would be:
///
///       {1:[2,7], 3:[2,5], 5:[5,11]}  in JSON
///
///       Its c++ type is std::map< IntType, std::vector<IntType> >
///        
///\author Mike Summers.           
///                                                           

#include <cstddef>
#include <limits>
#include <math.h>
#include <stdexcept>

#include <iostream>
#include <sstream>

#include <string>
#include <map>
#include <vector>
#include "AbstractRat.h"
#include "PSIMAGAssert.h"

namespace psimag {
  
  // This belongs in a general string utility package
inline  std::vector<std::string> split(const std::string input, const char sep) {
    
    std::vector<std::string> result;
    size_t start  = 0;
    size_t len    = input.length();
    while(true) {
      size_t sepPos = input.find(sep, start);
      if (sepPos == std::string::npos) {
	result.push_back(input.substr(start, len-start));
	return result;
      }
      else {
	result.push_back(input.substr(start, sepPos-start));
	start = sepPos + 1;
      }
    }
  }
    
  template<typename IntType> 
inline  IntType getInt(const std::string& intString) {
    if (intString.length() == 0 ||
	intString.find_first_not_of(" ") == std::string::npos) 
      return IntType(1);
    IntType result;
    std::istringstream reader(intString, std::istringstream::in);
    reader >> std::ws;
    reader >> result;
    return result;
  }
  
  template<typename IntType> 
std::pair<IntType, IntType> numStringToFraction(std::string decimalString) {
    
    std::pair<IntType, IntType> result(IntType(1), IntType(1));

    size_t len            = decimalString.length();
    size_t decimalPos     = decimalString.find('.');

    if(std::string::npos == decimalPos) {
      // Just an integer
      result.first = psimag::getInt<IntType>(decimalString);
      return result;
    }

    size_t places         = len - (decimalPos+1);
    std::string numString = decimalString.erase(decimalPos,1);
    result.first = psimag::getInt<IntType>(numString);
    
    std::string denString(places+1,'0');
    denString[0] = '1';
    result.second = psimag::getInt<IntType>(denString);
    
    return result;
  }

  template<typename IntType> 
  class FieldParser {
    
  private:
    static std::string   sqrtString;
    AbstractRat<IntType> abstractRat;
    
  public:
    
    std::string  input;

    FieldParser(const std::string value): input(value) { parseInput(); }
    FieldParser(const char* value):       input(value) { parseInput(); }
    
    void parseInput() {
      std::vector<std::string> terms = split( input, '+');
      for (size_t i=0; i<terms.size(); i++)
	getTerm(terms[i]);
    }
    
    const AbstractRat<IntType> parse() const {
      return abstractRat;
    }
   
	float parseToFloat() { return abstractRat.getFloat(); }
 
    double parseToDouble() {
      return abstractRat.getDouble();
    }

  protected:
    
    void getTerm(const std::string& termStr){
      
      std::vector<std::string> factors = split( termStr, '*');
      int prime(1);
      std::pair<IntType, IntType> intPair(IntType(1), IntType(1));

      for (size_t i=0; i<factors.size(); i++) {
	
	int aPrime = getPrime(factors[i]);
	if (aPrime != -911) { 
	  prime = aPrime;
	  continue;
	}
	// If not a prime factor then a intPair factor
	intPair = getIntPair(factors[i]);
      }
      abstractRat.setTerm(prime, intPair);
    }
    
    std::pair<IntType, IntType> getIntPair(const std::string& inString) {
      
      std::vector<std::string> numStrs = split(inString, '/');
      std::pair<IntType, IntType>  result(IntType(1), IntType(1));
      
      if (numStrs.size() > 2)
	throw std::range_error("sqrtExtendedRational.getIntPair() to many '/'");
      
      if (numStrs.size() == 0)
	return result;
      
      std::pair<IntType, IntType> numPair = numStringToFraction<IntType>(numStrs[0]);
      if (numStrs.size() == 1) {
	return numPair;
      }
      
      std::pair<IntType, IntType> denomPair = numStringToFraction<IntType>(numStrs[1]);
      
      result.first  = numPair.first  * denomPair.second;
      result.second = numPair.second * denomPair.first;
      
      return result;
    }
      
    IntType getPrime(std::string& term){

      size_t sqrtPos = term.find(sqrtString);
      if (sqrtPos == std::string::npos)
	return IntType(-911);

      size_t endPos = term.find(')');
      if (endPos == std::string::npos)
	throw std::ios::failure("No ending ')' provided in sqrtExtendedRational.getPrime");
      
      size_t      startPos    = sqrtPos + sqrtString.length();
      std::string primeString = term.substr(startPos, endPos-startPos);
      
      return getInt<IntType>(primeString);
    }

  };

  template<typename IntType>
  std::string FieldParser<IntType>::sqrtString("sqrt(");


  template<typename IntType>
  std::ostream& operator<<(std::ostream& os, const FieldParser<IntType>& fp)
  {
    os << "FieldParser('" << fp.input <<"')" << std::endl;
    os << "  abstractRat = (" << fp.parse() <<")" << std::endl;
    return os;
  }
} /* namespace psimag */

#endif //PSIMAG_FieldParser_H
