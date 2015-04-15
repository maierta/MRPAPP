//-*-C++-*-

#ifndef PSIMAG_AbstractRat_H
#define PSIMAG_AbstractRat_H

///
///\file AbstractRat.h 
///
///\brief Contains a class to represent the output of the FieldParser.
///
///       The parser input strings can be of the form:
///
///       "2/7 + 2/5*sqrt(3) + 5/11*sqrt(5) + . . ."
///
///       The parsed result of the above example would be:
///
///       {1:[2,7], 3:[2,5], 5:[5,11]}  in JSON
///
///       Its c++ type is std::map< long, std::vector<long> >
///        
///\author Mike Summers.           
///                                                           

#include <cstddef>
#include <limits>
#include <cassert>
#include <iostream>
#include <sstream>

#include <string>
#include <map>
#include <vector>
#include "Vec.h"
#include "PSIMAGAssert.h"

namespace psimag {
  
  /* *
   * \brief A class to be used as a common representation of a number
   *        by FieldParser, rational and sqrtExtendedRational
   *
   *
   * \note The class maps the prime number (1, 2, 3, 5 . . ) to a long
   *       vector of size 2, representing the two ints of a rational.
   */

  template<typename Rat>
  Rat rsqrt(int sq) {
    return Rat(sqrt(sq));
  }

  // Testing shows problems, some may be with Rat::IntegerType
  // overflow in gcd.
  template<typename Rat>
  Rat rsqrt_2(int sq) {
    
    assert(sq >= 0);
    if (sq == 0 ) return 0;
    
    typedef typename Rat::IntegerType Int;
  
    Int numerator(1), denominator(1), s(sq);

    size_t count = 0;
    while(true) {
      Int n   = numerator + denominator * s;
      Int d   = numerator + denominator;
      if (n <= 0 || d <= 0 ) 
	break;
      numerator = n;
      denominator = d;
      if (count > 100) break;
    }
    return Rat(numerator, denominator);
  }


  template<typename IntType> class AbstractRat: 
  public std::vector< std::pair<IntType, IntType> >
  {
  private:
    static const int    primes[]; //= {2,3,5,7,11,13,17,19};
    static const size_t n_roots;  //= 8
  public:

    AbstractRat(): std::vector< std::pair<IntType, IntType> >
		   (1, 
		    std::pair<IntType, IntType>(IntType(0), IntType(1)))  {}

    int getPrime(const size_t idx) const {
      int b=1;
      for(size_t i=0; i < n_roots; i++)
	if(idx & (1<<i)) b*=primes[i];
      
      return b;
    }

    size_t getPrimeIndex(const int prime) const {
      if (prime == 1) return 0;
      for(size_t i=0; i<n_roots; i++)
	if (prime == primes[i])
	  return i+1;
      throw std::out_of_range("sqrtExtendedRational.getPrimeIndex");
    }

    void setTerm(const int prime, const std::pair<IntType, IntType>& intPair) {
      size_t index = getPrimeIndex(prime);
      if (index >= this->size())
	this->resize(index+1, std::pair<IntType, IntType>(IntType(0),IntType(1)));
      (*this)[index] = intPair;
    }

    double getDouble() {
      
      double result(0.0);
      
      for (size_t i=0; i < (*this).size(); i++) {

	double term(0.0);

	IntType numerator    = (*this)[i].first;
	IntType denominator  = (*this)[i].second;

	term = (double) numerator / (double) denominator;
	
	if (i > 0) 
	  term *= sqrt(static_cast<double>(getPrime(i)));

	result += term;
      }

      return result;
    }
	float getFloat() { return float(getDouble()); }
    
    template<typename Rat>
    Rat& getRational(Rat& result) {

      typedef typename Rat::IntegerType IntegerType;

      result = Rat(IntegerType(0));

      for (size_t i=0; i < this->size(); i++) {
	
	Rat r( (*this)[i].first, 
	       (*this)[i].second );
	if (i > 0 ) 
	  r *= rsqrt<Rat>(this->getPrime(i));
	
	result += r;
      }
      return result;
    }

  };
  
  template<typename IntType>
  const int AbstractRat<IntType>::primes[]={2,3,5,7,11,13,17,19};

  template<typename IntType>
  const size_t AbstractRat<IntType>::n_roots=8;


  template<typename IntType>
  std::ostream& operator<<(std::ostream&               os, 
			   const AbstractRat<IntType>& abstractRat) {
    
    os << "{ ";
    for (size_t i=0; i< abstractRat.size(); i++) {
      os << "  " << abstractRat.getPrime(i) 
	 << ": [ "
         << abstractRat[i].first << "," 
	 << abstractRat[i].second
	 << " ], ";
    }
    os << "}" << std::endl;
    return os;
  }
} /* namespace psimag */

#endif 
