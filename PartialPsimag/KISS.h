//-*-C++-*-

#ifndef PSIMAG_KISS_H
#define PSIMAG_KISS_H

#include "Real.h"

namespace psimag {

/**  \class KISS
 *   \brief Keep It Simple Stupid random number generator.
 *
 *   Originally by George Marsaglia and Arif Zaman of Florida State University.  
 *   This random number generator was designed specifically to have extremely 
 *   long periods (2^96) and have good qualities for low as well as high bits.
 *   Based on 32 bit integers, so portability is an issue.
 *
 *   \ingroup RNG
 *   \Author Gregory Brown and Hwee Kuan Lee
 */
class KISS {

  public:

    typedef Real   result_type;
    // Because of bit operatirions, seed_type must be 32-bit integer
    // For some systems, need to use a regular int
#ifdef KISS_INT
    typedef unsigned int seed_type; 
#else
    typedef unsigned long int seed_type;  
#endif
    static const bool has_fixed_range = false;
    static result_type min() { return 0; }
    static result_type max() { return 1; }

    explicit KISS(seed_type s=1) : _small(217), _fixed(20030711)
    { 
      seed(s);
    }
   
    /// Seed the random number generator  
    void seed(seed_type s) 
    { 
      _setState(_small,s,_fixed);
    }

    /** Choose a family of random number generators.
     *  The second argument can be any unsigned 32-bit
     *  integer, but the first must be less than 2580
     *  to guarantee a long-period psuedorandom sequence.
     *
     *  These values will not affect the state of the
     *  random number generator until the next time
     *  seed(s) is called.
     *
     *  In some sense, KISS really takes three seeds.
     *  To stay consistent with the idea that a RNG is
     *  seeded by a single integer, we've split the
     *  seeds into seed(s) and the two specified here.
     */ 
    void family(seed_type a, seed_type b)
    { 
      if(a>2580)
      {
	std::ostringstream msg;
        msg << "__FILE__::__LINE__:"
	    << " psimag::KISS::setFamily(a,b)" << std::endl
	    << " a>2580 does not guarantee a long-period sequence" << std::endl
	    << " specified value is a=" << a << std::endl <<std::flush;
	throw std::logic_error(msg.str());
      }
      _small = a;
      _fixed = b;
    }

    /// Return random number uniformly distributed on [0,1) 
    result_type operator ()()
    {
      // This code is taken directly from the reprint
      j =   j ^ (j<<17);
      k = ( k ^ (k<<18) ) & 0x7FFFFFFF;
      i = 69069*i + 23606797 + (j^=(j>>15)) + (k^=(k>>13));
      // This is changing it to a floating point number
      return static_cast<Real>(i)/static_cast<Real>(4294967295UL);
    }

    /// return random number on [0,n)    
    template<class T>
    T operator () (T n)
    {
      return static_cast<T>( (*this)() * n);
    }

    /// Returns true if states are identical
    friend bool operator == (const KISS& a, const KISS& b)
    { 
      return (a.i == b.i) && (a.j == b.j) && (a.k == b.k); 
    }

  private:

    void _setState(seed_type a, seed_type b, seed_type c)
    { i=a; j=b; k=c; }

  private:

    seed_type i,j,k;  // these constitute the state of the generator

    seed_type _small; // to guarantee long period 1 < _smallSeed < 2580
    seed_type _fixed; // this seed can have any value

  };
 
 
} /* namespace psimag */


#endif /* PSIMAG_KISS_H */
