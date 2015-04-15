// -*- C++ -*-

#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>
#include <sstream>
#include "AbstractRat.h"
#include "FieldParser.h"

template<class Int> class rational
{
 private:

  Int num, den; // numerator, denominator

  void setValue(Int u, Int v) {
    Int t=gcd(u,v); 

    num=u/t; 
    den=v/t; 
    if(den<0) {
      den=-den;
      num=-num;
    } 
  }

 public:

  typedef Int            IntegerType;
  typedef rational<Int>  ThisType;

  rational<Int>(Int u, Int v) {setValue(u,v);}
  rational<Int>(Int u) {num=u; den=Int(1);}
  rational<Int>() {num=Int(0); den=Int(1);}
  rational<Int>(const rational<Int> &a){num=a.num; den=a.den;}

  // mss
  rational<Int>(psimag::AbstractRat<long>& abstractRat)
  {
    abstractRat.getRational(*this);
  }

  operator psimag::AbstractRat<Int>() {

    psimag::AbstractRat<Int> result(num, den);
    return result;
  }

  // mss
  explicit
  rational<Int>(const double val) {

    std::ostringstream buffer;
    buffer << val;
    std::pair<Int, Int> numAndDen = psimag::numStringToFraction<Int>(buffer.str());

    setValue(numAndDen.first, numAndDen.second);
  }

  rational<Int> operator*(rational<Int> v)
    {
      rational<Int> w; Int d1, d2;

      d1=gcd(num,v.den);
      d2=gcd(den,v.num);
      w.num=(num/d1)*(v.num/d2);
      w.den=(den/d2)*(v.den/d1);
      /*
      w.num=num*v.num;
      w.den=den*v.den;
      d=gcd(w.num,w.den);
      w.num/=d; w.den/=d;
      */
      return w;
    }
  rational<Int> operator*=(rational<Int> v)
    {
      Int d1, d2;

      d1=gcd(num,v.den);
      d2=gcd(den,v.num);
      num=(num/d1)*(v.num/d2);
      den=(den/d2)*(v.den/d1);
      
      return *this;
    }
  rational<Int> operator*(Int a)
    {
      rational<Int> w; Int d;

      w.den=den;
      d=gcd(a,den);
      if(d==Int(1))
        w.num=num*a;
      else
        {
          w.num=num*(a/d);
          w.den/=d;
        }
      return w;
    }
  rational<Int> operator*=(Int a)
    {
      Int d;

      d=gcd(a,den);
      if(d==Int(1))
        num*=a;
      else
        {
          num*=a/d;
          den/=d;
        }
      return *this;
    }

  rational<Int> operator/(rational<Int> v)
    {
      rational<Int> w; Int d1,d2;

      d1=gcd(num,v.num);
      d2=gcd(den,v.den);
      w.num=(num/d1)*(v.den/d2);
      w.den=(den/d2)*(v.num/d1);
      /*
      w.num=num*v.den;
      w.den=den*v.num;
      d=gcd(w.num,w.den);
      w.num/=d; w.den/=d;
      */
      if(w.den<0){w.den=-w.den;w.num=-w.num;}
      return w;
    }
  rational<Int> operator/=(rational<Int> v)
    {
      Int d1,d2;

      d1=gcd(num,v.num);
      d2=gcd(den,v.den);
      num=(num/d1)*(v.den/d2);
      den=(den/d2)*(v.num/d1);

      if(den<0){den=-den;num=-num;}
      return *this;
    }

  rational<Int> operator+(rational<Int> v)
    {
      rational<Int> w; Int d1,d2,t;

      d1=gcd(den,v.den);
      if(d1==Int(1))
        {
          w.num=num*v.den+den*v.num;
          w.den=den*v.den;
        }
      else
        {
          t=num*(v.den/d1)+v.num*(den/d1);
          d2=gcd(t,d1);
          w.num=t/d2;
          w.den=(den/d1)*(v.den/d2);
        }
      /*
      w.num=num*v.den+den*v.num;
      w.den=den*v.den;
      d=gcd(w.num,w.den);
      w.num/=d; w.den/=d;
      */
      if(w.den<0){w.den=-w.den;w.num=-w.num;}
      return w;
    }
  rational<Int> operator+=(rational<Int> v)
    {
      Int d1,d2,t;

      d1=gcd(den,v.den);
      if(d1==Int(1))
        {
          num=num*v.den+den*v.num;
          den=den*v.den;
        }
      else
        {
          t=num*(v.den/d1)+v.num*(den/d1);
          d2=gcd(t,d1);
          num=t/d2;
          den=(den/d1)*(v.den/d2);
        }

      if(den<0){den=-den;num=-num;}
      return *this;
    }

  rational<Int> operator-(rational<Int> v)
    {
      rational<Int> w; Int d1,d2,t;

      d1=gcd(den,v.den);
      if(d1==Int(1))
        {
          w.num=num*v.den-den*v.num;
          w.den=den*v.den;
        }
      else
        {
          t=num*(v.den/d1)-v.num*(den/d1);
          d2=gcd(t,d1);
          w.num=t/d2;
          w.den=(den/d1)*(v.den/d2);
        }
      /*
      w.num=num*v.den-den*v.num;
      w.den=den*v.den;
      d=gcd(w.num,w.den);
      w.num/=d; w.den/=d;
      */
      if(w.den<0){w.den=-w.den;w.num=-w.num;}
      return w;
    }

  rational<Int> operator-=(rational<Int> v)
    {
      Int d1,d2,t;

      d1=gcd(den,v.den);
      if(d1==Int(1))
        {
          num=num*v.den-den*v.num;
          den=den*v.den;
        }
      else
        {
          t=num*(v.den/d1)-v.num*(den/d1);
          d2=gcd(t,d1);
          num=t/d2;
          den=(den/d1)*(v.den/d2);
        }

      if(den<0){den=-den;num=-num;}
      return *this;
    }

  rational<Int> abs() { rational<Int> w; w.num=abs(num); w.den=abs(den); return w; }

  bool operator==(rational<Int> v) { return (num==v.num) && (den==v.den); }
  bool operator!=(rational<Int> v) { return (num!=v.num) || (den!=v.den); }
  bool operator>(rational<Int> v) { return (num*v.den) > (den*v.num); }
  bool operator<(rational<Int> v) { return (num*v.den) < (den*v.num); }

  Int numerator(){return num;}
  Int denominator(){return den;}
  };

template<class Int>
rational<Int> abs(rational<Int> a)
  {return rational<Int>(abs(a.numerator()),abs(a.denominator()));}

template<class Int>
std::ostream& operator<<(std::ostream& s, rational<Int> a)
{
  return s<<a.numerator()<<"/"<<a.denominator();
}

// Find the greatest common divisor using Euclid's algorithm ( Knuth 4.5.2, algorithm A)
template<class Int> Int gcd(Int a, Int b)
{
  Int r,u,v;

  u=abs(a); v=abs(b);
  while(v!=Int(0))
    {
      r=u%v;
      u=v;
      v=r;
    }
  return(u);
}

#endif
