//-*-C++-*-
#ifndef  Psimag_LatticeCoordinates
#define  Psimag_LatticeCoordinates

/** \ingroup extendedMatrices */
/*@{*/

/** \file LatticeCoordinates.h
 *
 * Contains a class for implementing LatticeCoordinates (SeitzTranslation subclass) objects.
 */

#include "SeitzVector.h"
#include "CellTranslation.h"
#include "CellPosition.h"

namespace psimag {

  template<typename Field, size_t DIM> class CellTranslation;
  template<typename Field, size_t DIM, typename Algorithms> class CellPosition;
  template<typename Field, size_t DIM, int IND, typename Algorithms> bool isNegative(const SeitzVector<Field,DIM,IND>& v);

  /** \ingroup extendedMatrices
   *
   * \brief A marker class indicating that the Translation indicates a direction from within a Cell.
   *
   * These can be converted to  and from MillerIndices for a given Cell.
   */
  template<size_t DIM> 
  class LatticeCoordinates: public SeitzVector< int, DIM, 0 > {

  public:

    enum{ IND=0 };
    typedef LatticeCoordinates<DIM>   ThisType;
    typedef SeitzVector<int,DIM,IND>  BaseType;

    /** The Default Constructor produces the zero translation. */
    LatticeCoordinates(): BaseType() {}

    /** Construct a translation whose components are set to the given value. */
    LatticeCoordinates(int val): BaseType(val) {}

    /** Copy Construct a translation. */
    LatticeCoordinates(const ThisType& v): BaseType(v) {}

    /** Copy Construct a translation. */
    LatticeCoordinates(const BaseType& v): BaseType(v) {}

    template<typename Field>
    CellTranslation<Field,DIM> cellTranslation() const {
      CellTranslation<Field,DIM> result;
      for (size_t i=0; i< DIM; i++)
	result[i] = (*this)[i];
      return result;
    }

    template<typename Field, typename Algorithms>
    CellPosition<Field,DIM,Algorithms> cellPosition() const {
      CellPosition<Field,DIM,Algorithms> result;
      for (size_t i=0; i< DIM; i++)
	result[i] = (*this)[i];
      return result;
    }


    template<typename Algorithms>
    LatticeCoordinates<DIM> normalize() const {
      
      if(isNegative<int,DIM,0,Algorithms>(*this)) {
	LatticeCoordinates<DIM> result(*this);
	result *= -1; //make it positive
	return result;
      }
      return *this;
    }    

//       if(result[0] > 0)
// 	return result;

//       if (result[0] == 0)
// 	if (result[1] >= 0)
// 	  return result;

    ThisType& operator = (const ThisType& v)
    { 
      const BaseType& bv = v;
      BaseType& bt = (*this);
      bt = bv;
      return *this; 
    }


  };

  /**
   * The less operator.
   */
  template<size_t DIM>
  bool operator< (const LatticeCoordinates<DIM>& lhs, 
		  const LatticeCoordinates<DIM>& rhs) {
    
    for (size_t i = 0; i<DIM; i++) {
      int d1 = lhs[i] - rhs[i];
      if (d1 == 0) continue;
      return d1 < 0.0;
    }
    
    return false;
  }

  template<typename Field>
  Field getSlope(const LatticeCoordinates<2> coord) {
    return convert<Field>(coord[1])/convert<Field>(coord[0]);
  }

/* ATTENTION/FIXME: Header files should not contain implementations
        (because it causes multiple definitions when using more than one compilation unit)
        quick fix: I added the keyword inline (Gonzalo) */

  inline bool isVertical(LatticeCoordinates<2> coord) {
    return (coord[0] == 0 && coord[1] != 0);
  }
/* ATTENTION/FIXME: Header files should not contain implementations
        (because it causes multiple definitions when using more than one compilation unit)
        quick fix: I added the keyword inline (Gonzalo) */

  inline LatticeCoordinates<2> latticeCoordinate(int t0, int t1) {
    LatticeCoordinates<2> result;
    result[0] = t0;
    result[1] = t1;
    return result;
  }

  template<typename Field, typename Algorithms>
  LatticeCoordinates<2> latticeCoordinate(Field slope) {
    static Field zero(0);
    static int   netMax(20);

    if (Algorithms::close(zero,slope))
      return latticeCoordinate(1,0);

    for (int m=1; m <netMax; m++) {
      for (int n=-netMax; n <netMax; n++) {
	Field ratio = Field(n)/Field(m);
	if (Algorithms::close(ratio,slope))
	  return latticeCoordinate(m,n);
	if (ratio > slope)
	  break;
      }
    }
    std::ostringstream buff;
    buff << "latticeCoordinate(" << slope << ") failed!";
    throw std::logic_error(buff.str());
  }

  //====================================================================== 

  template<typename Field, typename Algorithms>
  LatticeCoordinates<2> latticeCoordinate(const CellTranslation<Field,2>& trans) {

    static Field zero(0);
    if (Algorithms::close(trans[0],zero)) {
      if(Algorithms::close(trans[1],zero)) {
	return latticeCoordinate(0,0);
      } else {
	if (trans[1] > 0) {
	  return latticeCoordinate(0,1);
	} else {
	  return latticeCoordinate(0,-1);
 	}
     }
   }
    //slopes only work in cartesian?
    LatticeCoordinates<2> result = latticeCoordinate<Field,Algorithms>(slope<Field,Algorithms>(trans));

    // Restore the direction if need be.
    if( (isNegative<Field,2,0,Algorithms>(trans)) && 
	!(isNegative<int,2,0,Algorithms>(result)) )
      result *= -1;

    return result;
  }

  //====================================================================== 
  /*
   *
   * \brief Throws a range_error if trans is less than its lattice
   *        coordinate, otherwise it returns the maximal
   *        latticeCoordinate contained in trans.
   *
   * \note If trans is bad in some way it returns a logic_error.
   */
  template<typename Field, typename Algorithms>
  LatticeCoordinates<2> maxContainedLatticeCoordinate(const CellTranslation<Field,2>& trans) {
    
    // pick which index to work with 0 or 1
    int                   testIndex = 0;
    LatticeCoordinates<2> minResult = latticeCoordinate<Field,Algorithms>(trans);
    if (minResult[0] == 0) {
      if (minResult[1] == 0)
	throw std::range_error("maxContainedLatticeCoordinate given a zero translation!");
      testIndex = 1;
    }
    
    // If necessary flip the vectors so that we are comparing with positive numbers.
    int flip = 1;
    if (minResult[testIndex] < 0) 
      flip = -1;
    Field flipField = Field(flip);
      
    // Try different factors until the result is to big, 
    // then return the previous one (if i> 1)
    int i = 1;
    for (; i< 5; i++) {
      LatticeCoordinates<2> result(minResult * i);
      if(Algorithms::close(trans[testIndex], Field(result[testIndex])))
	return result;
      if (Field(result[testIndex] * flip) > (trans[testIndex] * flipField)) { // then result is to big
	if (i==1) {
	  throw std::range_error("The given translation does not contain it's latticeCoordinate.");
	} else {
	  return minResult * (i-1);
        }
      }
    }
    throw std::logic_error("maxContainedLatticeCoordinate failed!");
  }


//     static Field zero(0);
//     static Field one(1);
//     static Field two(2);
//     static Field three(3);
//     static Field minusOne(-1);
//     static Field minusTwo(-2);
//     static Field minusThree(-3);

//     Field theSlope = latticeCoordinate(slope(trans))
    
//     for (int l=1; l<3; l++) {
      
//       LatticeCoordinates<DIM> result;
//       bool                    isIntegral = true;
      
//       for (size_t i=0; i<DIM; i++) {
	
// 	result[i] = -911;
// 	Field value = trans[i] * Field(l);

// 	if (Algorithms::close(value,zero))       result[i] =  0;
// 	if (Algorithms::close(value,one))        result[i] =  1;
// 	if (Algorithms::close(value,two))        result[i] =  2;
// 	if (Algorithms::close(value,three))      result[i] =  3;
// 	if (Algorithms::close(value,minusOne))   result[i] = -1;
// 	if (Algorithms::close(value,minusTwo))   result[i] = -2;
// 	if (Algorithms::close(value,minusThree)) result[i] = -3;
// 	if (result[i] == -911) { 
// 	  isIntegral = false; 
// 	  break;
// 	} 
//       }
//       if (isIntegral)
// 	return result;
//     }
//     std::ostringstream buff;
//     buff << "latticeCoordinate(" << trans << ") failed!";
//     throw std::logic_error(buff.str());
//   }
  
  //====================================================================== 
/* ATTENTION/FIXME: Header files should not contain implementations
        (because it causes multiple definitions when using more than one compilation unit)
        quick fix: I added the keyword inline. (Gonzalo) */

  inline LatticeCoordinates<3> latticeCoordinate(int t0, int t1, int t2) {
    LatticeCoordinates<3> result;
    result[0] = t0;
    result[1] = t1;
    result[2] = t2;
    return result;
  }

  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<size_t DIM> 
  Tag toXML(const LatticeCoordinates<DIM> lTrans,
	    std::string name="LatticeCoordinates") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << lTrans[i] ;
    return result;
  }


} /* namspace psimag */

#endif 
/*@}*/
