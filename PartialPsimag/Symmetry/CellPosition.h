//-*-C++-*-

#ifndef PSIMAG_CellPosition_H
#define PSIMAG_CellPosition_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file  CellPosition.h
 *  Contains the class definition for objects representing positions in cells.
 */
 
#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Real.h"
#include "Tag.h"

#include "PSIMAGAssert.h"
#include "SeitzVectors.h"
#include "SeitzPosition.h"
#include "CellTranslation.h"

namespace psimag {

  /** \ingroup latticeConcepts
   *
   *  \brief For the representation of positions within cells.
   *
   *  Positions are given in the coordinate system defined by the
   *  cells translation vectors. This is sometimes reffered to as the
   *  fractional coordinates.
   *
   *  The constructors are coded so that all cell positions (which are
   *  SeitzPosition<Field,DIM>s) will have components in the interval
   *  [0,1]. If a cell constructor is given a component value outside
   *  of [0,1] the translationally equivalent cell position,
   *  identified by the snap member function will be used instead.
   *
   *  The cell position's operator< has a tolerence value built into
   *  it so that cell positions appear to be the same position if they are
   *  close enough according to the given tolerence value.
   *
   *  Setting up the operator< in this way allows pattern objects to
   *  match approximately. This is made use of by the pattern's
   *  satifies(symmetryOperation) member function.
   *
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM: Dimensionality of the lattice being
   *             represented. Legal values are 1,2,3.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  class CellPosition: public SeitzPosition<Field,DIM>  {

  private:

  public:
    typedef CellPosition<Field,DIM,Algorithms>  CellPositionType;
    typedef CellPositionType                    ThisType;
    typedef            SeitzPosition<Field,DIM> BaseType;
    typedef Field      ArrayType[DIM];

    /**
     * Default Constructor
     *
     */
    CellPosition() : BaseType () {}

    /**
     * Construct a cell position from a dataset
     *
     * \note When new CellPositions are created they are snapped
     *       back into the lattice's cell if the given values represent positions outside of the Cell.  
     */
    CellPosition(const ArrayType& values): BaseType(values)
    {
      //snap(); This is done externally now 
    }
    
    /**
     * Construct a cell position from a SeitzVector
     *
     * \note When new CellPositions are created they are snapped
     *       back into the lattice's cell if the given values represent positions outside of the Cell.  
     */
    CellPosition(const SeitzVector<Field,DIM,1>& svec): BaseType(svec)
    {
      //snap(); This is done externally now 
    }
    
    ThisType& operator = (const ThisType& v)
    { 
      const BaseType& bv = v;
      BaseType& bt = (*this);
      bt = bv;
      return *this; 
    }

    //============================================================

    /** 
     * Review and if neccesary change all of the positions so that
     * their coefficients are in the range [0,1).
     *
     **/
    ThisType& normalize()  {       // &*&*&*&* this is used a lot probably should do a foreach2 here
      for (size_t i=0; i< DIM; i++) 
	(*this)[i] = Algorithms::modulus((*this)[i],Field(1));
      return (*this);
    }
    
    /** 
     * Determine if a given cell position is close to this position.
     **/
    bool closeTo (const CellPositionType& other) const {
      return CLOSE<BaseType,BaseType,
	typename BaseType::Traits,
	typename BaseType::Traits,
	Algorithms>::EXEC((*this), other);
    }

    template<typename IN_TYPE>
    static
    ThisType make(IN_TYPE x, IN_TYPE y, IN_TYPE z) {
      assert(DIM == 3);
      ThisType result;
      result[0] =  convert<Field>(x);
      result[1] =  convert<Field>(y);
      result[2] =  convert<Field>(z);
      return result;
    }
    
    template<typename IN_TYPE>
    static
    ThisType make(IN_TYPE x, IN_TYPE y) {
      assert(DIM == 2);
      ThisType result;
      result[0] = convert<Field>(x);
      result[1] = convert<Field>(y);
      return result;
    }

  };

  //====================================================================== 


  template<typename Field, typename Algorithms>
  static
  CellPosition<Field,3,Algorithms> cellPosition(Field x, Field y, Field z) {
    CellPosition<Field,3,Algorithms> result;
    result[0] = x;
    result[1] = y;
    result[2] = z;
    return result;
  }

  template<typename Field, typename Algorithms>
  CellPosition<Field,2,Algorithms> cellPosition(Field x, Field y) {
    CellPosition<Field,2,Algorithms> result;
    result[0] = x;
    result[1] = y;
    return result;
  }


  template<typename Field, typename IN_TYPE, typename Algorithms> 
  CellPosition<Field,2,Algorithms> cellPosition(IN_TYPE t0, IN_TYPE t1) {
    CellPosition<Field,2,Algorithms> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    return result;
  }

  template<typename Field, typename IN_TYPE, typename Algorithms> 
  CellPosition<Field,3,Algorithms> cellPosition(IN_TYPE t0, IN_TYPE t1, IN_TYPE t2) {
    CellPosition<Field,3,Algorithms> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    result[2] = convert<Field>(t2);
    return result;
  }
  
  //====================================================================== 

  /**
   * Return the CellTranslation from the difference of two positions.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  CellTranslation<Field,DIM> operator-(const CellPosition<Field,DIM,Algorithms>& cellPosition1, 
				       const CellPosition<Field,DIM,Algorithms>& cellPosition2) {
    const SeitzVector<Field,DIM,1>& v1 = cellPosition1;
    const SeitzVector<Field,DIM,1>& v2 = cellPosition2;

    CellTranslation<Field,DIM> result(v1-v2);
    return result;
  }

  /**
   * Return the CellPosition from a position and a translation.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  CellPosition<Field,DIM,Algorithms> operator+ (const CellPosition<Field,DIM,Algorithms>& cellPosition, 
					        const CellTranslation<Field,DIM>& cellTranslation) {
    const SeitzVector<Field,DIM,1>& p = cellPosition;
    const SeitzVector<Field,DIM,0>& t = cellTranslation;
    return CellPosition<Field,DIM,Algorithms>(t + p);
  }
  /**
   * The less operator.
   *
   *      This less operator takes into account a tolerance
   *      for slightly different positions.
   *
   *  Don't think this is used like this anymore
   *
   *\note This operator is used by map/set to:
   *      - determine if two keys are the same, and
   *      - to determine the order of elements in the map/set. 
   */
  template<typename Field, size_t DIM, typename Algorithms>
  bool operator< (const CellPosition<Field,DIM, Algorithms>& lhs, 
		  const CellPosition<Field,DIM, Algorithms>& rhs) {
      
    if (lhs.closeTo(rhs)) return false;

    for (size_t i = 0; i<DIM; i++) {
      Field d1 = lhs[i] - rhs[i];
      if (d1 == 0.0) continue;
      return d1 < 0.0;
    }

    return false;
  }

  /** \ingroup ostream
   *
   * SeitzVector output stream operator 
   **/
  template<typename Field,size_t DIM, int IND, typename Algorithms>
  std::ostream& operator << (std::ostream& os, const CellPosition<Field,DIM, Algorithms>& p) {
    
    typedef CellPosition<Field,DIM, Algorithms> CellPositionType;

    MAT_PRINT<CellPositionType,ColMajorTraits<Field,DIM,1> >::EXEC(p);
    os << "CellPosition";
    return os;
  }

  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<typename Field, size_t DIM, typename Algorithms> 
  Tag toXML(const CellPosition<Field,DIM,Algorithms> pos, 
	    std::string name="CellPosition") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << pos[i] ;
    return result;
  }

  //   /** \ingroup ostream
  //    *
  //    * Output stream operator 
  //    */
  //   template<typename Field,size_t DIM>
  //   std::ostream& operator << (std::ostream& os, const CellPositionType& t) {
    
  //     os << "[ ";
  //     // should pass off << to Field overload
  //     // could print SymOp symbol here
  //     os.setf(std::ios_base::fixed, std::ios_base::floatfield);
  //     os.precision(6);
  //     for(size_t i=0; i< DIM; ++i) 
  //       if (i<DIM-1)
  // 	os << t[i] << ", ";
  //       else
  // 	os << t[i] << " ";

  //     os << "]";
    
  //     return os;
  //   };
  
  //   /**
  //    * Set the value of the default tolerance.
  //    */
  //   template<typename Field,size_t DIM>
  //   const Field CellPosition<Field, DIM>::defaultTolerance = Field(0.0000001);
 
  //   /**
  //    * Set the initial value of the tolerance.
  //    */
  //   template<typename Field,size_t DIM>
  //   Field CellPosition<Field, DIM>::tolerance = Field(0.0000001);
 
} /* namespace psimag */

#endif // PSIMAG_CellPosition_H

/*@}*/

