//-*-C++-*-
#ifndef  Psimag_CellTranslation
#define  Psimag_CellTranslation

/** \ingroup extendedMatrices */
/*@{*/

/** \file CellTranslation.h
 *
 * Contains a class for implementing CellTranslation (SeitzTranslation subclass) objects.
 */

#include "SeitzVector.h"
#include "LatticeCoordinates.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Algorithms> class CellPosition;
  
  /** \ingroup extendedMatrices
   *
   * \brief A marker class indicating that the Translation indicates a direction from within a Cell.
   *
   * These can be converted to  and from MillerIndices for a given Cell.
   */
  template<typename Field, size_t DIM> 
  class CellTranslation: public SeitzVector< Field, DIM, 0 > {
    
  public:

    typedef          CellTranslation<Field, DIM> ThisType;
    typedef          SeitzVector<Field,DIM,0>    BaseType;
    typedef typename BaseType::ArrayType         ArrayType;

    // Constructors are never inherited :-(

    /** The Default Constructor produces the zero translation. */
    CellTranslation(): BaseType() {}

    /** Construct a translation whose components are set to the given value. */
    //template<typename IN_TYPE> CellTranslation(const IN_TYPE& val): SeitzVector<Field, DIM, 0>(val) {}

     /** Construct a translation whose components are set to the given value. */
    CellTranslation(const Vec<Field, DIM>& v):  BaseType(v) {}

     /** Construct a translation whose components are set to the given value. */
    CellTranslation(const BaseType& v):  BaseType(v) {}

    /** Construct a translation whose components are set to the given value. */
    CellTranslation(const ThisType& v): BaseType(v) {}

    /** Construct a translation whose components are set to the given value array. */
    CellTranslation(const ArrayType& vals):  BaseType(vals) {}

    /** Construct a translation whose components are set to the given value array. */
    template<typename Algorithms>
    explicit CellTranslation(const CellPosition<Field,DIM,Algorithms> p):  BaseType(p) {}

    ThisType& operator = (const ThisType& v)
    { 
      const BaseType& bv = v;
      BaseType& bt = (*this);
      bt = bv;
      return *this; 
    }
  };

  /*
   * \brief If the translation is on the negative side of the y-z plane return true.
   *
   */
  template<typename Field, size_t DIM, int IND, typename Algorithms>
  bool isNegative(const SeitzVector<Field,DIM,IND>& v) {
    static Field                      zero(0);
    for(size_t i = 0; i<DIM; i++) 
      if(!Algorithms::close(v[i],zero))
	return v[i] < zero;
    return false; // all zeros
  }

  template<typename Field, size_t DIM, int IND, typename Algorithms>
  bool closeTo(const SeitzVector<Field,DIM,IND>& v1,
	       const SeitzVector<Field,DIM,IND>& v2) {
    return CLOSE<SeitzVector<Field,DIM,IND>, SeitzVector<Field,DIM,IND>,
      typename SeitzVector<Field,DIM,IND>::Traits,
      typename SeitzVector<Field,DIM,IND>::Traits,
      Algorithms>::EXEC(v1, v2);
  }

  /* \brief Subtract the largest net translation that the given cell
   *        translation extends.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  CellTranslation<Field,DIM> normalizeTranslation(CellTranslation<Field,DIM>& trans) {
    static Field one(1);
    for (size_t i=0; i< DIM; i++) 
      trans[i] = Algorithms::modulus(trans[i],one);
    return trans;
  }

  /* \brief Subtract the largest net translation that the given cell
   *        translation extends.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  void netReduce(CellTranslation<Field,DIM>& trans) {
    

    try {
      LatticeCoordinates<DIM> latticeCoordinate = maxContainedLatticeCoordinate<Field,Algorithms>(trans);
      COPY<CellTranslation<Field,DIM>,
	LatticeCoordinates<DIM>,
	typename CellTranslation<Field,DIM>::Traits,
	typename LatticeCoordinates<DIM>::Traits>::EXEC_MINUS(trans,latticeCoordinate);

    }
    catch (std::range_error& e) {
      // If there is a range_error then the trans is already normalized!
    }
    
  }

  /* \brief Subtract the largest net translation that the given cell
   *        translation extends.
   */
  template<typename Field, size_t DIM, typename Algorithms>
  void normalizeGlide(CellTranslation<Field,DIM>& trans) {

    netReduce<Field,DIM,Algorithms>(trans);

    if(isNegative<Field,DIM,0,Algorithms>(trans)) {

      LatticeCoordinates<DIM> netDirection = latticeCoordinate<Field,Algorithms>(trans);
      COPY<CellTranslation<Field,DIM>,
	LatticeCoordinates<DIM>,
	typename CellTranslation<Field,DIM>::Traits,
	typename LatticeCoordinates<DIM>::Traits>::EXEC_MINUS(trans,netDirection);

    }
  }

  //======================================================================
  /**
   * \brief The vector times a scalor.
   */
  template<typename Field,size_t DIM>
  CellTranslation<Field,DIM>  operator * (const CellTranslation<Field,DIM>& lhs, Field scalar){
    CellTranslation<Field,DIM> result;
    for (size_t i=0; i< DIM;i++)
      result[i] = lhs[i] * scalar;
    return result;
  }

  template<typename Field, typename IN_TYPE> 
  CellTranslation<Field,2> cellTranslation(IN_TYPE t0, IN_TYPE t1) {
    CellTranslation<Field,2> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    return result;
  }

  template<typename Field, typename IN_TYPE> 
  CellTranslation<Field,3> cellTranslation(IN_TYPE t0, IN_TYPE t1, IN_TYPE t2) {
    CellTranslation<Field,3> result;
    result[0] = convert<Field>(t0);
    result[1] = convert<Field>(t1);
    result[2] = convert<Field>(t2);
    return result;
  }
  
  //====================================================================== 

  /** \ingroup xml
   *
   * XML Output function for CartesianPositions.
   */
  template<typename Field, size_t DIM> 
  Tag toXML(const CellTranslation<Field,DIM> t, std::string name="CellTranslation") {
      
    Tag result(name);
    for (size_t i=0; i<DIM; i++)
      result.content <<  " " << t[i] ;
    return result;
  }


}  /* namspace psimag */

#endif // Psimag_Miller_Direction
/*@}*/
