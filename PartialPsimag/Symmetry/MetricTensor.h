//-*-C++-*-

#ifndef PSIMAG_MetricTensor_H
#define PSIMAG_MetricTensor_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file MetricTensor.h Contains a class definition for MetricTensor objects. */

#include <cstddef>
#include <limits>
#include <list>
#include <vector>
#include <algorithm>

#include "Vec.h"
#include "Mat.h"
#include "Real.h"
#include "Lattice.h"
#include "SeitzVectors.h"

#include "PSIMAGAssert.h"

namespace psimag {

  //======================================================================
  
  /** Forward Declaration of class Lattice. **/
  template <typename Field, size_t DIM, typename Algorithms> class Lattice;

  /** Forward Declaration of class MetricTensor. **/
  template <typename Field, size_t DIM> class MetricTensor;

  //======================================================================

  /** \ingroup latticeConcepts
   * \brief Used to sort basis vectors according to a canonical order  (see 2002 ITC, pg 750).
   */
  template <typename Field, size_t DIM> 
  class MetricTensorSortHelper {
  private:
    const MetricTensor<Field,DIM>& metricTensor;
  public:
    MetricTensorSortHelper(const MetricTensor<Field,DIM>& mt): metricTensor(mt) {}

    bool operator() ( const size_t& i, const size_t& j) {
      return metricTensor.compareAxisPositions(i,j);
    }
  };

  //====================================================================== General

  /** \ingroup latticeConcepts
   *
   * \brief A symmetric matrix composed of the inner products of the
   *        CartesianTranslation vectors of a given cell.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM: Dimensionality of the lattice being
   *             represented. Legal values are 1,2,3.
   */
  template <typename Field, size_t DIM>
  class MetricTensor: public Mat<Field,DIM,DIM>  {
    MetricTensor(): Mat<Field,DIM,DIM>() {}
  };
    
  //====================================================================== 3D

  /** \ingroup latticeConcepts
   *
   * \brief A 3D specialization of symmetric matrix composed of the
   *        inner products of the CartesianTranslation vectors of a
   *        given cell.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM: Dimensionality of the lattice being
   *             represented. Legal values are 1,2,3.
   */
  template <typename Field>
  class MetricTensor<Field,3>: public Mat<Field,3,3>  {
    
  public:
    
    /* Don't call the default constructor, use the metricTensor generic function. */
    MetricTensor(): Mat<Field,3,3>() {}

    /** Return the A component of this MetricTensor. (see 2002 ITC, pg 753) */
    Field  A() const { return (*this)(0,0);  }

    /** Return the B component of this MetricTensor.  */
    Field  B() const { return (*this)(1,1); }

    /** Return the C component of this MetricTensor. */
    Field  C() const { return (*this)(2,2); }

    /** Return the D component of this MetricTensor. */
    Field  D() const { return (*this)(1,2); }

    /** Return the E component of this MetricTensor. */
    Field  E() const { return (*this)(0,2); }

    /** Return the F component of this MetricTensor. */
    Field  F() const { return (*this)(0,1); }

    /** Return the cosine of the angle between the b and c basis vectors. */
    Field  cos_alpha() const { return D() / (sqrt(B())*sqrt(C()));}

    /** Return the the angle between the b and c basis vectors. */
    Field  alpha()     const { return acos(cos_alpha())*180.0/PI; }

    /** Return the cosine of the angle between the a and c basis vectors. */
    Field  cos_beta () const { return E() / (sqrt(C())*sqrt(A()));}
    
    /** Return the the angle between the a and c basis vectors. */
    Field  beta()      const { return acos(cos_beta())*180.0/PI; }
       
    /** Return the cosine of the angle between the a and b basis vectors. */
    Field  cos_gamma () const { return F() / (sqrt(A())*sqrt(B()));}

    /** Return the the angle between the a and b basis vectors. */
    Field  gamma()      const { return acos(cos_gamma())*180.0/PI; }

    /**
     * Return true if the two axes given by the two input indexes are in canonical order.
     */
    bool compareAxisPositions ( const size_t& i, const size_t& j) const {
      
      if ( i == j ) return false;
      
      if ( (*this)(i,i) < (*this)(j,j) )
	return true;
      
      if ( (*this)(i,i) > (*this)(j,j) )
	return false;
      
      // *this(i,i) == *this(j,j)
      size_t k = otherIndex(i, j );
      return (*this)(j,k) < (*this)(i,k);

      return false;
    }

    /**
     * Given two indicies return the other possible index.
     * (e.g. given 0 and 2 return 1)
     */
    size_t otherIndex(size_t i, size_t j) {
      
      assert( i != j);
      
      if ( i == 0 ) {
	if ( j == 1) return 2;
	if ( j == 2) return 1;
      }
      if ( i == 1 ) {
	if( j == 0) return 2;
	if( j == 2) return 1;
      }
      if ( i == 2 ) {
	if( j == 0) return 1;
	if( j == 1) return 0;
      }
      return 0; // Should never get here because of asserts
    }

    /** \ingroup xml
     *
     * XML Output function for 3D CellParameters.
     */
    Tag toXML() {
      
      Tag result("MetricTensor");
      result["DIM"]    = 3;
      result.content << ( Mat< Field, 3, 3  > ) (*this);
      return result;
    }

  };
  //======================================================================

  /** \ingroup latticeConcepts
   *
   * \brief A 2D specialization of symmetric matrix composed of the
   *        inner products of the CartesianTranslation vectors of a
   *        given cell.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM: Dimensionality of the lattice being
   *             represented. Legal values are 1,2,3.
   */
  template <typename Field>
  class MetricTensor<Field,2>: public Mat<Field,2,2>  {
    
  public:

    typedef  Mat<Field,2,2> MatType;
    
    /* Don't call the default constructor, use the metricTensor generic function. */
    MetricTensor(): Mat<Field,2,2>() {}

    /** Return the A component of this MetricTensor. (see 2002 ITC, pg 753) */
    Field  A() const { return (*this)(0,0);  }

    /** Return the B component of this MetricTensor.  */
    Field  B() const { return (*this)(1,1); }

    /** Return the C component of this MetricTensor. */
    Field  C() const { return (*this)(0,1); }

    /** Return the cosine of the angle between the a and b basis vectors.  */
    Field  cos_alpha () const { return C() / (sqrt(A())*sqrt(B()));}

    /** Return the the angle between the a and b basis vectors.   */
    Field  alpha()      const { return acos(cos_alpha())*180.0/PI; }

    /** Return true if the two axes given by the two input indexes are in canonical order. */
    bool compareAxisPositions ( const size_t& i, const size_t& j) const {
      
      if ( i == j ) return false;
      
      if ( (*this)(i,i) < (*this)(j,j) )
	return true;
      
      if ( (*this)(i,i) > (*this)(j,j) )
	return false;
      
      return false;
    }

  };

  //======================================================================
    
  /** 
   *  \brief Return a vector of ints indicating the canonical order
   *         of the translations of a Lattice corresponding to this
   *         MetricTensor.
   *   
   *   This is generally the order of the length of the translations
   *   with special case rules to be used when the vectors are of
   *   equal length.
   *
   *   Care should be taken when using this order since the ordering produced
   *   may be numerically sensitive to the Lattice's basis vectors, parameters, etc.
   *
   */
  template<typename Field, size_t DIM>
  Vec<size_t, DIM > axisOrder(MetricTensor<Field,DIM> metric) {
    
    Vec<size_t, DIM> order;
    for (size_t i = 0; i < DIM; i++)
      order[i] = i;

    MetricTensorSortHelper<Field,DIM> cmp(metric);

    std::sort(order.begin(), order.end(), cmp); 
      
    return order;
  }

  //======================================================================
  
  /**
   * Construct a MetricTensor that corresponds to a given Lattice.
   */
  template<typename Field,size_t DIM,typename Algorithms>
  MetricTensor<Field,DIM> metricTensor(Lattice<Field,DIM,Algorithms>& cell) {
    
    typedef Lattice<Field,DIM,Algorithms>     LatticeType;
    typedef MetricTensor<Field,DIM>           MetricTensorType;
    typedef typename LatticeType::MatType     LatticeMatType;
    typedef typename LatticeMatType::Traits   LatticeMatTraits;
    typedef typename LatticeType::BasisVectors BasisVectors;

    MetricTensorType result;

    BasisVectors translations(DIM);
    
    // Unpack the cell (as a matrix) into the basis vectors 
    COPY<BasisVectors, LatticeMatType,ReverseDoubleIndexTraits<Field,DIM,DIM>,LatticeMatTraits>::EXEC(translations, cell);

    // Make a template for this operation?
    for(size_t i=0; i < DIM; i++) {
      for(size_t j=0; j <= i; j++) {
	result(i,j) = translations[i] * translations[j];
	if (i == j) continue;
	result(j,i) = result(i,j);
      }
    }
    return result;
  }

  //======================================================================
  
  /** \ingroup xml
   *
   * XML Output function for 2D CellParameters.
   */
  template<typename Field, size_t DIM> 
  Tag toXML(MetricTensor<Field,DIM> metric,
	    std::string name = "MetricTensor") {
    
    Tag result(name);
    result["DIM"]    = 2;
     
    typedef  Mat<Field,2,2> MatType;

    std::ostringstream buff;
    MAT_PRINT<MatType,typename MatType::Traits>::JUST_NUMBERS( metric, buff);
    result.content << buff.str();
    
    return result;
  }

  //======================================================================
  
  template<typename Field, size_t DIM> 
  std::ostream& operator << (std::ostream& os, 
			     MetricTensor<Field,DIM> metric) {
    os << "------------------------------------ MetricTensor<" << DIM << ">:" << std::endl;
    os <<  ( Mat< Field, DIM, DIM  > ) metric << std::endl;
//     os << "------------------------------------     MetricTensor<" << DIM << ">";
//     os << " A()         = " << metric.A() << std::endl;
//     os << " B()         = " << metric.B() << std::endl;
//     os << " C()         = " << metric.C() << std::endl;
//     os << " cos_alpha() = " << metric.cos_alpha() << std::endl;
//     os << " cos_alpha() = " << metric.cos_alpha() << std::endl;
    os << "------------------------------------ End MetricTensor<" << DIM << ">";
    return os;
  }

} /* namespace psimag */

#endif
/*@}*/

