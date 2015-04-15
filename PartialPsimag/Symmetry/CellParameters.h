//-*-C++-*-

#ifndef PSIMAG_CellParameters_H
#define PSIMAG_CellParameters_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file CellParameters.h Class definitions for CellParameter objects. */ 

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Mat.h"
#include "Real.h"
#include "Tag.h"

#include "SeitzVectors.h"
#include "MetricTensor.h"
#include "Reciprocal.h"

#include "PSIMAGAssert.h"

namespace psimag {

  template<typename Field>
  static Field radians(Field degrees) {
    return degrees/Field(180)* PI;
  }

  template<typename Field>
  static Field degrees(Field radians) {
    return (radians/PI)*Field(180);
  }

  /** \ingroup latticeConcepts
   *
   *  \brief Primary Template for CellParameter classes. Don't use
   *         this template, it is meant to be specialized.
   * 
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   *  \param DIM: Dimensionality of the lattice being
   *              represented. Legal values are 1,2,3.
   */
  template <typename Field, size_t DIM, typename Algorithms> class CellParameters {};
    
  //======================================================================
  // 1D case
  //
  /** \ingroup latticeConcepts
   *
   *  \brief A parameter that can be used to define a 1 dimensional
   *         Cell.
   * 
   *  The parameter is the lengths of the cell's translation vector.
   *
   *  \note Given a CellParameters object a canonical set of cell translation vectors
   *        can be constructed. 
   *
   *        -# The a vector is put on the x axis.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   */
  template <typename Field, typename Algorithms>
  class CellParameters<Field,1,Algorithms> {
    
  public:


    Field a;          /**< The length of the shortest translation vector. **/

    /**
     * Default Constructor, all parameters set to zero.
     * Most of the time this is not what you want to use.
     */
    CellParameters():
      a    (Field(1))
    {}

    /**
     * Construct CellParameters from a given set of parameters.
     * Perform type translations as necessary.
     */
    template<typename In_Type>
      CellParameters(In_Type _a) :
	a    (convert<Field>(_a))
    {}

    /**
     * Construct CellParameters from a given MetricTensor.
     */
    CellParameters(MetricTensor< Field, 1>& metric) {
      a = sqrt(metric.A());
    }

    std::string typeString() const {
      std::string result = "{ TDB }";
      result += "}";
      return result;
    }

    /**
     * \brief Setup the basis vectors of the given cell to reflect the
     *        parameters of this cell.
     */
    template<typename LAlgorithms>
    void makeBasisFor(Lattice<Field,1,LAlgorithms>& cell) {
      
      cell(0,0) = a;
    }

   };

  //====================================================================== 
  
  /** \ingroup xml
   *
   * XML Output function for 3D CellParameters.
   */
  template<typename Field, typename Algorithms >
  Tag toXML(const CellParameters< Field, 1, Algorithms >& p,
	    std::string name="CellParameters") {
      
      Tag result(name);
      result["a"]     = p.a;

      return result;
    }

  /** \ingroup ostream
   *
   * Output operator for CellParameters.
   */
  template<typename Field,typename Algorithms>
  std::ostream& operator << (std::ostream& os, const CellParameters< Field, 1, Algorithms >& p) {
    
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    
    os << "[ " << p.a << " ]" << std::endl;
    
    return os;
  }

  //======================================================================
  // 2D case
  //
  /** \ingroup latticeConcepts
   *
   *  \brief Three parameters that can be used to define a 2 dimensional
   *         Cell up to a change of orientation.
   * 
   *  The parameters are the lengths of the cell's translation vectors,
   *  a and b and the angle between them, alpha.
   *
   *  \note Given a CellParameters object a canonical set of cell translation vectors
   *        can be constructed. 
   *
   *        -# The a vector is put on the x axis.
   *        -# The b vector is put in the x-y plane at an angle of alpha to the a vector.
   *
   *  \note Small changes in either the CellParameters or in the basis
   *        vectors can result in large changes in reduced cell and the
   *        BravaisType of the lattice.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   */
  template <typename Field, typename Algorithms>
  class CellParameters<Field,2, Algorithms> {
    
  public:

    static bool close(Field v1, Field v2) {
      return Algorithms::close(v1,v2);
    }

    Field a;          /**< The length of the      shortest translation vector. **/
    Field b;          /**< The length of the next shortest translation vector. **/
    Field alpha;      /**< The angle between b and c. **/
    Field sin_alpha;  /**< sin(alpha) **/
    Field cos_alpha;  /**< cos(alpha) **/


    /**
     * Default Constructor, all parameters set to zero.
     * Most of the time this is not what you want to use.
     */
    CellParameters():
	a    (Field(1)),
	b    (Field(1)),
	alpha(Field(0)),
	sin_alpha(Field(0)),
	cos_alpha(Field(1))
    {}

//     /**
//      * Default Constructor, all parameters set to zero.
//      * Most of the time this is not what you want to use.
//      */
//     CellParameters(const CellParameters<Field,2>& p):
//       a(p.a),
//       b(p.b),
//       alpha(p.alpha)
//     {}

    /**
     * Construct CellParameters from a given set of parameters.
     * Perform type translations as necessary.
     */
    template<typename In_Type>
      CellParameters(In_Type _a, In_Type _b, In_Type _alpha) :
	a        (convert<Field>(_a)),
	b        (convert<Field>(_b)),
	alpha    (convert<Field>(_alpha)),
	sin_alpha(sin(radians(alpha))),
	cos_alpha(cos(radians(alpha)))
    {
      assert(alpha < convert<Field>(180));
    }

    /**
     * Construct CellParameters from a given MetricTensor.
     */
    CellParameters(MetricTensor< Field, 2 >& metric):
      a        (sqrt(metric.A())),
      b        (sqrt(metric.B())),
      alpha    (metric.alpha()),
      sin_alpha(sin(radians(alpha))),
      cos_alpha(cos(radians(alpha)))
    {}

    bool areSquare()        const { return  close(a,b) &&  close(alpha,90.0); }
    bool areRectangular()   const { return !close(a,b) &&  close(alpha,90.0); }
    bool areOblique()       const {
      if (close(a,b))
	return (!close(alpha,90.0) && !close(alpha,60.0) && !close(alpha,120.0));
      else
	return !close(alpha,90.0); 
    }
    bool areHexagonal()     const { return  close(a,b) && (close(alpha,60.0) || close(alpha,120.0)); }

    bool areRhomboherdal()  const { return bSidedCenteredRectangle() || aSidedCenteredRectangle(); }

    bool aSidedCenteredRectangle() const {
      return ( close(a,cos_alpha*2.0*b)  );
    }

    bool bSidedCenteredRectangle() const {
      return ( close(b,cos_alpha*2.0*a)  );
    }

    std::string typeString() const {
      std::string result = "{";
      if (areSquare())       result += "Square,";
      if (areRectangular())  result += "Rectangular,";
      if (areOblique())      result += "Oblique,";
      if (areHexagonal())    result += "Hexagonal,";
      if (areRhomboherdal()) result += "Rhomboherdal,";
      result += "}";
      return result;
    }

    /**
     * \brief Setup the basis vectors of the given cell to reflect the
     *        parameters of this cell.
     */
    template<typename LAlgorithms>
    void makeBasisFor(Lattice<Field,2,LAlgorithms>& cell) const {
      
      // column 1 is basis vector 1
      cell(0,0) = a;
      cell(1,0) = Field(0);

      // column 2 is basis vector 1
      cell(0,1) = cos_alpha * b;
      cell(1,1) = sin_alpha * b;

    }

  };

  //====================================================================== 
  
  /** \ingroup xml
   *
   * XML Output function for 3D CellParameters.
   */
  template<typename Field, typename Algorithms >
  Tag toXML(const CellParameters< Field, 2, Algorithms >& p,
	    std::string name="CellParameters") {
      
      Tag result(name);
      result["a"]     = p.a;
      result["b"]     = p.b;
      result["alpha"] = p.alpha;

      return result;
    }

  /** \ingroup ostream
   *
   * Output operator for CellParameters.
   */
  template<typename Field, typename Algorithms >
  std::ostream& operator << (std::ostream& os, const CellParameters< Field, 2 , Algorithms >& p) {
    
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    
    os << "[ " 
       << p.a << ","
       << p.b << ","
       << p.alpha 
       << " ]" ;
    
    return os;
  }

  //======================================================================
  // 3D case
  //
  /** \ingroup latticeConcepts
   *
   *  \brief Six parameters that can be used to define a 3 dimensional
   *         Cell up to a change of orientation.
   * 
   *  The parameters are the lengths of the cell's translation vectors,
   *  a, b and c and the angles between them alpha, beta, gamma.
   *
   *  \note Given a CellParameters object a canonical set of cell translation vectors
   *        can be constructed. 
   *
   *        -# The a vector is put on the x axis.
   *        -# The b vector is put in the x-y plane.
   *        -# The c vector's cartesian coordinate are determined by the
   *           placement of the a and b vectors.
   *
   *  \note Small changes in either the CellParameters or in the basis
   *        vectors can result in large changes in reduced cell and the
   *        BravaisType of the lattice.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   */
  template <typename Field, typename Algorithms>
  class CellParameters<Field,3, Algorithms> {
    
  public:

    Field a;          /**< The length of the      shortest translation vector. **/
    Field b;          /**< The length of the next shortest translation vector. **/
    Field c;          /**< The length of the       longest translation vector. **/
    Field alpha;      /**< The angle between b and c. **/
    Field beta;       /**< The angle between a and c. **/
    Field gamma;      /**< The angle between a and b. **/

    Field sin_alpha;  /**< sin(alpha) **/
    Field sin_beta;   /**< sin(beta) **/
    Field sin_gamma;  /**< sin(gamma) **/
    Field cos_alpha;  /**< cos(alpha) **/
    Field cos_beta;   /**< cos(beta) **/
    Field cos_gamma;  /**< cos(gamma) **/

    /**
     * Default Constructor, 
     * Most of the time this is not what you want to use.
     */
    CellParameters():
	a        (Field(1)),
	b        (Field(1)),
	c        (Field(1)),
	alpha    (Field(90)),
	beta     (Field(90)),
	gamma    (Field(90)),
	sin_alpha(sin(radians(alpha))),
	cos_alpha(cos(radians(alpha))),
	sin_beta (sin(radians(beta))),
	cos_beta (cos(radians(beta))),
	sin_gamma(sin(radians(gamma))),
	cos_gamma(cos(radians(gamma)))
    {}

    /**
     * Construct CellParameters from a given set of parameters.
     * Perform type translations as necessary.
     */
    template<typename In_Type>
      CellParameters(In_Type _a,    In_Type _b,    In_Type _c, 
		    In_Type _alpha, In_Type _beta, In_Type _gamma) :
	a    (convert<Field>(_a)), 
	b    (convert<Field>(_b)), 
	c    (convert<Field>(_c)),
	alpha(convert<Field>(_alpha)),
	beta (convert<Field>(_beta)),
	gamma(convert<Field>(_gamma)),
	sin_alpha(sin(radians(alpha))),
	cos_alpha(cos(radians(alpha))),
	sin_beta (sin(radians(beta))),
	cos_beta (cos(radians(beta))),
	sin_gamma(sin(radians(gamma))),
	cos_gamma(cos(radians(gamma)))
    {
      assert(alpha < convert<Field>(180));
      assert(beta  < convert<Field>(180));
      assert(gamma < convert<Field>(180));

    }

    /**
     * Construct CellParameters from a given MetricTensor.
     */
    CellParameters(MetricTensor< Field, 3 >& metric):
      a    (sqrt(metric.A())), 
      b    (sqrt(metric.B())), 
      c    (sqrt(metric.C())),
      alpha(metric.alpha()),
      beta (metric.beta()),
      gamma(metric.gamma()),
      sin_alpha(sin(radians(alpha))),
      cos_alpha(cos(radians(alpha))),
      sin_beta (sin(radians(beta))),
      cos_beta (cos(radians(beta))),
      sin_gamma(sin(radians(gamma))),
      cos_gamma(cos(radians(gamma)))
    {}

    std::string typeString() const {
      std::string result = "{ TDB }";
      result += "}";
      return result;
    }

    /**
     * \brief Setup the basis vectors of the given cell to reflect the
     *        parameters of this cell.
     */
    template<typename LAlgorithms>
    void makeBasisFor(Lattice<Field,3,LAlgorithms>& cell) {
      
      cell(0,0) = a;
      cell(0,1) = Field(0);
      cell(0,2) = Field(0);
      cell(1,0) = cos_gamma * b;
      cell(1,1) = sin_gamma * b;
      cell(1,2) = Field(0);
      cell(2,0) = c * cos_beta;
      cell(2,1) = - c * sin_beta * cos_alpha;
      cell(2,2) = 1 / c;
    }
 
   };

  //====================================================================== 
  
  /** \ingroup xml
   *
   * XML Output function for 3D CellParameters.
   */
  template<typename Field, typename Algorithms >
  Tag toXML(const CellParameters< Field, 3 , Algorithms >& p,
	    std::string name="CellParameters") {
      
      Tag result(name);
      result["a"]     = p.a;
      result["b"]     = p.b;
      result["c"]     = p.c;
      result["alpha"] = p.alpha;
      result["beta"]  = p.beta;
      result["gamma"] = p.gamma;

      return result;
    }

  //====================================================================== 

  /** \ingroup ostream
   *
   * Output operator for CellParameters.
   */
  template<typename Field, typename Algorithms >
  std::ostream& operator << (std::ostream& os, const CellParameters< Field, 3 , Algorithms >& p) {
    
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    
    os << "[ " 
       << p.a << ","
       << p.b << ","
       << p.c << ","
       << p.alpha << ","
       << p.beta << ","
       << p.gamma << ","
       << " ]" ;
    
    return os;
  }

} /* namespace psimag */

#endif

/*@}*/

