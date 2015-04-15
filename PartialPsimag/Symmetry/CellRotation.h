//-*-C++-*-

#ifndef  Psimag_CellRotation
#define  Psimag_CellRotation

/** \ingroup extendedMatrices */
/*@{*/

/** \file CellRotation.h
 *
 *  Contains a class for Matricies that represent the rotation part of
 *  Symmetry Operations, et. al.
 */
 
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

#include "SeitzVectors.h"
#include "Lattice.h"
#include "CellDirection.h"
#include "CellPosition.h"
#include "CartesianRotation.h"
#include "Mat.h"
#include "PSIMAGAssert.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Algorithms> class Lattice;
  
  /** \ingroup extendedMatrices 
   *
   *\brief A subclass of the integer matrix, Mat<Field, DIM, DIM>.
   *       which is used to represent the rotation part of symmetry
   *       operations.  (See 2002 ICT Section 5.2)
   *
   * These rotations are always expressed relative to the coordinate
   * system of some cell. They sometimes also are limited to the representation
   * of rotations about axes which are CellDirections with integer
   * cooeficients. This obviously includes the CellDirections
   * corresponding to the cells basis vectors.
   *
   * For the corresponding CartesianRotation object call
   * cartesianRotation(cell).
   *
   * CellRotation objects keep track of their inverses, produced by lazy
   * evaluation, so they can be reused if necessary and these inverses
   * destroyed with the rotation.
   *
   *  see r3_rotation.h in cctbx
   *
   */
  template<typename Field, size_t DIM>
  class CellRotation: public Mat< Field, DIM, DIM > {
    
  public:

    typedef Mat< Field, DIM, DIM > RotationClass;

    typedef enum { sixfoldInv=-6, fourfoldInv=-4, threefoldInv=-3, twofoldInv=-2, inversion=-1, identity=1, 
		   twofold=2, threefold=3,  fourfold=4, sixfold=6 } RotationType;
  private:
    
    CellRotation<Field, DIM>* inv;
    bool        invSet;

    CellDirection<Field, DIM> axis;

  protected:

    /** 
     * \brief Perform checks to see if this matrix is a valid rotation
     *        expressed in cell coordinates.
     *
     *  - Check to see if the determinant == 1 or -1;
     *  - Check to see if the trace is in the right range.
     *  - Possibly check to see it we are an integer matrix
     *
     * \throw Throw a range_error if the check fails.
     *
     * \note This should only be called by constructors.
     **/
    void check() {
      
      Field det = Det((*this));
      if ( det != 1 and det != -1 ) {
	std::string message("ERROR CellRotation: ");
	message += " Tried to construct a CellRotation with a Det != +/- 1.";
	throw std::range_error(message);
      }
      
      Field tr  = Trace((*this));
      if ( tr < -3 or tr > 3 ) {
	std::string message("ERROR SymmetryOperation: ");
	message += " Tried to construct a CellRotation with a Trace out of range.";
	throw std::range_error(message);
      }
    }

  public:
    
    /** The destructor handles deleting 'embedded' inverses. **/
    ~CellRotation() {
      if (invSet)
	delete inv;
    }

    // ================================================== Constructors
    
    /** 
     * The default constructor results in an identity Matrix. 
     *
     **/
    CellRotation() : 
      invSet(false),               
      RotationClass( 0 ) 
    {
      MakeIdentity((RotationClass)*this);
    }

    /** 
     * Construct a rotationType and an extended matrix from a by DIMxDIM data array.
     *
     * The non-extended part is set to the given data. Data types are
     * converted as necessary.
     *
     **/
    CellRotation(RotationType t, int sense, const CellDirection<Field, DIM>& axisDirection): invSet(false) {

      throw "CellRotation(RotationType t, MillerIndices<DIM>& axisDirection): not coded";
      // no need for check();
    }

    /** 
     * Construct from a DIM by DIM  array.
     *
     **/
    template<typename In_Type>
    CellRotation(const In_Type data [DIM*DIM]) : 
      invSet(false),
      RotationClass(data)
    {
      check();
    }

    /** 
     * Construct from a DIM by DIM  array.
     *
     **/
    template<typename In_Type>
    CellRotation(const In_Type& val) : 
      invSet(false),
      RotationClass(val)
    {
      check();
    }

    //=================================================== Conversions

    /**
    * \brief Given the Lattice in whose coordinate system this rotation
    *        is defined in, Return the corresponding rotation in
    *        cartesian coordinates.
    *
    */
    template<typename In_Field, typename Algorithms>
    CartesianRotation<Field, DIM> cartesianRotation(Lattice<In_Field,DIM,Algorithms> cell) {
      
      CartesianRotation<Field, DIM> result;
      throw "Rotation::cartesianRotation(Lattice lattice) not coded yet.";
      return result;
    }

    // ================================================== Utility functions
    /**
    * \brief Return an integer representing the type of rotation
    *        contained in this SymmetryOperation.
    *
    * \throw If the det and the trace do not together represent a
    *        valid combination throw a range_error.
    */
    RotationType type() {

      static std::string message("ERROR SymmetryOperation:rotationPartType bad Det/Trace combination.");

      Field det = Det((*this));
      Field tr  = Trace((*this));

      if (tr == Field(-3)) {
	if (det == Field(-1))
	  return -1;
	if (det == Field(1))
	  throw std::range_error(message);
      }
      if (tr == Field(-2)) {
	if (det == Field(-1))
	  return -6;
	if (det == Field(1))
	  throw std::range_error(message);
      }
      if (tr == -1) {
	if (det == Field(-1))
	  return Field(-4);
	if (det == Field(1))
	  return Field(2);
      }
      if (tr == 0) {
	if (det == Field(-1))
	  return Field(-3);
	if (det == Field(1))
	  return Field(3);
      }
      if (tr == 1) {
	if (det == Field(-1))
	  return Field(-2);
	if (det == Field(1))
	  return Field(4);
      }
      if (tr == 2) {
	if (det == Field(-1))
	  throw std::range_error(message);
	if (det == Field(1))
	  return Field(6);
      }
      if (tr == 3) {
	if (det == Field(-1))
	  throw std::range_error(message);
	if (det == Field(1))
	  return Field(6);
      }
      throw std::range_error(message);
    }

    // ================================================== get axis

    /** 
     * Return a Vec containing the axis of rotation of this CellRotation.
     */
    CellDirection<Field, DIM> getAxis() {

      CellDirection<Field, DIM> result;
      return result;
    }

    // ================================================== Inverse
    /** 
     * Return the determinan of the rotation part.
     */
    const Field determinant() {
      return Det(*this);
    }

    // ================================================== Inverse
    /** 
     * Return by lazy evaluation a reference to the inverse of this CellRotation.
     */
    const CellRotation<Field, DIM>& inverse() {

      if (invSet) 
	return *inv;

      inv    = new CellRotation<Field, DIM>();
      invSet = true;

      Inverse<Field>(*this, inv);

      return *inv;
    }


//     // ================================================== Operators
//     /**
//      * \brief  The operator for multiplying two CellRotations
//      */
//     CellRotation<Field, DIM> operator*(const CellRotation<Field, DIM>& S) {
//       CellRotation<Field, DIM>  result;
//       Multiply<Field, DIM+1, DIM+1 >((*this), S, result);
//       return result;
//     }
    
    /**
     * \brief  The multiplication operator Translation = CellRotation * CellDirection
     */
    template<typename In_Field>
    CellDirection<Field,DIM> operator*(const CellDirection<In_Field,DIM> t) {
      CellDirection<Field,DIM> result;
      Multiply((RotationClass) (*this), 
	       t, 
	       (RotationClass) result);
      return result;
    }
  
    /**
     *\brief  The multiplication operator for SeitzPosition = CellRotation * CellPosition
     */
    template<typename In_Field, typename Algorithms>
    CellPosition<Field,DIM,Algorithms>  operator * (const CellPosition<In_Field,DIM,Algorithms>& p) const {
      
      CellPosition<Field,DIM,Algorithms> result;
      Multiply((RotationClass)(*this), 
	       p, 
	       (RotationClass)result);
      return result;
    }

  }; /* class CellRotation */

  template<typename Field, size_t DIM>
  Field Det(const CellRotation<Field, DIM>& cellRotation) {
    return Det( (typename CellRotation<Field, DIM>::RotationClass) cellRotation);
  }

  template<typename Field, size_t DIM>
  Field Trace(const CellRotation<Field, DIM>& cellRotation) {
    return Trace( (typename CellRotation<Field, DIM>::RotationClass) cellRotation);
  }

  template<typename Field, size_t DIM>
  void Inverse(const CellRotation<Field, DIM>& cellRotation,
	       CellRotation<Field, DIM>& inCellRotation) {

    typedef typename CellRotation<Field, DIM>::RotationClass RC;

    const RC& m = cellRotation;
    RC&     inv = inCellRotation;
    Inverse(m, inv);
  }

} /* namespace psimag */

#endif //Psimag_CellRotation
/*@}*/
