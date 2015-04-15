//-*-C++-*-

#ifndef  Psimag_Cartesian_CartesianRotation
#define  Psimag_Cartesian_Rotation

/** \ingroup extendedMatrices */
/*@{*/

/** \file CartesianRotation.h
 *
 *  Contains a class that represent rotations in cartesian space.
 *
 */
 
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

#include "SeitzVectors.h"
#include "CellDirection.h"
#include "Mat.h"
#include "PSIMAGAssert.h"

namespace psimag {
  
  /** \ingroup extendedMatrices 
   *
   *\brief A subclass of the matrix, Mat<Field, DIM, DIM>.
   *       which is used to represent a rotation is cartesian space.
   *
   * &*&*&* This class is not yet elaborated. 
   *
   */
  template<typename Field, size_t DIM>
  class CartesianRotation: public Mat< Field, DIM, DIM > {
    
  private:
    
    CartesianRotation<Field, DIM>* inv;
    bool        invSet;

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
      
      int det = Det((*this));
      if ( det != 1 and det != -1 ) {
	std::string message("ERROR CartesianRotation: ");
	message += " Tried to construct a CartesianRotation with a Det != +/- 1.";
	throw std::range_error(message);
      }
      
      int tr  = Trace((*this));
      if ( tr < -3 or tr > 3 ) {
	std::string message("ERROR SymmetryOperation: ");
	message += " Tried to construct a CartesianRotation with a Trace out of range.";
	throw std::range_error(message);
      }
    }

  public:
    
    /** The destructor handles deleting 'embedded' inverses. **/
    ~CartesianRotation() {
      if (invSet)
	delete inv;
    }

    // ================================================== Constructors
    
    /** 
     * The default constructor results in an identity Matrix. 
     *
     **/
    CartesianRotation() : 
      invSet(false),               
      Mat<Field, DIM, DIM > ( 0 ) {
      setToIdentity();
    }

    /** 
     * Construct from a DIM by DIM int array.
     *
     **/
    template<typename IN_TYPE>
    CartesianRotation(IN_TYPE data [DIM][DIM]) : invSet(false) {
      
      size_t index = 0;
      for (size_t i=0; i<DIM; i++) 
	for(size_t j=0; j<DIM; j++) 
	  (*this)(i,j) = convert<Field>(data[i][j]);
      check();
    }

    /** 
     * Construct from a DIM*DIM single dimension int array.
     *
     **/
    template<typename IN_TYPE>
    CartesianRotation(IN_TYPE data [DIM*DIM]) : invSet(false) {
      
      size_t index = 0;
      for (size_t i=0; i<DIM; i++) 
	for(size_t j=0; j<DIM; j++) 
	  (*this)(i,j) = convert<Field>(data[index++]);
      check();
    }

    /**
     * Set the diagonal of the non-extended part of the matrix to Field(1) .
     */
    void setToIdentity() {
      
      for (size_t i=0; i<DIM; i++) 
	for(size_t j=0; j<DIM; j++) 
	  if (i==j)
	    (*this)(i,j) = Field(1);
	  else
	    (*this)(i,j) = Field(0);
    }

    // ================================================== 
    /** 
     * Return the determinan of the rotation part.
     */
    const Field determinant() {
      return Det(*this);
    }

    // ================================================== Inverse
    /** 
     * Return by lazy evaluation a reference to the inverse of this Rotation.
     */
    const CartesianRotation<Field, DIM>& inverse() {

      if (invSet) 
	return *inv;

      inv    = new CartesianRotation<Field, DIM>();
      invSet = true;

      Inverse<Field>(*this, inv);

      return *inv;
    }

    // ================================================== Operators

    /**
     * \brief  The operator for multiplying two CartesianMatrices
     */
    CartesianRotation<Field,DIM> operator*(const CartesianRotation<Field,DIM>& S) {
      CartesianRotation<Field,DIM>  result;
      Multiply<Field, DIM+1, DIM+1 >((*this), S, result);
      return result;
    }
    
    /**
     * \brief  The multiplication operator Translation = CartesianRotation * CartesianTranslation
     */
    template<typename Trans_Field>
    CartesianTranslation<Field,DIM> operator*(const CartesianTranslation<Trans_Field,DIM> t) {
      CartesianTranslation<Field,DIM> result;
      // This statement is only needed if a field conversion is needed add a specialization?
      CartesianTranslation<Field,DIM> temp(t); 
      Multiply<Field,DIM+1,DIM+1>((*this), temp, result);
      return result;
    }
  
    /**
     *\brief  The multiplication operator for CartesianPosition = CartesianRotation * CartesianPosition
     */
    CartesianPosition<Field,DIM>  operator * (const CartesianPosition<Field,DIM>& p) const {
      
      CartesianPosition<Field,DIM> result;
      Multiply<Field,DIM,DIM>((*this), p, result);
      return result;
    }

    /**
     * \brief  Add a translation into this CartesianRotation 
     */
    template<typename IN_TYPE>
    CartesianRotation<Field,DIM> operator += (const CartesianTranslation<IN_TYPE,DIM>& T) {
      CartesianRotation<Field,DIM> result((*this));
      for(size_t i=0; i<DIM; ++i) 
	result(i,DIM) += convert<Field>(T[i]);
      return result;
    }
    
  }; /* class CartesianRotation: public Mat< Field, DIM+1, DIM+1 > */

} /* namespace psimag */

#endif //Psimag_Cartesian_Rotation
/*@}*/
