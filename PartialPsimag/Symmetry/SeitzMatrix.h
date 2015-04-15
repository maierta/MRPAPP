//-*-C++-*-

#ifndef  Psimag_SeitzMatrix
#define  Psimag_SeitzMatrix

/** \ingroup extendedMatrices */
/*@{*/

/** \file SeitzMatrix.h
 *
 *  Contains a class for Matricies that represent Symmetry Operations, et. al.
 */
 
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

#include "Vec.h"
#include "Mat.h"
#include "Tag.h"

#include "PSIMAGAssert.h"
#include "SeitzVector.h"
#include "SeitzMatrixTraits.h"

namespace psimag {


  /** \ingroup extendedMatrices 
   *
   *\brief A template class representng affine matricies and is used
   *       for representing symmetry operations.  (See 2002 ICT
   *       Section 5.2)
   *
   *\note These matricies are also refered to as extended
   *      matricies. If the physical dimension is DIM the matrix is of
   *      dimension DIM+1.
   */
  template<typename Field, size_t DIM, typename RotationType_=Mat<Field,DIM,DIM> >
  class SeitzMatrix
  {

  public:
    
    //STATIC_CHECK(DIM==DIM2);

    template<typename IN_TYPE>
    class TypeComputer {
    public:
      typedef IN_TYPE DataArray2D[DIM+1][DIM+1];
      typedef IN_TYPE RotArray2D[DIM][DIM];
      typedef IN_TYPE DataArray[DIM+1*DIM+1];
    };

    typedef SeitzMatrix<Field,DIM,RotationType_>  ThisType;
    typedef SeitzMatrixTraits<Field, DIM>         Traits;
    typedef RotationType_                         RotationType;
    typedef typename RotationType::Traits         RotationTraits;
    typedef Vec<Field, DIM>                       TranslationType;
    typedef Vec<Field, DIM+1>                     LastRowType;

    //    template<size_t ROW, size_t COL>
    //friend class SeitzMatrixTraits<Field, DIM>::REF<ThisType,ROW,COL>::RotRetriever;

    RotationType             rotation;
    TranslationType          translation;
    LastRowType              lastRow;
    
    RotationType             invRotation;
    RotationType             transRotation;
    TranslationType          invTranslation;
    bool                     invSet;
    bool                     transSet;
    
  public:
    
    // ================================================== Essential Operator
    
    /**
     * \brief Array access operator.
     *
     * Handle the 'extended' part and defer to the rotation and
     * translation mebers for the rest.
     *
     * \note For operations like =, ==, *, etc which are implemented
     *       (directly or indirectly) on top of algorithms such as
     *       MatForEach and MatReduce this operation will not be
     *       called. These algorithms access the proper value at
     *       compile time via SeitzMatrizTraits.
     */
    Field& operator () (const size_t& row, const size_t& col) {

      // handle the extended part
      if (row == DIM) {
	return lastRow[col];
      }
      // The value is in the translation part
      if ( col == DIM )
	return translation[row];

      // The value is in the rotation part.
      return rotation(row,col);
    }
    
    /**
     * \brief Array access operator.
     *
     * Handle the 'extended' part and defer to the rotation and
     * translation mebers for the rest.
     *
     * \note For operations like =, ==, *, etc which are implemented
     *       (directly or indirectly) on top of algorithms such as
     *       MatForEach and MatReduce this operation will not be
     *       called. These algorithms access the proper value at
     *       compile time via SeitzMatrizTraits.
     */
    const Field& operator () (const size_t& row, const size_t& col) const {

      // handle the extended part
      if (row == DIM) {
	return lastRow[col];
      }
      // The value is in the translation part
      if ( col == DIM )
	return translation[row];

      // The value is in the rotation part.
      return rotation(row,col);
    }
    
    // ================================================== Constructors
    
    /** 
     * The default constructor results in an identity Matrix. 
     *
     **/
    SeitzMatrix() : 
      rotation(Field(0)),
      translation(Field(0)),
      lastRow(Field(0)),
      invSet(false),               
      transSet(false)           
    {
      MakeIDENTITY<RotationType, RotationTraits>::EXEC(rotation);
      lastRow[DIM] = Field(1);
    }

     /** 
     * The copy constructor. 
     *
     **/
    SeitzMatrix(const ThisType& other) : 
      rotation(other.rotation),
      translation(other.translation),
      lastRow(Field(0)),
      invRotation(other.invRotation),
      transRotation(other.transRotation),
      invTranslation(other.invTranslation),
      invSet(other.invSet),               
      transSet(other.transSet)              
    {
      lastRow[DIM] = Field(1);
    }

   /** 
     * \breif Construct an extended matrix with the cooficients of
     *        it's non-extended part set to the given Field value.
     *
     **/
    SeitzMatrix(const Field& a): 
      rotation(Field(a)),
      translation(Field(a)),
      lastRow(Field(0)),
      invSet(false),        
      transSet(false)  
    { lastRow[DIM] = Field(1); }

//     /** 
//      * Construct an extended matrix from a by DIMxDIM+1 data array.
//      *
//      * The non-extended and translation parts are set by the given
//      * data. Data types are converted as necessary.
//      *
//      **/
//     template<typename IN_TYPE>
//     SeitzMatrix(const IN_TYPE data [DIM+1][DIM+1]): 
//       invSet(false),  
//       lastRow(0)
//     {
//       lastRow[DIM] = Field(1);
//       COPY<typename TypeComputer<IN_TYPE>::DataArray2D, ThisType, RowMajorTraits<IN_TYPE, DIM+1>, Traits >(data, (*this));
//     }

//     /** 
//      * Construct an extended matrix from a by DIMxDIM data array.
//      *
//      * The non-extended part is set to the given data. Data types are
//      * converted as necessary.
//      *
//      **/
//     template<typename IN_TYPE>
//     SeitzMatrix(const IN_TYPE data [DIM][DIM], SeitzTranslation<IN_TYPE, DIM>& t): 
//       invSet(false),  
//       rotation(data),
//       translation(t),
//       lastRow(0)
//     {
//       lastRow[DIM] = Field(1);
//       // this->rotation = data;
//       //COPY<TypeComputer<In_TYPE>::RotArray2D, RotationType, RowMajorTraits<IN_TYPE, DIM>, RotationTraits >(data, This->rotation);
//       //setTranslation(t);
//     }

    /** 
     * Construct an extended matrix from a by DIMxDIM data array.
     *
     * The non-extended part is set to the given data. Data types are
     * converted as necessary.
     *
     **/
    template<typename IN_TYPE>
    SeitzMatrix(const IN_TYPE data [DIM-1*DIM-1] ): 
      rotation(data),
      translation(Field(0)),
      lastRow(Field(0)),
      invSet(false),  
      transSet(false)
    {
      lastRow[DIM] = Field(1);
      // this->rotation = data;
      //COPY<TypeComputer<In_TYPE>::RotArray2D, RotationType, RowMajorTraits<IN_TYPE, DIM>, RotationTraits >(data, This->rotation);
      //setTranslation(t);
    }

    /** 
     * Construct an extended matrix from a DIM*DIM single dimension data array.
     *
     * The non-extended part of the matrix is determined by the
     * data. Data types are converted as necessary.
     *
     **/
    template<typename IN_TYPE>
    SeitzMatrix(const IN_TYPE data [(DIM-1)*(DIM-1)], IN_TYPE tdata [DIM-1]) : 
      rotation(data),
      translation(tdata),
      invSet(false),  
      transSet(false)
    {
      lastRow[DIM] = Field(1);
    }

    /** 
     * Construct an extended matrix from a DIM*DIM single dimension data array.
     *
     * The non-extended part of the matrix is determined by the
     * data. Data types are converted as necessary.
     *
     **/
    template<typename IN_TYPE>
    SeitzMatrix(const IN_TYPE* data, const IN_TYPE* tdata) : 
      rotation(data),
      translation(tdata),
      invSet(false),  
      transSet(false)
    {
      lastRow[DIM] = Field(1);
    }

//     /** 
//      * Construct an extended matrix from a given translation. 
//      *
//      * The non-extended part is set to the
//      * identity and the rightmost column is set to the
//      * SeitzTranslation. Convert types as necessary.
//      *
//      **/
//     template<typename IN_TYPE, 
// 	     template<typename, size_t> class VecType> 
//     SeitzMatrix(const VecType<IN_TYPE, DIM>& t) : 
//       invSet(false),
//       rotation(Field(0)),
//       translation(Field(0)),
//       lastRow(Field(0))
//     {
//       MakeIDENTITY<RotationType, RotationTraits>::EXEC(rotation);
//       lastRow[DIM] = Field(1);
//       setTranslation(t);
//     }

    /** 
     * Construct an extended matrix from a given rotation and translation. 
     *
     * The non-extended part is set to the
     * identity and the rightmost column is set to the
     * SeitzTranslation. Convert types as necessary.
     *
     **/
    SeitzMatrix(const RotationType&    m,
		const TranslationType& t) : 
      rotation(m),
      lastRow(Field(0)), 
      invSet(false),
      transSet(false)
    {
      setTranslation(t);
      lastRow[DIM] = Field(1);
    }

    /** 
     * Construct an extended matrix from a given rotation and translation. 
     *
     * The non-extended part is set to the
     * identity and the rightmost column is set to the
     * SeitzTranslation. Convert types as necessary.
     *
     **/
    SeitzMatrix(const TranslationType& t) : 
      rotation(Field(0)),
      lastRow(Field(0)),
      invSet(false),
      transSet(false)
    {
      MakeIDENTITY<RotationType, RotationTraits>::EXEC(rotation);
      setTranslation(t);
      lastRow[DIM] = Field(1);
    }

    /** 
     * Construct an extended matrix from a given rotation and translation. 
     *
     * The non-extended part is set to the
     * identity and the rightmost column is set to the
     * SeitzTranslation. Convert types as necessary.
     *
     **/
    SeitzMatrix(const SeitzVector<Field,DIM,0>& t) : 
      rotation(Field(0)),
      lastRow(Field(0)), 
      invSet(false),
      transSet(false)
    {
      MakeIDENTITY<RotationType, RotationTraits>::EXEC(rotation);
      setTranslation(t);
      lastRow[DIM] = Field(1);
    }

    /** 
     * Construct an extended matrix from a given rotation and translation. 
     *
     * The non-extended part is set to the
     * identity and the rightmost column is set to the
     * SeitzTranslation. Convert types as necessary.
     *
     **/
    SeitzMatrix(const RotationType&    m) : 
      rotation(m),
      translation(Field(0)),
      lastRow(Field(0)),
      invSet(false),
      transSet(false)
    {
      lastRow[DIM] = Field(1);
    }

    /** 
     * Construct an extended matrix from a given rotation and translation. 
     *
     * The non-extended part is set to the
     * identity and the rightmost column is set to the
     * SeitzTranslation. Convert types as necessary.
     *
     **/
    template<typename IN_TYPE, 
	     template<typename, size_t> class MatType,
	     template<typename, size_t> class VecType> 
    SeitzMatrix(const MatType<IN_TYPE, DIM>& m,
		const VecType<IN_TYPE, DIM>& t) : 
      lastRow(Field(0)),
      invSet(false),
      transSet(false)
      //rotation(m),
      //translation(0),
    {
      //setTranslation(t);
      lastRow[DIM] = Field(1);
    }

     // ================================================== Utility functions
    
    // ================================================== get Translation

    /** 
     * Return a Vec containing the translation components of this SeitzMatrix.
     */
    const TranslationType& getTranslation() {
      return translation;
    }

    // ================================================== set Translation

    /** 
     * Set the translation components of this SeitzMatrix according to the values of the given SeitzTranslation.
     *
     * Convert types as necessary.
     */
    template<typename IN_TYPE, 
	     template<typename, size_t, size_t>   class TraitsTemplate,
	     template<typename, size_t, typename> class VecTemplate> 

    void setTranslation(const VecTemplate<IN_TYPE, DIM, TraitsTemplate<IN_TYPE, DIM,1> >& t) {
      COPY<TranslationType, 
	VecTemplate<IN_TYPE, DIM, TraitsTemplate<IN_TYPE, DIM,1> >, 
	typename TranslationType::Traits,
	TraitsTemplate<IN_TYPE, DIM,1> >::EXEC(translation, t);
      //for(size_t i=0; i< DIM; i++)
      //	this->translation[i] = convert<Field>(t[i]);
    }

    /** 
     * Set the translation components of this SeitzMatrix according to the values of the given SeitzTranslation.
     *
     * Convert types as necessary.
     */
    template<typename IN_TYPE, 
	     template<typename, size_t, size_t>   class TraitsTemplate> 
    void setTranslation(const SeitzVector<IN_TYPE, DIM, 0, TraitsTemplate<IN_TYPE, DIM+1, 1> >& t) {
      COPY<TranslationType, 
	SeitzVector<IN_TYPE, DIM, 0, TraitsTemplate<IN_TYPE, DIM+1, 1> >, 
	typename TranslationType::Traits,
	TraitsTemplate<IN_TYPE, DIM, 1> >::EXEC(translation, t);
      //for(size_t i=0; i< DIM; i++)
      //	this->translation[i] = convert<Field>(t[i]);
    }

    /** 
     * Set the translation components of this SeitzMatrix to the given value.
     *
     * Convert types as necessary.
     */
    template<typename IN_TYPE>
    void setTranslationValue(const IN_TYPE& value) {
      for(size_t i=0; i< DIM; i++)
	this->translation[i] = convert<Field>(value);
    }

    /** 
     * Set the translation components of this SeitzMatrix according to the given array of values.
     *
     * Convert types as necessary.
     */
    template<typename IN_TYPE>
    void setTranslation(const IN_TYPE values[DIM]) {
      for(size_t i=0; i< DIM; i++)
	this->translation[i] = convert<Field>(values[i]);
    }

    /** 
     * Set the translation components of this SeitzMatrix according to the values of the given Vec.
     *
     * Convert types as necessary.
     */
    template<typename IN_TYPE>
    void setTranslation(const Vec<IN_TYPE, DIM> t) {
      for(size_t i=0; i< DIM; i++)
	this->translation[i] = convert<Field>(t[i]);
    }

    // ================================================== RotationMat

    // ================================================== Determinant
    /** 
     * Return the determinan of the rotation part.
     */
    Field determinant() const {
      return Det(rotation);
    }

    // ================================================== Magnitude
    /** 
     * Return the magnitude of the rotation and translation parts.
     */
    Field magnitude() const {
      Field rMag = Magnitude(this->rotation);
      Field tMag = L1Norm(this->translation);
      return  max(rMag,tMag);
    }

    // ================================================== Transpose
    /** 
     * Compute the transpose
     */
    const SeitzMatrix<Field, DIM> transpose()  {

      if (transSet) {
	SeitzMatrix<Field, DIM> result(transRotation, translation);
	return result;
      }
      
      Transpose(rotation, transRotation);
      SeitzMatrix<Field, DIM> result(transRotation, translation);
      transSet = true;
      return result;
    }

    // ================================================== Const Transpose
    /** 
     * Compute the transpose
     */
    const SeitzMatrix<Field, DIM> transpose()  const{

      if (transSet) {
	SeitzMatrix<Field, DIM> result(transRotation, translation);
	return result;
      }
      
      RotationType             transRotation_;
      Transpose(rotation, transRotation_);
      SeitzMatrix<Field, DIM> result(transRotation_, translation);
      return result;
    }

   // ================================================== Inverse
    /** 
     * Compute the inverse
     */
    const SeitzMatrix<Field, DIM> inverse()  {

      if (invSet) {
	SeitzMatrix<Field, DIM> result(invRotation, invTranslation);
	return result;
      }
      
      Inverse(rotation, invRotation);
      Multiply(invRotation, translation, invTranslation);
      //Multiply<Field,DIM,DIM,DIM, RotationType, Mat, CMT, Mat, CMT,Mat>  (invRotation, translation, invTranslation);

      invTranslation *= Field(-1);
      invSet = true;

      SeitzMatrix<Field, DIM> result(invRotation, invTranslation);
      return result;
    }

    // ================================================== Const Inverse
    /** 
     * Compute the inverse const ly
     */
    const SeitzMatrix<Field, DIM> inverse() const {

      if (invSet) {
	SeitzMatrix<Field, DIM> result(invRotation, invTranslation);
	return result;
      }
      
      RotationType             invRotation_;
      TranslationType          invTranslation_;

      Inverse(rotation, invRotation_);
      Multiply(invRotation_, translation, invTranslation_);
      invTranslation_ *= Field(-1);

      SeitzMatrix<Field, DIM> result(invRotation_, invTranslation_);
      return result;
    }

    // ================================================== Const Inverse
    /** 
     * Compute the inverse const ly
     */
    const SeitzMatrix<Field, DIM> constInverse() const {

      if (invSet) {
	SeitzMatrix<Field, DIM> result(invRotation, invTranslation);
	return result;
      }
      
      RotationType             invRotation_;
      TranslationType          invTranslation_;

      Inverse(rotation, invRotation_);
      Multiply(invRotation_, translation, invTranslation_);
      invTranslation_ *= Field(-1);

      SeitzMatrix<Field, DIM> result(invRotation_, invTranslation_);
      return result;
    }

    // ================================================== Operators

    /**
     * \brief  The operator for multiplying a SeitzMatrix times another matrix
     */
    SeitzMatrix<Field,DIM> operator*(const SeitzMatrix<Field,DIM>& S) {
      SeitzMatrix<Field,DIM>  result;
      Multiply((*this), S, result);
      return result;
    }
    
    /**
     * \brief  The multiplication operator SeitzVector = SeitzMatrix * SeitzVector
     */
    template<typename Trans_Field>
    SeitzVector<Field,DIM> operator*(const SeitzVector<Trans_Field,DIM> t) {
      SeitzVector<Field,DIM> result;
      // This statement is only needed if a field conversion is needed add a specialization?
      SeitzVector<Field,DIM> temp(t); 
      Multiply((*this), temp, result);
      return result;
    }
  
    /**
     * \brief  Add a translation into this SeitzMatrix 
     */
    template<typename IN_TYPE>
    SeitzMatrix<Field,DIM> operator += (const SeitzVector<IN_TYPE,DIM>& T) {
      SeitzMatrix<Field,DIM> result((*this));
      for(size_t i=0; i<DIM; ++i) 
	result(i,DIM) += convert<Field>(T[i]);
      return result;
    }

    /**
     * \brief  Swap the rows of the rotation 
     */
    void swapRows(size_t i, size_t j) { rotation.swapRows(i,j);}
    
    /**
     * \brief  Swap the cols of the rotation 
     */
    void swapCols(size_t i, size_t j) { rotation.swapCols(i,j);}
    
    /**
     * \brief  Negate the row of the rotation 
     */
    void negateRow(size_t i) { rotation.negateRow(i);}
    
    /**
     * \brief  Negate the cols of the rotation 
     */
    void negateCol(size_t i) { rotation.negateCol(i);}


  }; /* class SeitzMatrix: public Mat< Field, DIM+1, DIM+1 > */


  // ================================================== Algorithms

  /**
   * \brief Algorithm to calculate the inverse of a Seitz matrix.
   *
   * see 2002 ITC, page 84 inverse = (rot_inv, -(rot_inv)t)
   */
  template<class Field, size_t DIM>
  void Inverse(SeitzMatrix<Field, DIM>& m, SeitzMatrix<Field, DIM>& inv) {
    inv = m.inverse();
  }
  
  /**
   * \brief Algorithm to calculate the Determinant of a Seitz matrix.
   *
   * see 2002 ITC, page 84 inverse = (rot_inv, -(rot_inv)t)
   */
  template<class Field, size_t DIM>
inline  Field Det(const SeitzMatrix<Field, DIM>& m) {
    return Det(m.rotation);
  }
  
  /**
   * \brief Algorithm to calculate the Trace of a Seitz matrix.
   *
   * see 2002 ITC, page 84 inverse = (rot_inv, -(rot_inv)t)
   */
  template<class Field, size_t DIM>
inline  Field Trace(const SeitzMatrix<Field, DIM>& m) {
    return Trace(m.rotation) + 1;
  }
  
  /** 
   *
   *  Equality operator
   */
  template<typename Field, size_t DIM>
  bool operator == (const SeitzMatrix<Field,DIM>& lhs, 
		    const SeitzMatrix<Field,DIM>& rhs) {
    return (lhs.rotation == rhs.rotation) && 
      (lhs.translation == rhs.translation);
  }
  
//   /** 
//    *
//    *  Equality operator
//    */
//   template<typename Field, size_t DIM, typename ValType>
//   bool operator == (const SeitzMatrix<Field,DIM>& lhs, 
// 		    const ValType&                val) {
//     //return (lhs.rotation == val) && (lhs.translation == val);

//     // Look like a 3x4 matrix instead of a 4x4 matrix
//     typedef SeitzMatrixTraits<Field,DIM,DIM,DIM+1>      SMT;

//     return EQUAL_VAL<SeitzMatrix<Field,DIM>, SMT, ValType>::EXEC(lhs,val);
//   }

  /**  Multiplication of a SeitzMatrix by SeitzVector */
  template< class T, size_t DIM, int IND>
  void Multiply(const SeitzMatrix<T,DIM>&            m,
		const SeitzVector<T,DIM,IND>&  v,
		SeitzVector<T,DIM,IND>&        result) { 
    
    MatMultVec<SeitzMatrix<T,DIM>, 
      typename SeitzMatrix<T,DIM>::Traits,
      SeitzVector<T,DIM,IND>, 
      typename SeitzVector<T,DIM,IND>::Traits,
      SeitzVector<T,DIM,IND>, 
      typename SeitzVector<T,DIM,IND>::Traits >::EXEC(m,v,result);
  }

  /** 
   * Multiplication of a SeitzMatrix by a SeitzMatrix
   *  
   */
  template<typename Field, size_t DIM>
  void Multiply(const SeitzMatrix<Field,DIM>& lhs,
		const SeitzMatrix<Field,DIM>& rhs,
		SeitzMatrix<Field,DIM>&       result) { 
    
    MatMultVec<SeitzMatrix<Field,DIM>, 
      typename SeitzMatrix<Field,DIM>::Traits, 
      SeitzMatrix<Field,DIM>, 
      typename SeitzMatrix<Field,DIM>::Traits,
      SeitzMatrix<Field,DIM>, 
      typename SeitzMatrix<Field,DIM>::Traits >::EXEC(lhs,rhs,result);
  }

  //======================================================================

  /** \ingroup ostream
   *
   * Output stream operator 
   */
  template<typename Field, size_t DIM>
  SeitzMatrix<Field,DIM> operator - (const SeitzMatrix<Field,DIM>& lhs,
				     const SeitzMatrix<Field,DIM>& rhs) {
    
    typedef typename SeitzMatrix<Field,DIM>::RotationType RT;
    typedef typename RT::Traits                           RTT;
    typedef typename SeitzMatrix<Field,DIM>::TranslationType TT;
    typedef typename TT::Traits                           TTT;
    
    SeitzMatrix<Field,DIM> result;
    DIFFERENCE<RT,RT,RT,RTT,RTT,RTT>::EXEC(lhs.rotation, rhs.rotation, result.rotation);
    DIFFERENCE<TT,TT,TT,TTT,TTT,TTT>::EXEC(lhs.translation, rhs.translation, result.translation);

    return result;
  }
  
  /** 
   * Maximum value of a SeitzMatrix
   *  
   */
  template< class Field, size_t DIM>
  Field Max(const SeitzMatrix<Field,DIM>& m) {
    return 
      MAX<SeitzMatrix<Field,DIM>,
      SeitzMatrix<Field,DIM>,
      SeitzMatrixTraits<Field,DIM,DIM,DIM+1>,
      SeitzMatrixTraits<Field,DIM,DIM,DIM+1> >::EXEC(m,m);

  //      max(Max(m.rotation), Max(m.translation));
  }

  //======================================================================
  //======================================================================

  /** \ingroup xml
   *
   * XML Output function for 2D CellParameters.
   */
  template<typename Field, size_t DIM>
  Tag toXML(const SeitzMatrix<Field,DIM>& m,
	    std::string name="SeitzMatrix") {
      
    typedef SeitzMatrix<Field,DIM>   MatType; 
    typedef typename MatType::Traits MatTraitsType;

    Tag result(name);
    std::ostringstream buff;
    MAT_PRINT<MatType,MatTraitsType>::JUST_NUMBERS( m, buff);
    result.content << buff.str();
    
    return result;
  }
  
  //======================================================================

  /** \ingroup ostream
   *
   * Output stream operator 
   */
  template<typename Field, size_t DIM>
  inline  std::ostream& operator << (std::ostream& os, const SeitzMatrix<Field,DIM>& s) {
    
      MAT_PRINT<SeitzMatrix<Field,DIM>, typename SeitzMatrix<Field,DIM>::Traits>::EXEC(s, os);

    return os;
  }
  
} /* namespace psimag */

#endif // Psimag_SeitzMatrix
/*@}*/
