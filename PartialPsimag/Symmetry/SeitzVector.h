//-*-C++-*-
#ifndef  Psimag_Seitz_Vector
#define  Psimag_Seitz_Vector

/** \ingroup extendedMatrices */
/*@{*/

/** \file SeitzVector.h
 *
 * Contains classes for implementing a seitz vector and the
 * outputstream operator for that class.
 */

#include "Vec.h"
#include <iostream>
#include "FieldConvert.h"
#include "Mat.h"

namespace psimag {
  
  /** \ingroup extendedMatrices
   *
   * \brief Base class for implementing seitz position and translation vectors.
   *
   * Such vectors have an extra component indicating whether they are
   * translation vectors ( component == 0) or position vectors (
   * component == 1).
   *
   * The dimension of such vectors is one greater than the
   * working/physical dimension.
   *
   * \note We keep the vector stuff private because some of the
   *       sematics is different. For example the negation operator.
   */
  template<typename Field, 
	   size_t   DIM,
	   int      IND=0,
	   typename TRAITS=ColMajorTraits<Field,DIM+1,1> > 
  class SeitzVector: protected Vec< Field, DIM+1, TRAITS > {

    typedef Vec<Field, DIM+1, TRAITS>            BaseType;
    typedef ColMajorTraits<Field,DIM,1>          LimitedType;
    typedef SeitzVector<Field, DIM, IND, TRAITS> ThisType;

  public:

    template<typename T> class TypeComputer {public: typedef T Result[DIM]; };

    typedef typename TypeComputer<Field>::Result ArrayType;
    typedef TRAITS                               Traits;

    /** 
     * The Default Constructor produces the zero translation. 
     */
    SeitzVector(): BaseType(Field(0)) 
    { (*this)[DIM] = Field(IND);  }
    
    /** 
     * Construct a translation whose components are all the given value. 
     *
     * Convert data types as necessary.
     */
    //template<typename IN_TYPE>
    //SeitzVector(const IN_TYPE& val): BaseType(convert<Field>(val)) 
    explicit SeitzVector(const Field& val): BaseType(val) 
    {
      (*this)[DIM] = Field(IND); 
    }
    
    /** 
     * Construct a vector whose components are from the given SeitzVector. 
     *
     * \note The vector may not have the same IND. If so converting
     *       from position to translation or visa versa
     *
     * \note Implicit Field type conversions made possible.
     */
    template<typename IN_TYPE, int otherIND, typename OtherTraits>
    explicit SeitzVector(const SeitzVector<IN_TYPE, DIM, otherIND, OtherTraits>& val)
    {
      COPY<ThisType, 
	SeitzVector<IN_TYPE, DIM, otherIND, OtherTraits>,
	Traits,
	OtherTraits >::EXEC(*this, val);
    }
    
    /** 
     * Construct a vector whose components are from the given value. 
     *
     * \note I tried templatizing this Vec -> VecType template arg but
     *       I had trouble with this template matching when I did not
     *       want it to. This can probably be done if this issue is
     *       worked thorugh, Later. MSS
     */
    template<typename IN_TYPE, typename VTraits>
    SeitzVector(const Vec<IN_TYPE, DIM+1, VTraits>& val)
    {
      typedef Vec<IN_TYPE, DIM+1, VTraits>  OtherType;

      COPY<ThisType, OtherType, LimitedType, LimitedType>::EXEC(*this, val);
      (*this)[DIM] = Field(IND); 
    }
    
    /** 
     * Construct a vector whose components are from the given value. 
     *
     * \note I tried templatizing this Vec -> VecType template arg but
     *       I had trouble with this template matching when I did not
     *       want it to. This can probably be done if this issue is
     *       worked thorugh, Later. MSS
     */
    template<typename IN_TYPE, typename VTraits>
    explicit
    SeitzVector(const Vec<IN_TYPE, DIM, VTraits>& val)
    {
      // Only work with the first three components of the SeitzVector
      typedef Vec<IN_TYPE, DIM, VTraits>  OtherType;
      COPY<ThisType, OtherType, LimitedType, VTraits >::EXEC(*this, val);
      (*this)[DIM] = Field(IND); 
    }
    
    /** 
     * Construct a translation whose components are set to the given value array. 
     *
     * Convert data types as necessary.
     */
    SeitzVector(const ArrayType&  vals) { 
      
      // Only work with the first three components of the SeitzVector
      typedef LimitedType ArrayTraits;
      
      COPY<ThisType, ArrayType, LimitedType, ArrayTraits>::EXEC(*this, vals);
      
      (*this)[DIM] = Field(IND);   // Defaults to a translation vector
    }

    //======================================================================
    /** Returns the L2Norm of this vector
     */
    Field l2norm ()  { 
      return L2NORM<ThisType, ThisType,LimitedType,LimitedType>::EXEC((*this),(*this));
  }

    // Inheriting and Wrapping BaseType operators

    /** 
     * The [] operator is inherited from BaseType. 
     */
    Field& operator [] (size_t i) {  
      return BaseType::operator[](i);
    }

    /** 
     * The [] operator is inherited from BaseType. 
     */
    const Field& operator [] (size_t i) const {  
      return BaseType::operator[](i);
    }

    /**
     * The assignment operator.
     *
     * \note Translations and Positions cannot be assigned to each other.
     */
    template<typename Field2, typename OtherTraits>
    ThisType& operator = (const SeitzVector<Field2,DIM,IND, OtherTraits>& v)
    { 
      typedef SeitzVector<Field2,DIM,IND, OtherTraits> OtherType;
      COPY<ThisType,OtherType,Traits,OtherTraits>::EXEC((*this),v);
      return *this; 
    }
  
    /**
     * The translation/position add translation into operator.
     *
     * \note Adding a position into a translation or position is not defined.
     */
    template<typename Field2, typename OtherTraits>
    ThisType& operator += (const SeitzVector<Field2,DIM,0, OtherTraits>& v)
    { 
      typedef SeitzVector<Field2,DIM,IND, OtherTraits> OtherType;
      COPY<ThisType,OtherType,Traits,OtherTraits>::EXEC_PLUS((*this),v);
      return *this; 
    }
  
    /**
     * The translation/position subtract into operator.
     *
     * \note Subtracting a position into a translation or position is not defined.
     */
    template<typename Field2, typename OtherTraits>
    SeitzVector<Field,DIM,IND>& operator -= (const SeitzVector<Field2,DIM,0>& v)
    { 
      typedef SeitzVector<Field2,DIM,IND, OtherTraits> OtherType;
      COPY<ThisType,OtherType,Traits,OtherTraits>::EXEC_MINUS((*this),v);
      return *this; 
    }
  

    //======================================================================
    /**
     * \brief Add a scalar value into this vector.
     *
     */
    template<typename Field2>
    ThisType& operator += (const Field2& val){
      // Replace with a loopless generic later.
      for(size_t i=0; i< DIM; i++)
	(*this)[i] += convert<Field>(val);
      return (*this);
    }

    //======================================================================
    /**
     * \brief Subtract a scalar value from this vector.
     *
     */
    template<typename Field2>
    ThisType& operator -= (const Field2& val){
      // Replace with a loopless generic later.
      for(size_t i=0; i< DIM; i++)
	(*this)[i] -= convert<Field>(val);
      return (*this);
    }

    //======================================================================
    /**
     * \brief Multiply a scalar value into this vector.
     *
     */
    template<typename Field2>
    ThisType& operator *= (const Field2& val){
      // Replace with a loopless generic later.
      for(size_t i=0; i< DIM; i++)
	(*this)[i] *= convert<Field>(val);
      return (*this);
    }

  };

  //======================================================================
  //======================================================================
  //======================================================================

  /**
   * \brief Multiply: MatType<T,DIM,DIM> * SeitzVector<..,0>
   *
   * \note we do not apply matrices to positions, just translations.
   */
  template<size_t DIM, typename VecField, typename MatField,  template<typename,size_t,size_t> class MatType>
  void Multiply(MatType<MatField,DIM,DIM> m,
		SeitzVector<VecField,DIM,0> v,
		SeitzVector<VecField,DIM,0> result) {
    
    typedef MatType<MatField,DIM,DIM>      MatrixType;
    typedef typename MatrixType::Traits    MatrixTraits;
    typedef SeitzVector<VecField,DIM,0>    VectorType;
    typedef ColMajorTraits<VecField,DIM,1> LimitedType;
    
    MatMultVec<MatrixType,MatrixTraits,VectorType,LimitedType,VectorType,LimitedType>::EXEC(m,v,result);
  }

  /**
   * \brief Compute the length of a SeitzTranslation vector
   *        represented in the given metric.
   *
   * That is, compute the length of a translation vector defined in
   * the coordinate system owhich has the given MetricTensor.
   *
   */
  template<size_t DIM, typename VecField, typename MatField,  template<typename,size_t,size_t> class MatType>
  MatField length(SeitzVector<VecField,DIM,0> v, MatType<MatField,DIM,DIM> metric) {

    typedef MatType<MatField,DIM,DIM>    MetricType;
    typedef typename MetricType::Traits  MetricTraits;
    typedef SeitzVector<VecField,DIM,0>  VectorType;
    
    VectorType  temp;
    Multiply(metric,v,temp);
    return v * temp;
  }

  /**
   * \brief The negative of a translation is the inverse translation.
   *
   * \note The negative of a position is not defined.
   */
  template<typename Field,size_t DIM>
  SeitzVector<Field,DIM,0> operator - (const SeitzVector<Field,DIM,0>& lhs){
    typedef SeitzVector<Field,DIM,0> VT;
    VT result;
    // Replace with a loopless generic later.
    for(size_t i=0; i< DIM; i++)
      result[i] *= Field(-1);
    return result;
  }

  //======================================================================
  /**
   * \brief The Seitz Translation|Position Equality operator
   *
   * \note Translations are not comparable to positions.
   */
  template<typename Field1,typename Field2,size_t DIM, int IND>
  bool operator == (const SeitzVector<Field1,DIM,IND>& lhs,
		    const SeitzVector<Field2,DIM,IND>& rhs){
    typedef SeitzVector<Field1,DIM,IND>   V1;
    typedef SeitzVector<Field2,DIM,IND>   V2;
    typedef ColMajorTraits<Field1,DIM,1>  V1T;
    typedef ColMajorTraits<Field2,DIM,1>  V2T;
    return EQUAL<V1,V2,V1T,V2T>::EXEC(lhs,rhs);
  }

  //======================================================================
  /**
   * \brief The sum of two translations is the composit translation.
   */
  template<typename Field,size_t DIM>
  SeitzVector<Field,DIM,0> operator + (const SeitzVector<Field,DIM,0>& lhs,
				       const SeitzVector<Field,DIM,0>& rhs){
    typedef SeitzVector<Field,DIM,0> VT;
    typedef typename SeitzVector<Field,DIM,0>::Traits TTraits;
    VT result;
    SUM<VT,VT,VT,TTraits,TTraits,TTraits>::EXEC(lhs,rhs,result);
    return result;
  }

  //======================================================================
  /**
   * \brief An - operator for Seitz translation vectors.
   *
   */
  template<typename Field,size_t DIM>
  SeitzVector<Field,DIM,0> operator - (const SeitzVector<Field,DIM,0>& lhs,
				       const SeitzVector<Field,DIM,0>& rhs){
    typedef SeitzVector<Field,DIM,0> VT;
    typedef typename SeitzVector<Field,DIM,0>::Traits TTraits;
    VT result;
    DIFFERENCE<VT,VT,VT,TTraits,TTraits,TTraits>::EXEC(lhs,rhs,result);
    return result;
  }

  //======================================================================
  /**
   * \brief The difference between two positions is the translation
   *        vector taking the rhs to the lhs.
   *
   */
  template<typename Field,size_t DIM>
  SeitzVector<Field,DIM,0> operator - (const SeitzVector<Field,DIM,1>& lhs,
				       const SeitzVector<Field,DIM,1>& rhs){
    typedef SeitzVector<Field,DIM,0> VT;
    typedef SeitzVector<Field,DIM,1> VP;
    typedef typename SeitzVector<Field,DIM,0>::Traits TTraits;
    typedef typename SeitzVector<Field,DIM,1>::Traits PTraits;
    VT result;
    DIFFERENCE<VP,VP,VT,PTraits,PTraits,TTraits>::EXEC(lhs,rhs,result);
    return result;
  }

  //======================================================================
  /**
   * \brief The sum of a translation and a position is the translated
   *        position.
   */
  template<typename Field,size_t DIM>
  SeitzVector<Field,DIM,1> operator + (const SeitzVector<Field,DIM,0>& lhs,
				       const SeitzVector<Field,DIM,1>& rhs){
    typedef SeitzVector<Field,DIM,0> VT;
    typedef SeitzVector<Field,DIM,1> VP;
    typedef typename SeitzVector<Field,DIM,0>::Traits TTraits;
    typedef typename SeitzVector<Field,DIM,1>::Traits PTraits;
    VP result;
    SUM<VT,VP,VP,TTraits,PTraits,PTraits>::EXEC(lhs,rhs,result);
    return result;
  }

  //======================================================================
  /**
   * \brief The inner product of two translations.
   */
  template<typename LField, typename RField, size_t DIM>
  LField operator * (const SeitzVector<LField,DIM,0>& lhs,
		    const SeitzVector<RField,DIM,0>& rhs){
    typedef SeitzVector<LField,DIM,0>    LVT;
    typedef SeitzVector<RField,DIM,0>    RVT;
    typedef Mat<LField,1,1>              RT;
    typedef ColMajorTraits<LField,1,DIM> VTL;
    typedef ColMajorTraits<RField,DIM,1> VTR;
    typedef typename RT::Traits          RTT;
    RT result;
    MatMultVec<LVT,VTL,RVT,VTR,RT,RTT>::EXEC(lhs,rhs,result);
    return result[0];
  }

  //======================================================================
  /**
   * \brief The vector time a scalor.
   */
  template<typename Field,size_t DIM, int IND>
  SeitzVector<Field,DIM,IND>  operator * (Field scalar, const SeitzVector<Field,DIM,IND>& lhs){
    return operator * (lhs,scalar);
  }

  //======================================================================
  /**
   * \brief The vector times a scalor.
   */
  template<typename Field,size_t DIM, int IND>
  SeitzVector<Field,DIM,IND>  operator * (const SeitzVector<Field,DIM,IND>& lhs, Field scalar){
    SeitzVector<Field,DIM,IND> result;
    for (size_t i=0; i< DIM;i++)
      result[i] = lhs[i] * scalar;
    return result;
  }

  //======================================================================
  /**
   * \brief The vector divided by a scalor.
   *
   * change later &*&*&*&* use Mat.foreach
   */
  template<typename Field,size_t DIM>
  SeitzVector<Field,DIM,0>  operator / (const SeitzVector<Field,DIM,0>& lhs, Field scalar){
    SeitzVector<Field,DIM,0> result;
    for (size_t i=0; i< DIM;i++)
      result[i] = lhs[i] / scalar;
    return result;
  }

//   //======================================================================
//   /**
//    * \brief The vector times a scalor.
//    */
//   template<typename Field,size_t DIM>
//   SeitzVector<Field,DIM,0>  operator * (const SeitzVector<Field,DIM,0>& lhs, char* scalarStr){
//     Field scalar = convert<Field>(scalarStr);
//     return operator * (lhs,scalar);
//   }

  //======================================================================
  /**
   * \brief The vector times a scalor.
   */
  template<typename Field,size_t DIM, int IND>
  SeitzVector<Field,DIM,IND>  operator * (const SeitzVector<Field,DIM,IND>& lhs, std::string scalarStr){
    Field scalar = convert<Field>(scalarStr);
    return operator * (lhs,scalar);
  }

  //======================================================================

  /** Returns vector product of a and b 
   */
  template < class Field, int IND>
  SeitzVector<Field,3,IND> operator % (const SeitzVector<Field,3,IND>& a, const SeitzVector<Field,3,IND>& b) 
  { 
    SeitzVector<Field,3,IND> result;
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
  }

  //======================================================================

  /** Returns vector product of a and b 
   */
  template < class Field, int IND>
  SeitzVector<Field,3,IND> operator % (const SeitzVector<Field,2,IND>& a, const SeitzVector<Field,2,IND>& b) 
  { 
    const static Field zero = Field(0);
    SeitzVector<Field,3,IND> result;
    result[0] = zero;
    result[1] = zero;
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
  }

  //======================================================================

  /** Returns vector product of a and b 
   */
  template <class Field, int IND>
  Field Magnitude (const SeitzVector<Field,2,IND>& v)
  { 
    return L2NORM<SeitzVector<Field,2,IND>, 
      SeitzVector<Field,2,IND>,
      ColMajorTraits<Field,2,1>,
      ColMajorTraits<Field,2,1> >::EXEC(v,v);
  }

  //======================================================================
  
  /** \ingroup ostream
   *
   * SeitzVector output stream operator 
   **/
  template<typename Field,size_t DIM, int IND>
  std::ostream& operator << (std::ostream& os, const SeitzVector<Field,DIM,IND>& t) {
    
    bool flat = true;
    MAT_PRINT<SeitzVector<Field,DIM,IND>,ColMajorTraits<Field,DIM,1> >::EXEC(t, os, flat);
    if (IND==1)
      os << "POSIT";
    else
      os << "TRANS";
    return os;
  }

  //======================================================================
}  /* namspace psimag */

#endif // Psimag_Seitz_Vector
/*@}*/
