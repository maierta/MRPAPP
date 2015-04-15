//-*-C++-*-

/** \ingroup Data_Structures */
/*@{*/

/*! \file ValueTypeFunction.h  
 *
 */

#ifndef PSIMAG_ValueTypeFunction_H
#define PSIMAG_ValueTypeFunction_H

namespace psimag {
  
  // The problem addressed is that:
  // (const Matrix<double>)::valueType is double and not const double
  template <typename MatrixOrVectorLikeType>
  class ValueTypeFunction {
  public:
    typedef typename MatrixOrVectorLikeType::value_type ValueType;
    typedef typename MatrixOrVectorLikeType::value_type MaybeConstValueType;
  };
  
  template <typename MatrixOrVectorLikeType>
    class ValueTypeFunction<const MatrixOrVectorLikeType> {
    public:
    typedef typename MatrixOrVectorLikeType::value_type        ValueType;
    typedef const typename MatrixOrVectorLikeType::value_type  MaybeConstValueType;
  };


} // end namespace PSIMAG


/*@}*/
#endif
