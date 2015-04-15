//-*-C++-*-

/** \ingroup Data_Structures */
/*@{*/

/*! \file RealPart.h  
 *
 */

#ifndef PSIMAG_RealPart_H
#define PSIMAG_RealPart_H

namespace psimag {
  
  //====================================================================== Row Slice
  
  template<typename ComplexMatrixLikeType>
  class RealPart {
  public:
    
    typedef typename ComplexMatrixLikeType::value_type              complex_value_type;
    typedef  typename ComplexMatrixLikeType::value_type::value_type value_type;
    typedef value_type                                              FieldType;
    
    ComplexMatrixLikeType&       mat;
    
    RealPart(ComplexMatrixLikeType& m):
      mat(m)
    {}
    
    // ----- MatrixLike
    
    inline size_t n_col() const { return mat.n_col(); }
    inline size_t n_row() const { return mat.n_row(); }
    
    value_type& operator () (size_t rowIndex, size_t colIndex) {
      return real(mat(rowIndex,colIndex));
    }
    
    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      return real(mat(rowIndex,colIndex));
    }

    // ----- VectorLike

    size_t size()  const { return mat.size(); }
    
    const value_type& operator[] (size_t componentIndex) const {
      return real(mat[componentIndex]);
    }
    
    inline
    value_type operator[] (size_t componentIndex)  {
      return real(mat[componentIndex]);
    }
  };
  
} // end namespace PSIMAG


/*@}*/
#endif
