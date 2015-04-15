//-*-C++-*-

#ifndef PSIMAG_Mat_H
#define PSIMAG_Mat_H

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"
#include "Vec.h"

#include "MatCopy.h"
#include "MatDet.h"
#include "MatDifference.h"
#include "MatEqual.h"
#include "MatIdentity.h"
#include "MatInverse.h"
#include "MatMagnitude.h"
#include "MatMax.h"
#include "MatMult.h"
#include "MatPrint.h"
#include "MatSum.h"
#include "MatTrace.h"
#include "MatTraits.h"
#include "MatTranspose.h"

namespace psimag {

  //===================================================================// 
  //==================== Class Mat<TYPE,NROW,NCOL> ====================// 
  //===================================================================// 
  
  /** \brief A simple, fix-sized matrix. 
   *
   * \param Field:    The type of the element to keep in the matrix.
   * \param NROW: the number of rows in the matrix.
   * \param NCOL: the number of columns in the matrix.
   *
   * \note This class is intended to be used with Vec<Field,NROW>
   *
   * \author Thomas Schulthess, Michael Summers
   */
  template <typename Field, size_t NROW, size_t NCOL, 
	    typename TRAITS=ColMajorTraits<Field,NROW,NCOL> > 
  class Mat {

  private:

    enum { LENGTH=NCOL*NROW };
    typedef Field LocalType[LENGTH];
    LocalType d;

    //    Field d[NCOL*NROW];
    
  public:

    typedef TRAITS                           Traits;
    typedef ConstantTraits<Field,NROW,NCOL>  ConstantMatTraitsType;
    typedef Mat<Field,NROW,NCOL,TRAITS>      ThisType;

    /** 
     * \brief Default constructor, all elements set to zero.
     */
    Mat() { 
      static Field zero(0);
      COPY<LocalType, Field, Traits, ConstantMatTraitsType>::EXEC(d, zero);
    }
    
    /** 
     * \brief Simple constructor, all elements set to a.
     *
     * \note This constructor will not be available for use in type
     *       conversions.
     */
    explicit Mat(const Field& a) { 
      COPY<LocalType, Field, Traits, ConstantMatTraitsType>::EXEC(d, a);
    }

    /** 
     * \brief Construct a Mat<...> from the given array of Field. Use the
     *        static function Traits::index to determine the element
     *        order. It is assumed that the element order of the given
     *        array is the same as that of this mat<...>.
     *
     * \note If C++ were to allow us to provide a template argument to
     *       a templated function then we could templatize this so
     *       that we could associate different traits with the input
     *       array.q
     *
     */
    template<typename In_Type>
    Mat(const In_Type* const data, size_t ldim=NROW)
    {
      COPY<LocalType, const In_Type* const, Traits, Traits>::EXEC(d, data);
    }

    /** 
     * \brief Copy construct a Mat<...> from the given Mat<...>. Use
     *        the static function Traits::index to determine the
     *        element order of this matrix and the other
     *        matrix. Generalize later.
     */
    Mat(const ThisType& m) 
    { 
      COPY<LocalType, ThisType, Traits, Traits>::EXEC(d, m);
    }

    ~Mat() {}

    /// Element in row i and column j (using zero indexing).
    const Field& operator () (size_t i,size_t j) const
    {
      ASSERT(i<NROW,
	     std::range_error("i>=NROW in Mat<Field,NROW,NCOL>::operator()(size_t i,size_t j)")); 
      ASSERT(j<NCOL,
	     std::range_error("j>=NCOL in Mat<Field,NROW,NCOL>::operator()(size_t i,size_t j)")); 
      return d[Traits::index(i,j)];
    }

    /// Element in row i and column j
    Field& operator () (size_t i,size_t j) 
    {
      ASSERT(i<NROW,
	     std::range_error("i>=NROW in Mat<Field,NROW,NCOL>::operator()(size_t i,size_t j)")); 
      ASSERT(j<NCOL,
	     std::range_error("j>=NCOL in Mat<Field,NROW,NCOL>::operator()(size_t i,size_t j)")); 
      return d[Traits::index(i,j)];
    }

    /// Access elements in Traits order
    Field& operator[] (size_t i)
    {
      //      ASSERT(i<NROW*NCOL,
      //       std::range_error("i>=NROW*NCOL in Mat<Field,NROW,NCOL>::operator[](size_t i)"));
      return d[i];
    }

    const Field& operator[] (size_t i) const
    {
      //ASSERT(i<NROW*NCOL,
      //       std::range_error("i>=NROW*NCOL in Mat<Field,NROW,NCOL>::operator[](size_t i)"));
      return d[i];
    }

    /// Matrix assignment
    ThisType& operator = (const ThisType& m)
    {
      COPY<LocalType, ThisType, Traits, Traits>::EXEC(d,m);
      //      for(size_t i=0;i<NROW*NCOL;i++) *((*d)+i) = *((*m.d)+i);
      return *this;
    }

    /// Assignment of scalar to each element
    ThisType& operator = (const Field& a)
    {
      COPY<LocalType, Field, Traits, ConstantMatTraitsType>::EXEC(d, a);
      //for(size_t i=0;i<NROW*NCOL;i++) *((*d)+i) = a;
      return *this;
    }

    /**
     *\brief Swap the specified rows.
     */
    void swapRows(size_t r1, size_t r2) {
      //      ASSERT(r1 < NROW,
      //	     std::range_error("Mat.swapRows() row 1 out of range!"));
      //ASSERT(r2 < NROW,
      //	     std::range_error("Mat.swapRows() row 2 out of range!"));
      
      for(size_t col=0; col<NCOL; col++) 
	std::swap(d[Traits::index(r1,col)], d[Traits::index(r2,col)]);
    }
    
    /**
     *\brief Swap the specified rows.
     */
    void swapCols(size_t c1, size_t c2) {
      //      ASSERT(r1 < NROW,
      //	     std::range_error("Mat.swapRows() row 1 out of range!"));
      //ASSERT(r2 < NROW,
      //	     std::range_error("Mat.swapRows() row 2 out of range!"));
      
      for(size_t row=0; row<NROW; row++) 
	std::swap(d[Traits::index(row,c1)], d[Traits::index(row,c2)]);
    }
    
    /**
     *\brief Swap the specified rows.
     */
    void negateRow(size_t r) {
      //      ASSERT(r1 < NROW,
      //	     std::range_error("Mat.swapRows() row 1 out of range!"));
      //ASSERT(r2 < NROW,
      //	     std::range_error("Mat.swapRows() row 2 out of range!"));
      static Field minusOne(-1);
      for(size_t col=0; col<NCOL; col++) 
	d[Traits::index(r,col)] *= minusOne;
    }
    
    /**
     *\brief Negate the specified Column.
     */
    void negateCol(size_t c) {
      //      ASSERT(r1 < NROW,
      //	     std::range_error("Mat.swapRows() row 1 out of range!"));
      //ASSERT(r2 < NROW,
      //	     std::range_error("Mat.swapRows() row 2 out of range!"));
      static Field minusOne(-1);
      for(size_t row=0; row<NROW; row++) 
	d[Traits::index(row,c)] *= minusOne;
    }
    
  };

  //======================================================================

  /**
   * \brief The Mat<...>  >> operator template.
   *
   * This is the standard parser for Mat objects. It reads a sequence
   * of element from the input stream. The element type, T, must
   * provide a >> operator. 
   *
   * If the input stream does not provide enough elements, then the
   * operator will throw an error. The operator does not care if the
   * input stream is capable of providing more than enough elements.
   *
   * \note The parsing is buffered, the input is read into a local mat
   *       object and only copied to the output object if the read is
   *       successful.
   *
   */
  template <typename T, size_t NROW, size_t NCOL>
  std::istream& operator >> (std::istream& is, Mat<T,NROW,NCOL>& m) {

    Mat<T,NROW,NCOL>& buffer;

    size_t numRead = 0;

    while(is >> buffer[numRead++] && numRead <= NROW*NCOL );
    
    if( numRead ==  NROW*NCOL ) {
      m = buffer;
      return;
    }
    std::ostringstream outMessage;
    outMessage << " operator>>(Mat<...>) read only " << numRead 
	       << " elements, needed " << NROW*NCOL;
    throw std::range_error(outMessage.str());
  }
  
  /**
   * \brief The Mat<...>  << operator template.
   *
   * This is the standard printer for Mat objects. It writes a line,
   * then a series of rows and then a line.
   *
   */
  template<typename T, size_t NROW,size_t NCOL, typename Traits>
  std::ostream& operator << (std::ostream& os, const Mat<T,NROW,NCOL,Traits> &m) {
      
    os << "Mat: " << " ------------------------------" << std::endl;
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.precision(6);
    os.width(15);
    for(size_t i=0; i< NROW; ++i) {
      os.width(15);
      for(size_t j=0; j< NCOL; ++j) 
	os << m(i,j) << "  ";
      os << std::endl;
    } 
    os << "End Mat ------------------------------" << std::endl;
    return os;
  }

} /* namespace psimag */



#endif /* PSIMAG_Mat_H */
