//-*-C++-*-

#ifndef PSIMAG_MatPrint_H
#define PSIMAG_MatPrint_H

/** \file MatPrint.h
 *
 *  \brief Contains template functions providing print operations on Matrices of any sort.
 *
 */
 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"
#include "MatForEach.h"

namespace psimag {

  using std::setw;

  /**
   * \brief The primary template class for the Print operation.
   *
   * \param MatType:  The type of the matrix used on the LHS.
   * \param MatTypeO: The type of the matrix used on the RHS.
   * \param Traits:   A type (e.g. MatTraits<NROW,NCOL>) describng the matrix.
   * \param TraitO:   A type (e.g. MatTraits<NROW,NCOL>) describng the other matrix.
   *
   */
  template<typename MatTypeLHS, 
	   typename TraitsLHS>
  class PRINT_FUNCTOR {
  public:

    enum { NROW=TraitsLHS::NROW, NCOL=TraitsLHS::NCOL };
    
    class PrintFunctionType {
    public:
      
      typedef typename TraitsLHS::ElType LHSType;
      
      static void EXEC(const LHSType& lhs, size_t row, size_t col, std::ostream& os) {
	if (row != 0 && col == 0)
	  os << std::endl << " ";
	os << std::setw(10) << lhs << " ";
      }
    };

    class FlatPrintFunctionType {
    public:
      
      typedef typename TraitsLHS::ElType LHSType;
      
      static void EXEC(const LHSType& lhs, size_t row, size_t col, std::ostream& os) {
	os << std::setw(10) << lhs << " ";
      }
    };

    static void EXEC(const MatTypeLHS& lhs, std::ostream& os, bool flat=false) { 
      if (flat)
	return FOREACH<FlatPrintFunctionType, MatTypeLHS, TraitsLHS, std::ostream>::EXEC(lhs, os);
      else
	return FOREACH<PrintFunctionType, MatTypeLHS, TraitsLHS, std::ostream>::EXEC(lhs, os);
    }
  };

  //======================================================================

  template<typename MatType, typename Traits > 
  class MAT_PRINT {
  public:

    typedef typename Traits::ElType T;

    static void EXEC(const MatType& m, std::ostream& os, bool flat=false) {

      os.setf(std::ios_base::fixed, std::ios_base::floatfield);
      os.precision(6);

      os << "(" << Traits::name() << ")";
      if(!flat) os << std::endl;
      os << "[";
      PRINT_FUNCTOR<MatType, Traits>::EXEC(m, os, flat);
      os << "]";
    }

    static void JUST_NUMBERS(const MatType& m, std::ostream& os, bool flat=false) {

      os.setf(std::ios_base::fixed, std::ios_base::floatfield);
      os.precision(6);

      if(!flat) os << std::endl;
      PRINT_FUNCTOR<MatType, Traits>::EXEC(m, os, flat);
    }
  };

  //======================================================================

  template<typename T, 
	   template<typename, size_t, size_t> class MatType, 
	   size_t NROW, size_t NCOL>
  void Mat_Print(const MatType<T,NROW,NCOL>& m, std::ostream& os, bool flat=true) {
    return MAT_PRINT<MatType<T,NROW, NCOL>, typename MatType<T,NROW, NCOL>::Traits>::EXEC(m, os, flat); 
  }

#define Mat_Print_Template(T, DIM)					\
  void  Mat_Print(ColMajorTraits<T,DIM,DIM>::ConstRefType m, std::ostream& os, bool flat=true) {		\
    return MAT_PRINT< ColMajorTraits<T,DIM,DIM>::Type, ColMajorTraits<T,DIM,DIM> >::EXEC(m, os, flat); }

/**  \brief Overload to print a  2x2 double array. */
inline Mat_Print_Template(double,2)
  
/**  \brief Overload to print a 3x3 determinant of a double array. */
inline Mat_Print_Template(double,3)

/**  \brief Overload to print a 2x2 determinant of an int array. */
inline Mat_Print_Template(int,2)

/**  \brief Overload to print a 3x3 determinant of an int array. */
inline Mat_Print_Template(int,3)

  //====================================================================== << operator

  template<typename T, 
	   template<typename, size_t, size_t> class MatType, 
	   size_t NROW, size_t NCOL>
inline  std::ostream& operator << (std::ostream& os, const MatType<T,NROW,NCOL>& m) {
    MAT_PRINT<MatType<T,NROW, NCOL>, typename MatType<T,NROW, NCOL>::Traits>::EXEC(m, os); 
    return os;
  }

#define Mat_Print_OP_Template(T, DIM)					\
  std::ostream& operator << (std::ostream& os, ColMajorTraits<T,DIM,DIM>::ConstRefType m) { \
    MAT_PRINT< ColMajorTraits<T,DIM,DIM>::Type, ColMajorTraits<T,DIM,DIM> >::EXEC(m, os); \
    return os;								\
  }

/**  \brief Overload to print a  2x2 double array. */
inline Mat_Print_OP_Template(double,2)
  
/**  \brief Overload to print a 3x3 determinant of a double array. */
inline Mat_Print_OP_Template(double,3)

/**  \brief Overload to print a 2x2 determinant of an int array. */
inline Mat_Print_OP_Template(int,2)

/**  \brief Overload to print a 3x3 determinant of an int array. */
inline Mat_Print_OP_Template(int,3)




} /* namespace psimag */

#endif /* PSIMAG_MatPrint_H */
