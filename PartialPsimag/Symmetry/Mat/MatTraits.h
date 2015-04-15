//-*-C++-*-

#ifndef PSIMAG_MatTraits_H
#define PSIMAG_MatTraits_H

/** \file MatTraits.h
 *
 *  \brief Contains template classes which store information about a
 *         matrix and to perform type caculations.
 *
 */
 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"

namespace psimag {

  /**
   * \brief Template class representing row major matrix ordering.
   *
   * The class does not have any non-static members and is not
   * intended for instantiation. It provides constant definitions and
   * static member functions that determine the order of the matrix's
   * storage. (e.g. row-major vs. column major).
   *
   */
  template<typename T, size_t NROW_, size_t NCOL_=NROW_> 
  class RowMajorTraits
  {
  public:
    
    typedef T            ElType;
    typedef T            Type[NROW_*NCOL_];
    typedef Type&        RefType;
    typedef const Type   ConstType;
    typedef ConstType&   ConstRefType;  
    
    enum { NumElements = NROW_*NCOL_, 
	   NROW        = NROW_, 
	   NCOL        = NCOL_  };

    static char* name() {return "RowMajor";}

    template<size_t ROW_ , size_t COL_>
    class Index {public: enum { VAL=ROW_*NROW_+COL_ }; };

    static size_t index(size_t row, size_t col) {
      return row*NROW_+col;
    }
    
    template<typename MatType, size_t ROW_, size_t COL_>
    class REF
    {
    public:
      static T&  GET(MatType& m) {
	return m[Index<ROW_,COL_>::VAL];
      }
      static const T&  GETCONST(const MatType& m) {
	return m[Index<ROW_,COL_>::VAL];
      }
    };
  };

  //======================================================================
  //======================================================================

  /**
   * \brief Template class representing column major matrix ordering.
   *
   * The class does not have any non-static members and is not
   * intended for instantiation. It provides constant definitions and
   * static member functions that determine the order of the matrix's
   * storage. (e.g. row-major vs. column major).
   *
   */
  template<typename T, size_t NROW_, size_t NCOL_=NROW_> 
  class ColMajorTraits
  {
  public:
    
    typedef T            ElType;
    typedef T            Type[NROW_*NCOL_];
    typedef Type&        RefType;
    typedef const Type   ConstType;
    typedef ConstType&   ConstRefType;  
    
    enum { NumElements = NROW_*NCOL_, 
	   NROW        = NROW_, 
	   NCOL        = NCOL_  };

    static const char* name() {return "ColMajor";}\

    template<size_t ROW_ , size_t COL_>
    class Index { public: enum { VAL=COL_*NROW_+ROW_ }; };

    static size_t index(size_t row, size_t col) {
      return col*NROW_+row;
    }
    
    template<typename MatType, size_t ROW_, size_t COL_>
    class REF
    {
    public:
      static T&  GET(MatType& m) {
	return m[Index<ROW_,COL_>::VAL];
      }
      static const T&  GETCONST(const MatType& m) {
	return m[Index<ROW_,COL_>::VAL];
      }
    };
    
  };
  
  //======================================================================
  //======================================================================

  /**
   * \brief Template class representing A Constant Matrix
   *
   * The MatType is the same as the Element type.
   *
   */
  template<typename T, size_t NROW_, size_t NCOL_=NROW_> 
  class ConstantTraits
  {
  public:
    
    typedef T            ElType;
    typedef T            Type[NROW_*NCOL_];
    typedef Type&        RefType;
    typedef const Type   ConstType;
    typedef ConstType&   ConstRefType;  
    
    enum { NumElements = NROW_*NCOL_, 
	   NROW        = NROW_, 
	   NCOL        = NCOL_  };

    static char* name() {return "ConstantTraits";}\

    template<size_t ROW_ , size_t COL_>
    class Index { 
    public: 
      enum { VAL=0 }; 
    };

    static size_t index(size_t row, size_t col) {
      return 0;
    }
    
    template<typename MatType, size_t ROW_, size_t COL_>
    class REF
    {

    public:
      static T&  GET(MatType& m) {
	return m;
      }
      static const T&  GETCONST(const MatType& m) {
	return m;
      }
    };
    
  };
  
  //======================================================================

  /**
   * \brief Template class representing double index ordering.
   *
   * The class does not have any non-static members and is not
   * intended for instantiation. It provides constant definitions and
   * static member functions that determine the order of the matrix's
   * storage. (e.g. row-major vs. column major).
   *
   */
  template<typename T, size_t NROW_, size_t NCOL_=NROW_> 
  class DoubleIndexTraits
  {
  public:
    
    typedef T            ElType;
    
    enum { NumElements = NROW_*NCOL_, 
	   NROW        = NROW_, 
	   NCOL        = NCOL_  };

    static char* name() {return "DoubleIndex";}

    template<typename MatType, size_t ROW_, size_t COL_>
    class REF
    {
    public:
      static T&  GET(MatType& m) {
	return m[ROW_][COL_];
      }
      static const T&  GETCONST(const MatType& m) {
	return m[ROW_][COL_];
      }
    };
    
  };
  
  //======================================================================

  /**
   * \brief Template class representing double index ordering.
   *
   * The class does not have any non-static members and is not
   * intended for instantiation. It provides constant definitions and
   * static member functions that determine the order of the matrix's
   * storage. (e.g. row-major vs. column major).
   *
   */
  template<typename T, size_t NROW_, size_t NCOL_=NROW_> 
  class ReverseDoubleIndexTraits
  {
  public:
    
    typedef T            ElType;
    
    enum { NumElements = NROW_*NCOL_, 
	   NROW        = NROW_, 
	   NCOL        = NCOL_  };

    static char* name() {return "DoubleIndex";}\

    template<typename MatType, size_t ROW_, size_t COL_>
    class REF
    {
    public:
      static T&  GET(MatType& m) {
	return m[COL_][ROW_];
      }
      static const T&  GETCONST(const MatType& m) {
	return m[COL_][ROW_];
      }
    };
    
  };
  
} /* namespace psimag */

#endif /* PSIMAG_MatTraits_H */
