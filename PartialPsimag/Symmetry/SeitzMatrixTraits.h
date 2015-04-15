//-*-C++-*-

#ifndef PSIMAG_SeitzMatrixTraits_H
#define PSIMAG_SeitzMatrixTraits_H

/** \file SeitzMatrixTraits.h
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
#include "TypeManip.h"


namespace psimag {

  /**
   * \brief Template class representing how SeitzMatricies are
   *        constructed. This class is used by matrix algorithms to
   *        construct code that accesses the proper data.
   *
   * The class does not have any non-static members and is not
   * intended for instantiation. It provides constant definitions and
   * static member functions that determine the order of the matrix's
   * storage. (e.g. row-major vs. column major).
   *
   */
  template<typename T, size_t DIM, 
	   size_t NROW_=DIM+1, 
	   size_t NCOL_=DIM+1> 
  class SeitzMatrixTraits
  {
  public:

    enum { NumElements = NROW_*NCOL_, 
	   NROW        = NROW_, 
	   NCOL        = NCOL_ };

    typedef T ElType;

    static const char* name() {return "SeitzMatrix";}
    
    template<typename SeitzMatrixType, size_t ROW_, size_t COL_>
    class REF
    {
      typedef typename SeitzMatrixType::RotationType      RotationType_;
      typedef typename SeitzMatrixType::RotationTraits    RotationTraits_;
      typedef typename RotationTraits_::template REF<RotationType_, ROW_, COL_> RREF;
      
      class RotRetriever {
      public:
	static T& GET(SeitzMatrixType& m) {
	  return RREF::GET(m.rotation);
	}
	static const T& GETCONST(const SeitzMatrixType& m) {
	  return RREF::GETCONST(m.rotation);
	}
      };
      class TransRetriever {
      public:
	static T& GET(SeitzMatrixType& m) {
	  return m.translation[ROW_];
	}
	static const T& GETCONST(const SeitzMatrixType& m) {
	  return m.translation[ROW_];
	}
      };
      class ExtRetriever {
      public:
	static T& GET(SeitzMatrixType& m) {
	  return m.lastRow[COL_];
	}
	static const T& GETCONST(const SeitzMatrixType& m) {
	  return m.lastRow[COL_];
	}
      };
      
      static const bool BOTTOM_ROW = (ROW_ == DIM);
      static const bool LAST_COL   = (COL_ == DIM);
      
      typedef typename Select<LAST_COL,
			      TransRetriever,
			      RotRetriever >::Result  LastColRetriever;

      typedef typename Select<BOTTOM_ROW, 
			      ExtRetriever,
			      LastColRetriever> ::Result RETRIEVER;
    public:

      static T&  GET(SeitzMatrixType& m) {
	return RETRIEVER::GET(m);
      }

      static const T&  GETCONST(const SeitzMatrixType& m) {
	return RETRIEVER::GETCONST(m);
      }
    };

  };
  
} /* namespace psimag */

#endif /* PSIMAG_SeitzMatrixTraits_H */
