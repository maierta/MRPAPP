//-*-C++-*-

#ifndef PSIMAG_SearchTable_H
#define PSIMAG_SearchTable_H


/** \file SearchTable.h
 *
 *  \brief Contains generic template class for holding search state transition tables.
 *
 */
 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"


namespace psimag {

  //======================================================================
  /**
   *  \brief A generic template class for holding search state
   *         transition tables.
   *
   * \note This is the generic template it is not used directly only
   *       its specializations.
   */
  template<size_t DIM>
  class SearchTable {};

} /* namespace psimag */

#endif /* PSIMAG_SearchTable */
