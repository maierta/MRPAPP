//-*-C++-*-

#ifndef PSIMAG_ForEachCentering_H
#define PSIMAG_ForEachCentering_H


/** \file ForEachCentering.h
 *
 *  \brief Contains 'meta programming' template functions providing a
 *         ForEach operation over the Centerings. This operation
 *         applies a given function to centring for the specified
 *         dimension.
 *
 */
 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "PSIMAGAssert.h"
#include "MatTraits.h"

namespace psimag {


  //======================================================================
  /**
   * \brief The ForEachCentering_ template that calls
   *        Functor::ForCentering<Centering<DIM,NUM> For NUM, NUM-1 .. 0.
   */
  template<typename Field, size_t DIM, size_t NUM, 
	   template<typename, typename> class Functor, 
	   typename ArgType>
  class ForEachCentering_ {
  public:

    typedef Centering<Field,DIM,NUM> CenteringType;
    
    static bool EXEC(ArgType& arg) {
      bool success = Functor<CenteringType, ArgType>::ForCentering(arg);
      if (success) return success;
      return ForEachCentering_<Field, DIM, NUM-1, Functor, ArgType>::EXEC(arg);
    }
  };

  //======================================================================
  /**
   * \brief Terminating specialization.
   *
   */
  template<typename Field, size_t DIM,
	   template<typename, typename> class Functor, 
	   typename ArgType>
  class ForEachCentering_<Field,DIM,0,Functor,ArgType> {
  public:
    enum {NUM=0};
    typedef Centering<Field,DIM,NUM> CenteringType;

    static bool EXEC(ArgType& arg) {
      bool success = Functor<CenteringType, ArgType>::ForCentering(arg);
      return success;
    }
  };

  //======================================================================
  /**
   * \brief The FOREACH template that calls FOREACH_ to start the
   *        recursive computation of the aggregate operation at
   *        ROWS_LEFT =  NROW.
   *
   */
  template<typename Field, size_t DIM, 
	   template<typename, typename> class Functor,
	   typename ArgType>
  class ForEachCentering {
  public:
    
    enum { NUM=Centerings<2>::NumCenterings };

    static bool EXEC(ArgType& arg) {
      bool success = ForEachCentering_<Field,DIM,NUM,Functor,ArgType>::EXEC(arg);
      return success;
    }

  };

} /* namespace psimag */

#endif /* PSIMAG_Mat_ForEach_H */
