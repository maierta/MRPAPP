//-*-C++-*-

#ifndef PSIMAG_SymmetryOperations2D_H
#define PSIMAG_SymmetryOperations2D_H

/** \ingroup symmetryConcepts */
/*@{*/

/** \file  SymmetryOperations2D.h
 *
 *  Contains 2D symmetry operations that are used to classify 2D space groups.
 *  
 */
 
#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Mat.h"
#include "Real.h"

#include "SymmetryOperation.h"
#include "CellDirection.h"
#include "CellPosition.h"
#include "CellTranslation.h"
#include "CartesianTranslation.h"
#include "FieldConvert.h"

// #include "Mirror2D.h"
// #include "Glide2D.h"
// #include "TwoFold2D.h"
// #include "ThreeFold2D.h"
// #include "FourFold2D.h"
// #include "SixFold2D.h"
// #include "ThreeFoldN2D.h"
// #include "FourFoldN2D.h"
// #include "SixFoldN2D.h"


#include "PSIMAGAssert.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Algorithms>
  class SymmetryOperations {};

  template<typename Field, typename Algorithms>
  class SymmetryOperations<Field,2,Algorithms> {
  public:

    enum {DIM=2, NUMOPS=51};

    typedef SymmetryOperation<Field,2,Algorithms>                          SymmetryOperationType;
    typedef typename SymmetryOperationType::TranslationType                TranslationType;
    typedef Lattice<Field,DIM,Algorithms>                                  LatticeType;
    typedef CartesianTranslation<Field,2>                                  CartesianTranslationType;
    typedef LatticeCoordinates<DIM>                                        LatticeCoordinatesType;

    typedef std::vector<SymmetryOperationType>                             RatedOperationsType;

//     typedef Mirror<Field,2,Algorithms>                                     Mirror;
//     typedef Glide<Field,2,Algorithms>                                      Glide;

//     typedef TwoFold<Field,2,Algorithms>                                    TwoFold;
//     typedef ThreeFold<Field,2,Algorithms>                                  ThreeFold;
//     typedef FourFold<Field,2,Algorithms>                                   FourFold;
//     typedef SixFold<Field,2,Algorithms>                                    SixFold;

//     typedef ThreeFoldN<Field,2,Algorithms>                                 ThreeFoldN;
//     typedef FourFoldN<Field,2,Algorithms>                                  FourFoldN;
//     typedef SixFoldN<Field,2,Algorithms>                                   SixFoldN;

    //====================================================================== identity

    /**
     * \brief Returns the cartrsian SymmetryOperation for a mirror about the
     *        given axis.
     */
    static SymmetryOperationType ident() {
      SymmetryOperationType result;
      result.name         = "identity";
      result.element.name = "identity";
      return result;
    }						   
    
 //    /**
//      * \brief Returns the cartrsian centering translation SymmetryOperation.
//      */
//     static SymmetryOperationType centering() {
//       SymmetryOperationType result(cartesianTranslation<Field>("1/2","1/2"));
//       result.name         = "Centering";
//       result.element.name = "centering";
//       return result;
//     }						   

    static void load(const LatticeType& lattice, 
		     RatedOperationsType& result) {
      
      result.reserve(NUMOPS);

      result[0]  = ident();

      //      result[1]  = TwoFold() (lattice);
      //      result[2]  = TwoFold("0",   "1/2") (lattice);
      //      result[3]  = TwoFold("1/2", "0" ) (lattice);
      //      result[4]  = TwoFold("1/2", "1/2") (lattice);

//       result[5]  = ThreeFold("1/3","1/3") (lattice);
//       result[6]  = ThreeFold("2/3","2/3") (lattice);

//       result[7]  = ThreeFoldN("1/3","1/3") (lattice);
//       result[8]  = ThreeFoldN("2/3","2/3") (lattice);

//       result[9]  = FourFold() (lattice);
//       result[10] = FourFold("1/2","1/2") (lattice);

//       result[11] = FourFoldN() (lattice);
//       result[12] = FourFoldN("1/2","1/2") (lattice);

//       result[13] = SixFold() (lattice);
//       result[14] = SixFoldN() (lattice);
	
//       result[15] = Mirror::a() (lattice);
//       result[16] = Mirror::a("1/4") (lattice);
//       result[17] = Mirror::a("1/2") (lattice);
//       result[18] = Mirror::a("3/4") (lattice);
     
//       result[19] = Mirror::b() (lattice);
//       result[20] = Mirror::b("1/4") (lattice);
//       result[21] = Mirror::b("1/2") (lattice);
//       result[22] = Mirror::b("3/4") (lattice);

//       // A + B mirror     
//       result[23] = Mirror(latticeCoordinate(1,1), "a+b Mirror") (lattice);                              
//       result[24] = Mirror(latticeCoordinate(1,1), latticeCoordinate(1,0), "1/2", "a+b Mirror 1/2a->") (lattice); // a offset
//       result[25] = Mirror(latticeCoordinate(1,1), latticeCoordinate(0,1), "1/2", "a+b Mirror 1/2b->") (lattice); // b offset

//       // A - B mirror     
//       result[26] = Mirror(latticeCoordinate(1,-1), latticeCoordinate(0,1), "1/2", "a-b Mirror 1/2b->") (lattice); 
//       result[27] = Mirror(latticeCoordinate(1,-1), latticeCoordinate(0,1), "1",   "a-b Mirror b->") (lattice);                              
//       result[28] = Mirror(latticeCoordinate(1,-1), latticeCoordinate(0,1), "3/2", "a-b Mirror 3/2b->") (lattice); 

//       // -A + 2B mirror (Hexagonal)
//       result[29] = Mirror(latticeCoordinate(-1,2), latticeCoordinate(1,0), "1/2", "-a+2b Mirror 1/2a->") (lattice);
//       result[30] = Mirror(latticeCoordinate(-1,2), latticeCoordinate(1,0), "1",   "-a+2b Mirror a->") (lattice);   // Looks nicer here rather than at 0

//       // 2A - B mirror
//       result[31] = Mirror(latticeCoordinate(2,-1), latticeCoordinate(0,1), "1/2", "2a-b Mirror 1/2b->") (lattice);
//       result[32] = Mirror(latticeCoordinate(2,-1), latticeCoordinate(0,1), "1",   "2a-b Mirror b->") (lattice);   // Looks nicer here rather than at 0

//       result[33] = Glide::a()      (lattice);
//       result[34] = Glide::a("1/4") (lattice);
//       result[35] = Glide::a("1/2") (lattice);
//       result[36] = Glide::a("3/4") (lattice);
     
//       result[37] = Glide::b()      (lattice);
//       result[38] = Glide::b("1/4") (lattice);
//       result[39] = Glide::b("1/2") (lattice);
//       result[40] = Glide::b("3/4") (lattice);

//       // A + B glide     
//       result[41] = Glide(latticeCoordinate(1,1), convert<Field>("1/2"),
// 			 "a+b Glide") (lattice);                              
//       result[42] = Glide(latticeCoordinate(1,1), convert<Field>("1/2"),
// 			 latticeCoordinate(1,0), convert<Field>("1/2"), 
// 			 "a+b Glide 1/2a->") (lattice); // a offset
//       result[43] = Glide(latticeCoordinate(1,1), convert<Field>("1/2"), 
// 			 latticeCoordinate(0,1), convert<Field>("1/2"),
// 			 "a+b Glide 1/2b->") (lattice); // b offset

//       // A - B glide     
//       result[44] = Glide(latticeCoordinate(1,-1),  convert<Field>("1/2"),
// 			 latticeCoordinate(0,1),  convert<Field>("1/2"),
// 			 "a-b Glide 1/2b->") (lattice); 
//       result[45] = Glide(latticeCoordinate(1,-1),  convert<Field>("1/2"),
// 			 latticeCoordinate(0,1),  convert<Field>("1"),
// 			 "a-b Glide    b->") (lattice);                              
//       result[46] = Glide(latticeCoordinate(1,-1),  convert<Field>("1/2"),
// 			 latticeCoordinate(0,1),  convert<Field>("3/2"),
// 			 "a-b Glide 3/2b->") (lattice); 

//       // -A + 2B Glide (Hexagonal)
//       result[47] = Glide(latticeCoordinate(-1,2),  convert<Field>("1/2"),
// 			 latticeCoordinate(1,0),  convert<Field>("1/4"),
// 			 "-a+2b Glidew 1/4a->") (lattice);
//       result[48] = Glide(latticeCoordinate(-1,2),  convert<Field>("1/2"),
// 			 latticeCoordinate(1,0),  convert<Field>("3/4"),
// 			 "-a+2b Glide 3/4a->") (lattice);   // Looks nicer here rather than at 0

//       // 2A - B Glide
//       result[49] = Glide(latticeCoordinate(2,-1),  convert<Field>("1/2"),
// 			 latticeCoordinate(0,1),  convert<Field>("1/4"),
// 			 "2a-b Glide 1/4b->") (lattice);
//       result[50] = Glide(latticeCoordinate(2,-1),  convert<Field>("1/2"),
// 			 latticeCoordinate(0,1),  convert<Field>("3/4"),
// 			 "2a-b Glide 3/4b->") (lattice);   // Looks nicer here rather than at 0

      for( size_t i=0; i<NUMOPS; i++) {
	result[i].id = i;
      }
      
    }

  };


} /* namespace psimag */
#endif 
/*@}*/
