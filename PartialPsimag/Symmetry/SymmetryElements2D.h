//-*-C++-*-

#ifndef PSIMAG_SymmetryElements2D_H
#define PSIMAG_SymmetryElements2D_H

/** \ingroup symmetryConcepts */
/*@{*/

/** \file  SymmetryElements2D.h
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
#include "SymmetryElement.h"
#include "AppliedSymmetryElement.h"
#include "CellDirection.h"
#include "CellPosition.h"
#include "CellTranslation.h"
#include "CartesianTranslation.h"
#include "FieldConvert.h"
#include "LatticeWithPattern.h"
#include "ReducedLattice.h"

#include "IdentityElement2D.h"
#include "Mirror2D.h"
#include "Glide2D.h"
#include "TwoFold2D.h"
#include "ThreeFold2D.h"
#include "FourFold2D.h"
#include "SixFold2D.h"
#include "ThreeFoldN2D.h"
#include "FourFoldN2D.h"
#include "SixFoldN2D.h"


#include "PSIMAGAssert.h"

namespace psimag {
  
  template<typename Field, size_t DIM, typename Occupant,typename LatticeType,typename Algorithms>
  class AppliedSymmetryElement;

  template<typename Field, size_t DIM, typename Algorithms>
  class SymmetryElements {};

  template<typename Field, typename Algorithms>
  class SymmetryElements<Field,2,Algorithms> {
  public:

    enum {DIM=2};

    typedef SymmetryElement<Field,DIM,Algorithms>                          SymmetryElementType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                        SymmetryOperationType;
    typedef typename SymmetryOperationType::TranslationType                TranslationType;

    typedef CartesianTranslation<Field,DIM>                                CartesianTranslationType;
    typedef LatticeCoordinates<DIM>                                        LatticeCoordinatesType;

    typedef std::vector<SymmetryOperationType>                             RatedElementsType;

    typedef Mirror<Field,DIM,Algorithms>                                   MirrorType;
    typedef Glide<Field,DIM,Algorithms>                                    GlideType;

    typedef IdentityElement<Field,DIM,Algorithms>                          IdentityElementType;
    typedef TwoFold<Field,DIM,Algorithms>                                  TwoFoldType;
    typedef ThreeFold<Field,DIM,Algorithms>                                ThreeFoldType;
    typedef FourFold<Field,DIM,Algorithms>                                 FourFoldType;
    typedef SixFold<Field,DIM,Algorithms>                                  SixFoldType;

    typedef ThreeFoldN<Field,DIM,Algorithms>                               ThreeFoldNType;
    typedef FourFoldN<Field,DIM,Algorithms>                                FourFoldNType;
    typedef SixFoldN<Field,DIM,Algorithms>                                 SixFoldNType;

    /**
     *\brief This is the 2D version of this function.
     *
     * It's operation is different in 3D where an exhaustive search may not make sense.
     */
    template<typename Occupant, typename LatticeType>
    static
    std::vector< AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> >
    validElements(const  std::vector< AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> >& appliedEls) {
	
      typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
      typedef std::vector< AppliedSymmetryElementType >                         AppliedSymmetryElementsType;
      
      AppliedSymmetryElementsType result;
      static Field threshold = Algorithms::threshold();

      for (size_t i=0; i< appliedEls.size(); i++) {
	const AppliedSymmetryElementType& appliedEl = appliedEls[i];
	Field d = appliedEl.distance;

	if (d < threshold)
	  result.push_back(appliedEl);
      }
      
      return result;
    }
      
    //======================================================================
    //
    /**
    * \brief Return a vector of AppliedSymmetryElements which
    *        correspond to the SymmetryElements of this class.
    * 
    * \note While it might seem logical to put this method in the
    *       LatticeWithPattern class, it is a 2D method and we
    *       probably don't want the extra complexity of making 2d and
    *       3d versions of LatticeWithPattern.
    *
    */
    template<typename Occupant, typename LatticeType>
    static 
    std::vector< AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> > 
    allAppliedElements(const LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>& latpat) {
      
      typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
      typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>     LatticeWithPatternType;
      
      static const std::vector<SymmetryElementType*>  symmetryElements = elements_();
      std::vector<AppliedSymmetryElementType>  result;
      result.reserve(symmetryElements.size());

      for( size_t i=0; i<symmetryElements.size(); i++) {
	
	//const Lattice<Field,DIM,Algorithms>& lat = latpat;
	const SymmetryElementType& symEl = (*symmetryElements[i]);
	//SymmetryOperation<Field,DIM,Algorithms> so = symEl(lat);
	AppliedSymmetryElementType appliedElement(i, latpat, symEl);
	result.push_back(appliedElement);

      }
      
      return result;
    }

    /**
     *\brief This is the 2D version of this function.
     *
     * It's operation is different in 3D where an exhaustive search may not make sense.
     */
    template<typename Occupant, typename LatticeType>
    static
    std::vector< AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> >
    validElements(const LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>&             latpat,
		  std::vector< AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> > allAppliedElements_) {
      return validElements(allAppliedElements_);
    }

    static int getIndexFor(std::string name) {
       
      static const std::vector<SymmetryElementType*>  symmetryElements = elements_();
      
      for (size_t i=0; i<symmetryElements.size(); i++)  {

	SymmetryElementType* symEl = symmetryElements[i];
      
	if (name == symEl->name)
	  return i;
      }      
      std::ostringstream errorMsg;
      errorMsg << "SymmetryElements2D.getIndexFor could not find an index for '" << name << "'";
      printNames(errorMsg);
      throw std::range_error(errorMsg.str());
    }

    static void printNames(std::ostream& out) {

      static const std::vector<SymmetryElementType*>  symmetryElements = elements_();
      
      for (size_t i=0; i<symmetryElements.size(); i++)  {

	SymmetryElementType* symEl = symmetryElements[i];
	out << "----------------------------------------  " << symEl->name << std::endl;      
	out << toXML(*symEl) << std::endl;      
      }
    }      
      
      
  private:

    static void delElements() {

      std::vector<SymmetryElementType*> elements = elements_();

      for(size_t i=0; i< elements.size(); i++) {
	delete elements[i];
      }
    }

    static std::vector<SymmetryElementType*> elements_() {
      
      static std::vector<SymmetryElementType*> result;
      
      if (result.size() > 0) return result;
      
      result.push_back( new IdentityElementType() );

      result.push_back( new TwoFoldType());
      result.push_back( new TwoFoldType("0",   "1/2"));
      result.push_back( new TwoFoldType("1/2", "0" ));
      result.push_back( new TwoFoldType("1/2", "1/2"));

      result.push_back( new ThreeFoldType());
      result.push_back( new ThreeFoldType("1/3","1/3"));
      result.push_back( new ThreeFoldType("2/3","2/3"));

      result.push_back( new ThreeFoldNType());
      result.push_back( new ThreeFoldNType("1/3","1/3"));
      result.push_back( new ThreeFoldNType("2/3","2/3"));

      result.push_back( new FourFoldType());
      result.push_back( new FourFoldType("1/2","1/2"));

      result.push_back( new FourFoldNType());
      result.push_back( new FourFoldNType("1/2","1/2"));

      result.push_back( new SixFoldType());
      result.push_back( new SixFoldNType());
	
      result.push_back(MirrorType::a());
      result.push_back(MirrorType::a("1/4"));
      result.push_back(MirrorType::a("1/2"));
      result.push_back(MirrorType::a("3/4"));
     
      result.push_back(MirrorType::b());
      result.push_back(MirrorType::b("1/4"));
      result.push_back(MirrorType::b("1/2"));
      result.push_back(MirrorType::b("3/4"));

      // A + B mirror     
      result.push_back( new MirrorType(latticeCoordinate(1,1)));                              
      result.push_back( new MirrorType(latticeCoordinate(1,1), latticeCoordinate(1,0), "1/2")); // a offset

      // A - B mirror     
      result.push_back( new MirrorType(latticeCoordinate(1,-1), latticeCoordinate(1,0), "1/2")); 
      result.push_back( new MirrorType(latticeCoordinate(1,-1)));                              
      //result.push_back( new MirrorType(latticeCoordinate(1,-1), latticeCoordinate(0,1), "1"  ));                              
      //result.push_back( new MirrorType(latticeCoordinate(1,-1), latticeCoordinate(0,1), "3/2")); 

      // -A + 2B mirror (Hexagonal)
      result.push_back( new MirrorType(latticeCoordinate(-1,2)));
      result.push_back( new MirrorType(latticeCoordinate(-1,2), latticeCoordinate(1,0), "1/2"));
      //result.push_back( new MirrorType(latticeCoordinate(-1,2), latticeCoordinate(1,0), "1"  ));   // Looks nicer here rather than at 0
      //result.push_back( new MirrorType(latticeCoordinate(-1,2)));   // Looks nicer here rather than at 0

      // 2A - B mirror
      result.push_back( new MirrorType(latticeCoordinate(2,-1)));
      //result.push_back( new MirrorType(latticeCoordinate(2,-1), latticeCoordinate(0,1), "1/2"));
      //result.push_back( new MirrorType(latticeCoordinate(2,-1), latticeCoordinate(0,1), "1"  ));   // Looks nicer here rather than at 0

      result.push_back(GlideType::a());
      result.push_back(GlideType::a("1/4"));
      result.push_back(GlideType::a("1/2"));
      result.push_back(GlideType::a("3/4"));
     
      result.push_back(GlideType::b());
      result.push_back(GlideType::b("1/4"));
      result.push_back(GlideType::b("1/2"));
      result.push_back(GlideType::b("3/4"));

      // A + B glide     
      result.push_back( new GlideType(latticeCoordinate(1,1), convert<Field>("1/2"))); //    "a+b Glide"));                              
      result.push_back( new GlideType(latticeCoordinate(1,1), convert<Field>("1/2"),
				 latticeCoordinate(1,0), convert<Field>("1/2"))); //    "a+b Glide 1/2a->")); // a offset
//       result.push_back( new GlideType(latticeCoordinate(1,1), convert<Field>("1/2"), 
// 				 latticeCoordinate(0,1), convert<Field>("1/2"))); //    "a+b Glide 1/2b->")); // b offset
      
      // A - B glide     
      result.push_back( new GlideType(latticeCoordinate(1,-1), convert<Field>("1/2"),
				 latticeCoordinate(1,0),  convert<Field>("1/2"))); //   "a-b Glide 1/2a->")); 
      result.push_back( new GlideType(latticeCoordinate(1,-1), convert<Field>("1/2")));   //   "a-b Glide
//       result.push_back( new GlideType(latticeCoordinate(1,-1), convert<Field>("1/2"),
// 				 latticeCoordinate(0,1),  convert<Field>("3/2"))); //   "a-b Glide 3/2b->")); 
      
      // -A + 2B Glide (Hexagonal)
      result.push_back( new GlideType(latticeCoordinate(-1,2), convert<Field>("1/2"))); //   "-a+2b Glide
      result.push_back( new GlideType(latticeCoordinate(-1,2), convert<Field>("1/2"),
				      latticeCoordinate(1,0),  convert<Field>("1/4"))); //   "-a+2b Glidew 1/4a->"));
//       result.push_back( new GlideType(latticeCoordinate(-1,2), convert<Field>("1/2"),
// 				      latticeCoordinate(1,0),  convert<Field>("3/4"))); //   "-a+2b Glide 3/4a->"));   // Looks nicer here rather than at 0
//       result.push_back( new GlideType(latticeCoordinate(-1,2), convert<Field>("1/2"),
// 				      latticeCoordinate(1,0),  convert<Field>("5/4"))); //   "-a+2b Glide 3/4a->"));   // Looks nicer here rather than at 0
      
      // 2A - B Glide
      result.push_back( new GlideType(latticeCoordinate(2,-1), convert<Field>("1/2"))); //   "2a-b Glide 
      result.push_back( new GlideType(latticeCoordinate(2,-1), convert<Field>("1/2"),
				      latticeCoordinate(1,0),  convert<Field>("1/2"))); //   "2a-b Glide 1/2a->"));
//       result.push_back( new GlideType(latticeCoordinate(2,-1), convert<Field>("1/2"),
// 				 latticeCoordinate(0,1),  convert<Field>("3/4"))); //  "2a-b Glide 3/4b->"));       // Looks nicer here rather than at 0
//       result.push_back( new GlideType(latticeCoordinate(2,-1), convert<Field>("1/2"),
// 				 latticeCoordinate(0,1),  convert<Field>("5/4"))); //  "2a-b Glide 5/4b->"));       // Looks nicer here rather than at 0

      for (size_t i=0; i<result.size(); i++) 
	result[i]->id = i;
      
      return result;
    }

  };


} /* namespace psimag */
#endif 
/*@}*/
