//-*-C++-*-

#ifndef  Psimag_SpaceGroup2D
#define  Psimag_SpaceGroup2D

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SpaceGroup2D.h
 *
 *  1999, Grosse-Kunstleve, "Algorithms for deriving crystallographic space-group information"
 *
 *  some of this will be implemented in the cctbx where I may be able to reuse it.
 */  

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <bitset>
#include <iterator>

#include "SpaceGroupData2D.h"
#include "SymmetryOperation.h"
#include "SpaceGroup.h"
#include "CellPosition.h"
#include "LatticeWithPattern.h"
#include "Symmetry.h"
#include "Tag.h"
#include "SymmetryElements2D.h"
#include "AppliedSymmetryElementIterator.h"

namespace psimag {
  
  //====================================================================== 
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroup class specialization.
   */
  template<typename Field, typename Occupant, typename Algorithms>
  class SpaceGroup<Field,2,Occupant,Algorithms>   
  {
  public:

    enum {DIM=2, NUMGROUPS=17};

    typedef SpaceGroup<Field,2,Occupant, Algorithms>                          SpaceGroupType;
    typedef SymmetryElements<Field,DIM,Algorithms>                            SymmetryElementsType;  
    typedef CellPosition<Field,DIM,Algorithms>                                CellPositionType; 

    //    typedef std::bitset<NUMOPS> BitsType;

    bool                                    loaded;
    size_t                                  number;
    std::string                             name;
    bool                                    centered;
    Field                                   rating; 
    bool                                    selected;

    std::vector<int>                        operationIds;
    std::vector<CellPositionType>           verticies;

    SpaceGroup():
      loaded(false),
      name("noName"),
      centered(false),
      rating(-1),
      selected(false)
    {}

    SpaceGroup(const SpaceGroupType& sg):
      loaded  (sg.loaded),
      number  (sg.number),
      name    (sg.name),
      centered(sg.centered),
      rating  (-1),
      selected(false),
      operationIds (sg.operationIds)
    {}
    
    SpaceGroup& operator = (const SpaceGroup& sg) {
      loaded         = sg.loaded;
      name           = sg.name;
      centered       = sg.centered;
      number         = sg.number;
      rating         = sg.rating;
      selected       = sg.selected;
      operationIds   = sg.operationIds;
      return (*this);
    }

    template<template<typename, size_t, typename> class LatticeTemplate>
    std::vector<CartesianPosition<Field,DIM> >  cartesianAsymmetricUnitVerticies(const LatticeTemplate<Field,DIM,Algorithms>& lat) const {
      std::vector<CartesianPosition<Field,DIM> > result;
      for (size_t i=0; i < verticies.size(); i++) 
	result.push_back(lat.cartesian(verticies[i]));
      return result;
    }

    template<template<typename, size_t, typename> class LatticeTemplate>
    void initialize(size_t number_, 
		    const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry) {
      number = number_;
      loadSpaceGroupData();
      loaded = true;
      rating = getRating(symmetry);
    }

    size_t numOperations() const {
      return operationIds.size();
    }

    template<size_t SGNUM, size_t VARIENT>
    void loadGroupData() {

      typedef SpaceGroupData<Field,DIM,SGNUM,VARIENT,Algorithms> SpaceGroupDataType;

      SpaceGroupDataType::operations(operationIds);

      if (SpaceGroupDataType::Centered)
	centered = true;
      
      verticies = SpaceGroupDataType::verticies();
      name      = SpaceGroupDataType::name();

    }

    void loadSpaceGroupData() {

      switch (number) {
      case 1:  {loadGroupData<1, 1>(); break;}
      case 2:  {loadGroupData<2, 1>(); break;}
      case 3:  {loadGroupData<3, 1>(); break;}
      case 4:  {loadGroupData<4, 1>(); break;}
      case 5:  {loadGroupData<5, 1>(); break;}
      case 6:  {loadGroupData<6, 1>(); break;}
      case 7:  {loadGroupData<7, 1>(); break;}
      case 8:  {loadGroupData<8, 1>(); break;}
      case 9:  {loadGroupData<9, 1>(); break;}
      case 10: {loadGroupData<10,1>(); break;}
      case 11: {loadGroupData<11,1>(); break;}
      case 12: {loadGroupData<12,1>(); break;}
      case 13: {loadGroupData<13,1>(); break;}
      case 14: {loadGroupData<14,1>(); break;}
      case 15: {loadGroupData<15,1>(); break;}
      case 16: {loadGroupData<16,1>(); break;}
      case 17: {loadGroupData<17,1>(); break;}
      default: {
	throw std::range_error("loadSpaceGroupData recieved an out-of-range index");
      }
      }
    }

    template<template<typename, size_t, typename> class LatticeTemplate>
    Field getRating( const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry) {
      
      typedef LatticeTemplate<Field,DIM,Algorithms>                                         LatticeType;
      typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms>             AppliedSymmetryElementType;

      static Field zero(0);
      Field maxDistance = zero;

      //for (int i=0; i< symmetry.allAppliedSymmetryElements.size(); i++) {

      for (size_t i=0; i< operationIds.size(); i++) {
	
	int id = operationIds[i];

	const AppliedSymmetryElementType& op = symmetry.allAppliedSymmetryElements[id];
	
	if (op.distance > maxDistance)
	  maxDistance = op.distance;
      }
      
      return maxDistance;
    }

    template<template<typename, size_t, typename> class LatticeTemplate>
    AppliedSymmetryElementIterator<Field,DIM,Occupant,LatticeTemplate,Algorithms>
    beginOperations(const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry) {

      typedef AppliedSymmetryElementIterator<Field,DIM,Occupant,LatticeTemplate,Algorithms> AppliedSymmetryElementIteratorType;
      return AppliedSymmetryElementIteratorType(*this,symmetry.allAppliedSymmetryElements);
    }

    template<template<typename, size_t, typename> class LatticeTemplate>
    AppliedSymmetryElementIterator<Field,DIM,Occupant,LatticeTemplate,Algorithms>
    endOperations(const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry) {

      typedef AppliedSymmetryElementIterator<Field,DIM,Occupant,LatticeTemplate,Algorithms> AppliedSymmetryElementIteratorType;
      return AppliedSymmetryElementIteratorType::end(*this,symmetry.allAppliedSymmetryElements);
    }
    
    
  };

  //======================================================================

  /** \ingroup XML
   *
   * SpaceGroup  XML output
   **/
  template<typename Field, size_t DIM, typename Occupant, typename LatticeType, typename Algorithms> 
  Tag asymmetricUnitToXML(const SpaceGroup<Field,DIM,Occupant,Algorithms>& spaceGroup,
			  const LatticeWithPattern<Field,DIM,Occupant,LatticeType, Algorithms>& latpat) {

    const LatticeType& lat = latpat;
    std::vector< CartesianPosition<Field,DIM> >  verticies = spaceGroup.cartesianAsymmetricUnitVerticies(lat);
    
    Tag result("AsymmetricUnit");

    std::ostringstream buff;
    buff << "[";
    for (size_t i=0; i< verticies.size(); i++) {
      if (i!=0) buff << ",";
      buff << "(";
      for (size_t j=0; j<DIM; j++) {
	if (j!=0) buff << ",";
	buff << verticies[i][j];
      }
      buff << ")";
    }
    buff << "]";
    result["verticies"] = buff.str();
    return result;
  }

  //======================================================================

  /** \ingroup XML
   *
   * SpaceGroup  XML output
   **/
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms> 
  Tag toXML(const SpaceGroup<Field,DIM,Occupant,Algorithms>& spaceGroup,
	    std::string name="SpaceGroup") {
    
    Tag result(name);
    result["name"]     = spaceGroup.name;
    result["DIM"]      = DIM;
    result["number"]   = spaceGroup.number;
    result["numOps"]   = spaceGroup.numOperations();
    result["centered"] = spaceGroup.centered;
    result["rating"]   = spaceGroup.rating;
    result["selected"] = spaceGroup.selected;

    std::ostringstream buff;
    buff << "[";
    for (size_t i=0; i< spaceGroup.operationIds.size(); i++) {
      if (i!=0) buff << ",";
      buff << spaceGroup.operationIds[i];
    }
    buff << "]";
    result["operationIds"] = buff.str();
    return result;
  }

  //======================================================================

  /** \ingroup XML
   *
   * SpaceGroup  XML output
   **/
  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   template<typename, size_t, typename> class LatticeTemplate, 
	   typename Algorithms>
  Tag toXML(const SpaceGroup<Field,DIM,Occupant,Algorithms>& spaceGroup,
	    const Symmetry<Field,DIM,Occupant,LatticeTemplate,Algorithms>& symmetry,
	    std::string name="SpaceGroup") {
    
    Tag result = toXML(spaceGroup,name);
    
    if (spaceGroup.loaded) {
      
      typedef LatticeTemplate<Field,DIM,Algorithms>                             LatticeType;
      typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
      
      for (size_t i=0; i< symmetry.allAppliedSymmetryElements.size(); i++) {
	if (spaceGroup.opFlags != 1) continue;
	
	const AppliedSymmetryElementType& op = symmetry.allAppliedSymmetryElements[i];
	result.add(toXML(*op));
      }
    }
    return result;
  }

   //======================================================================

  template<typename Field, typename Occupant, typename LatticeType, typename Algorithms> 
  std::ostream& operator << (std::ostream& os, SpaceGroup<Field,2,Occupant,Algorithms> spaceGroup) {
    
    os << "2DSpaceGroup[" << spaceGroup.number << "]: \"" << spaceGroup.nam << "\" " 
       << "{numOps=" << spaceGroup.numOperations 
       << " numGens=" << spaceGroup.numGenerators
       << " centered=" << spaceGroup.centered << "}" << std::endl;
    
    return os;
  }

} /** namespace spimag **/

#endif // Psimag_SpaceGroup2D

/*@}*/
