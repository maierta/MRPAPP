//-*-C++-*-       

#ifndef FETCHGEOMETRY_HEADER_H
#define  FETCHGEOMETRY_HEADER_H

#include <iostream>                                  
#include <iomanip>

#include <sstream>                                  
#include <limits>      
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include "GroupActionFilter.h"
#include "DcaCrystal.h"
#include "BasicCrystalAlgorithms.h"

#include "CellParameters.h"
#include "Crystal.h"
#include "SuperCrystal.h"
#include "Occupant.h"
#include "SuperCrystal.h"
#include "SuperCrystalBuilder.h"
#include "PatternData.h"
#include "Lattice.h"
#include "Orbits.h"


namespace rpa {

  //====================================================================== fetchGeometry
  //======================================================================

  template<typename SymmetryType,
	   template<typename> class MatrixTemplate>
  void fetchGroupAction(const SymmetryType&        symmetry,
			MatrixTemplate<int>&       target_data, 
			GroupActionFilter::Type    filterType) 
  {
    switch (filterType) {
      
      // Use the full group action table
    case (GroupActionFilter::None): {
      target_data = symmetry.symmetryGroup.groupAction;
      return;
    }
      // Use the full group action table
    case (GroupActionFilter::Full): {
      target_data = symmetry.symmetryGroup.groupAction;
      return;
    }
      // Use just the point group operations from the group action table
      //
    case (GroupActionFilter::PointGroupOnly): {
      target_data = symmetry.pointGroup.groupAction;
      return;
    }
    case (GroupActionFilter::Unique): {
      symmetry.symmetryGroup.groupAction.uniqueActions(target_data);
      return;
    }
    }
    throw std::logic_error("In fetchGroupAction, unknown filterType given");
  }

  //----------------------------------------------------------------------

  //fetchGeometry2d(4,2,0,4,
  //	1,1,90,
  template<typename Field,
	   template<typename> class MatrixTemplate,
	   typename VectorLikeType>
  void fetchGeometry2D(Field latticeA_length, 
		       Field latticeB_length, 
		       Field latticeAngle,
		       int superlatticeA_x,
		       int superlatticeA_y,
		       int superlatticeB_x,
		       int superlatticeB_y,
		       GroupActionFilter::Type  filterType,
		       Crystal<Field,MatrixTemplate>& rCluster,
		       const VectorLikeType&          phase,
		       Crystal<Field,MatrixTemplate>& kCluster,
		       int                            parallelProcessingId=-911)
  {

    enum{DIM=2};

    typedef psimag::BasicCrystalAlgorithms     Algorithms;
  
    typedef psimag::CellParameters<Field,DIM,Algorithms>                            CellParametersType;
    typedef psimag::Lattice<Field,DIM,Algorithms>                                   LatticeType;
    typedef psimag::Crystal<Field,DIM, Occupant, Algorithms>                        CrystalType;
    typedef psimag::SuperCrystal<Field,DIM, Occupant, Algorithms>                   SuperCrystalType;
    typedef psimag::SuperCrystalBuilder<Field,DIM,Occupant,Algorithms>              SuperCrystalBuilderType;
    typedef psimag::PatternData<Field,DIM,Occupant,Algorithms>                      PatternDataType;
  
    //======================================== Build the crystal lattice
    LatticeType lat(CellParametersType(latticeA_length, latticeB_length, latticeAngle));

    //======================================== Build the crystal pattern
    PatternDataType patternData;
    patternData.insert( Occupant("atom1","yellow"), cellPosition<Field,Algorithms>(0.0,0.0) );
    CrystalType crystal(lat, patternData);
  
    //======================================== Build the K-point pattern
    PatternDataType rpatternData;
    rpatternData.insert( Occupant("kpoint1","red"),  cellPosition<Field,Algorithms>(phase[0],phase[1]) );
    //rpatternData.insert( Occupant("kpoint2","pink"), cellPosition<Field,Algorithms>(0.5,0.5) );

    //======================================== Build the SuperCrystal
    SuperCrystalBuilderType builder(superlatticeA_x, superlatticeA_y, superlatticeB_x, superlatticeB_y, rpatternData);
    SuperCrystalType superCryst(crystal, builder);
  
    //======================================== Write four XML files describing the supercrystal in detail

    if (parallelProcessingId == 0)
      superCryst.dumpXML(""); //xmlOutput/");

    //======================================== Fetch the rCluster lattice 
    lat.getBasisVectors(rCluster.latticeVectors);
  
    //======================================== Fetch the real sites cartesian coordinates 
    superCryst.pattern.setCartesianSites(rCluster.sites);
  
    //======================================== Fetch the kCluster lattice 
    superCryst.reciprocalCrystal.getBasisVectors(kCluster.latticeVectors);
  
    //======================================== Fetch the k sites cartesian coordinates 
    superCryst.reciprocalCrystal.pattern.setCartesianSites(kCluster.sites);
  
  
    //======================================== Fetch the real translation difference matricies
    superCryst.pattern.buildDiffIndex(rCluster.siteDiffToSite);

    //======================================== Fetch the reciprocal translation difference matricies
    superCryst.reciprocalCrystal.pattern.buildDiffIndex(kCluster.siteDiffToSite);

    //======================================== Fetch the reciprocal translation addition matricies
    superCryst.reciprocalCrystal.pattern.buildPlusIndex(kCluster.siteSumToSite);
  
    //======================================== Fetch the real group action matricies
    fetchGroupAction(superCryst.symmetry, rCluster.symmetryGroupAction, filterType);

    //======================================== Fetch the reciprocal group action matricies
    fetchGroupAction(superCryst.reciprocalCrystal.symmetry, kCluster.symmetryGroupAction, filterType);

    if (filterType == GroupActionFilter::PointGroupOnly) {
      //======================================== Fetch the real orbits
      rCluster.orbits = superCryst.symmetry.pointGroup.orbits;
    
      //======================================== Fetch the reciprocal orbits
      kCluster.orbits = superCryst.reciprocalCrystal.symmetry.pointGroup.orbits;
    }
    else {
      //======================================== Fetch the real orbits
      rCluster.orbits = superCryst.symmetry.symmetryGroup.orbits;
    
      //======================================== Fetch the reciprocal orbits
      kCluster.orbits = superCryst.reciprocalCrystal.symmetry.symmetryGroup.orbits;
    }
  }

  //======================================================================
  //======================================================================

  template<typename Field, 
	   size_t DIM, 
	   typename Algorithms>
  void loadPatternData(int numPoints, 
		       psimag::PatternData<Field,DIM,Occupant,Algorithms>& patternData ) {
    
    size_t nSide( static_cast<size_t>(sqrt(static_cast<Field>(numPoints))));

    Field  nSideF(nSide);
    
    for(size_t i=0; i<nSide;i++) 
      for(size_t j=0; j<nSide;j++) 

	patternData.insert( Occupant("mesh","blue"),   cellPosition<Field,Algorithms>(static_cast<Field>(i)/nSideF,
										      static_cast<Field>(j)/nSideF) );
  }

  // mesh fetch geometry
  template<typename Field,
	   template<typename> class MatrixTemplate>
  void fetchMeshGeometry2D(Field latticeA_length, 
				Field latticeB_length, 
				Field latticeAngle,
				int superlatticeA_x,
				int superlatticeA_y,
				int superlatticeB_x,
				int superlatticeB_y,
				int numPoints,
				MatrixTemplate<Field>& sitesKmesh)
  {
    enum{DIM=2};
    typedef psimag::BasicCrystalAlgorithms                                          Algorithms;
    typedef psimag::LatticeCoordinates<DIM>                                         SubLatticeVecCoordType;
    typedef std::vector<SubLatticeVecCoordType>                                     SubLatticeVecCoordsType;
    typedef psimag::CellParameters<Field,2,Algorithms>                              CellParametersType;
    typedef psimag::PatternData<Field,DIM,Occupant,Algorithms>                      PatternDataType;
    typedef psimag::Lattice<Field,DIM,Algorithms>                                   LatticeType;
    typedef psimag::LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>   LatticeWithPatternType;
    typedef psimag::FloodTiler<Field,DIM,Occupant,Lattice,Algorithms>               FloodTilerType;
    typedef psimag::ReciprocalLattice<Field,DIM,Lattice,Algorithms>                 ReciprocalLatticeType;

    typedef psimag::SuperCrystalBuilder<Field,DIM,Occupant,Algorithms>              SuperCrystalBuilderType;

    /** ====================================================================== Construct the given cell */

    PatternDataType patternData;
    loadPatternData(numPoints, patternData);

    LatticeType lat(CellParametersType(latticeA_length, latticeB_length, latticeAngle));

    LatticeWithPatternType crystalLatPat(lat, patternData);
  
    /** ====================================================================== Construct the k-point tile */
    PatternDataType rpatternData;
    loadPatternData(numPoints, rpatternData);

    /** ====================================================================== Construct a helper class */
    SuperCrystalBuilderType builder(superlatticeA_x,superlatticeA_y,superlatticeB_x,superlatticeB_y, rpatternData);
    SubLatticeVecCoordsType subLatVecCoords(builder.getSubLatticeVecCoords());

    /** ================================================== Construct the superCrystal's LatticeWithPattern  */
    LatticeWithPatternType superCrystalLatPat      (FloodTilerType(subLatVecCoords,crystalLatPat).getTiledLatticeWithPattern());

    /** ================================================== Construct the reciprocal lattice of the crystal that this is a SuperCrystal of. */
    //LatticeType            reciprocalCrystalLattice(ReciprocalLatticeType(crystalLatPat).getLattice());
  
    /**  ================================================== Construct the reciprocal lattice of this SuperCrystal's lattice. */
    LatticeType            reciprocalSuperLattice (ReciprocalLatticeType(superCrystalLatPat).getLattice()); 

    /**  ================================================== Construct the reciprocal  supe LatticeWithPattern. */
    LatticeWithPatternType reciprocalSuperLatPat  (reciprocalSuperLattice, builder.getReciprocalSuperCrystalPattern());

    /**  ================================================== Construct the reciprocal crystal's LatticeWithPattern. */
    //LatticeWithPatternType reciprocalCrystalLatPat(FloodTilerType(reciprocalCrystalLattice, reciprocalSuperLatPat)
    //						   .getTiledLatticeWithPattern());

    //reciprocalCrystalLatPat.pattern.loadCartesian(sitesKmesh);
    reciprocalSuperLatPat.pattern.loadCartesian(sitesKmesh);

    //  reciprocalSuperLatPat.pattern.loadCartesian(sitesRmesh);
  
  }

  template<typename Field>
  void dumpMeshXML(const std::string outputDir, const std::string name, const Matrix<Field>& mesh)  {
    
    std::ofstream xmlFile;

    xmlFile.open((outputDir+name+".xml").c_str());
    xmlFile << XMLHeading(XMLHeading::XHTML) << std::endl;
    xmlFile << toXML(mesh)  << std::endl;
    xmlFile.close();

  }


} // rpa namespace
  
#endif
		   
