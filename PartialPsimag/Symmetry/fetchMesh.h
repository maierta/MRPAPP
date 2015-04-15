//-*-C++-*-       
#ifndef FETCHMESH_HEADER_H
#define FETCHMESH_HEADER_H
#include <iostream>                                  
#include <iomanip>

#include <sstream>                                  
#include <limits>      
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include "Matrix.h"
#include "Tag.h"

#include "BasicCrystalAlgorithms.h"

#include "CellPosition.h"

#include "Occupant.h"
#include "PatternData.h"
#include "Orbits.h"

#include "CellParameters.h"
#include "Lattice.h"

#include "Crystal.h"
#include "SuperCrystal.h"
#include "SuperCrystal.h"
#include "SuperCrystalBuilder.h"

using namespace psimag;

namespace rpa {

  
  template<typename Field, size_t DIM, typename Algorithms>
  void loadPatternData(int numPoints, PatternData<Field,DIM,Occupant,Algorithms>& patternData ) {
    
    size_t nSide( static_cast<size_t>(sqrt(static_cast<Field>(numPoints))));

    Field  nSideF(nSide);
    
    for(size_t i=0; i<nSide;i++) 
      for(size_t j=0; j<nSide;j++) 

	patternData.insert( Occupant("mesh","blue"),   cellPosition<Field,Algorithms>(static_cast<Field>(i)/nSideF,
										       static_cast<Field>(j)/nSideF) );
  }

  // mesh fetch geometry
  template<typename Field,typename Algorithms>
  void fetchMeshGeometry2D(Field latticeA_length, Field latticeB_length, Field latticeAngle,
			   int superlatticeA_x,int superlatticeA_y,int superlatticeB_x,int superlatticeB_y,
			   int numPoints,
			   psimag::Matrix<Field>& sitesKmesh)
  {
    enum{DIM=2};
    typedef LatticeCoordinates<DIM>                                         SubLatticeVecCoordType;
    typedef std::vector<SubLatticeVecCoordType>                             SubLatticeVecCoordsType;
    typedef CellParameters<Field,2,Algorithms>                              CellParametersType;
    typedef PatternData<Field,DIM,Occupant,Algorithms>                      PatternDataType;
    typedef Lattice<Field,DIM,Algorithms>                                   LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>   LatticeWithPatternType;
    typedef FloodTiler<Field,DIM,Occupant,Lattice,Algorithms>               FloodTilerType;
    typedef ReciprocalLattice<Field,DIM,Lattice,Algorithms>                 ReciprocalLatticeType;

    typedef SuperCrystalBuilder<Field,DIM,Occupant,Algorithms>              SuperCrystalBuilderType;

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
    /*
   for (size_t i=0;i<sitesKmesh.n_row();i++) {
   	for (size_t j=0;j<sitesKmesh.n_col();j++) {
	if (j==0) 
    	sitesKmesh(i,j) = sitesKmesh(i,j) * Field(sqrt(numPoints))/Field(sqrt(numPoints)-1.0) - M_PI/double(superlatticeA_x);
	else
	sitesKmesh(i,j) = -sitesKmesh(i,j) * Field(sqrt(numPoints))/Field(sqrt(numPoints)-1.0) - M_PI/double(superlatticeA_x);
	
}} */
#ifndef USE_MPI    
    dumpMeshXML("xmlOutput/", "K_Mesh", sitesKmesh);
#endif    
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

    
  template<typename Field,typename Algorithms>
  void fetchMeshGeometry1D(int length1d,int numPoints, psimag::Matrix<Field>& sitesKmesh)
  {
    enum{DIM=1};
    //int tmp;
    std::vector<Field> tmpvector;
    tmpvector.resize(1);

    sitesKmesh.resize(numPoints,DIM);
    //======================================== Fetch the k sites cartesian coordinates 

    for (int i=0;i<numPoints;i++) 
      for (int j=0;j<DIM;j++) 
	sitesKmesh(i,j) =  2.0*M_PI*i/((numPoints-1.0)*length1d) -M_PI/Field(length1d);
  
  }

  //======================================================================
  template<typename Field,typename Algorithms>
  void fetchMeshGeometry(int dim,
			 Field latticeA_length, Field latticeB_length, Field latticeAngle,
			 int superlatticeA_x,int superlatticeA_y,int superlatticeB_x,int superlatticeB_y,
			 int numPoints,
			 psimag::Matrix<Field>& sitesKmesh)
  {
    if (dim==1) {
      fetchMeshGeometry1D<Field,Algorithms>(superlatticeA_x,numPoints, sitesKmesh);
      return;
    }
    if (dim==2) {
      fetchMeshGeometry2D<Field,Algorithms>(latticeA_length,  latticeB_length,  latticeAngle,
					    superlatticeA_x,  superlatticeA_y, superlatticeB_x, superlatticeB_y,
					    numPoints, sitesKmesh);
      return;
    }
    std::ostringstream buff;
    buff << "error in fetchGeometry! dim = " << dim << " not supported!";
    throw std::logic_error(buff.str());
  }    

}  
#endif
		   
