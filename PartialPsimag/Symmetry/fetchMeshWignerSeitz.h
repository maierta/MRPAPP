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

#include "CellParameters.h"
#include "Crystal.h"
#include "SuperCrystal.h"
#include "Occupant.h"
#include "BasicCrystalAlgorithms.h"
#include "XMLHeading.h"
#include "SuperCrystal.h"
#include "SuperCrystalBuilder.h"
#include "PatternData.h"
#include "CellPosition.h"
#include "Matrix.h"
#include "Lattice.h"
#include "Orbits.h"
#include "Tag.h"

using namespace psimag;

namespace rpa {

  //! generates a cubic mesh that is shifted 
  template<typename Field>
  void generateMeshPoints(int                    dim,
			  psimag::Matrix<Field>& meshPoints, 
			  size_t                 numPoints)
  {
    int x,y;
    std::ostringstream buff;
    int nx;
    Field shift;
	
    meshPoints.resize(numPoints,dim);
    for (size_t i=0;i<numPoints;i++) {
      if (dim==1) {
	nx = numPoints;
	shift = -M_PI+M_PI/static_cast<Field>(nx);
	meshPoints(i,0)=2.0*M_PI*i/static_cast<Field>(nx)+shift;
      } 
      else 
	if (dim==2) {
	  buff<<__FILE__<<": generateMeshPoints: numPoints="<<numPoints<<" is not a perfect square.\n";
	  if (!isPerfectSquare(numPoints,nx)) throw logic_error(buff.str());
	  int ny=nx;
	  y=int(i/nx);
	  x=i-y*nx;
	  shift = -M_PI+M_PI/static_cast<Field>(nx);
	  meshPoints(i,0)=2.0*M_PI*x/static_cast<Field>(nx)+shift;
	  meshPoints(i,1)=2.0*M_PI*y/static_cast<Field>(ny)+shift;
	} 
	else {
			
	  buff<<__FILE__<<": generateMeshPoints: dimension="<<dim<<" not implemented.\n"; 
	  throw std::logic_error(buff.str());
	}
    }
    
  }

  //! Distance between vectors v1 and v2 modulo 2pi
  template<class Field>
  Field meshDistance(std::vector<Field> const &v1,std::vector<Field> const &v2)
  {
    Field sum=0;
    std::vector<Field> v(v1.size());
    Field tmp;
    for (size_t i=0;i<v1.size();i++) {
      tmp = v1[i]-v2[i];
      if (tmp<-M_PI) tmp+=2.0*M_PI;
      if (tmp>M_PI) tmp-=2.0*M_PI;
      sum += tmp*tmp;
    }
    return sum;
  }

  //! Resturns true is thisMeshPoint is closest to the 0 k-vector than to any other k-vector
  //! and false otherwise. 
  template<class Field>
  bool isClosestToZero(std::vector<Field> const &thisMeshPoint,psimag::Matrix<Field> const &sitesK,int i)
  {
    size_t j;
    std::vector<Field> v(thisMeshPoint.size());
    for (j=0;j<v.size();j++) v[j]=0;
    Field d0=meshDistance(thisMeshPoint,v);
    Field tmp;
	
    for (j=0;j<sitesK.n_row();j++) {
      sitesK.getRow(j,v);
      tmp = meshDistance(thisMeshPoint,v);

      if (d0-tmp>1e-8)  {// if distance to zero > distance to other points,
	return false; // then return false
      }
    }
    return true;	
  }

  //! Filters meshPoints into sitesKmesh such that only the mesh points closests to sitesK==0 will be considered
  template<typename Field>
  void filterMeshPoints(psimag::Matrix<Field>& sitesKmesh,psimag::Matrix<Field> const &meshPoints,psimag::Matrix<Field> const &sitesK)
  {
    //int counter=0;
    size_t i;
    std::vector<int> array;
    std::vector<Field> v(meshPoints.n_col());
	
    for (i=0;i<meshPoints.n_row();i++) { // mesh point
      meshPoints.getRow(i,v);

      if (isClosestToZero(v,sitesK,i)) { // if it is closest to zero
	array.push_back(i); // then add this point
      }
    }
    if (array.size()<=0) {
      throw logic_error("filterMeshPoints: There are no points that survived the filtering process.\n");
    }
    size_t d;
    sitesKmesh.resize(array.size(),meshPoints.n_col());
    for (i=0;i<array.size();i++) {
      for (d=0;d<meshPoints.n_col();d++) sitesKmesh(i,d)=meshPoints(array[i],d);
    }
  }

  //! Sets sitesKmesh based on a minimum distance algorithm
  template<typename Field>
  void fetchMeshGeometry(int dim,
			 int numPoints,
			 psimag::Matrix<Field> const &sitesK,
			 psimag::Matrix<Field>& sitesKmesh)
  {
    psimag::Matrix<Field> meshPoints;
    generateMeshPoints(dim,meshPoints,numPoints);
    filterMeshPoints(sitesKmesh,meshPoints,sitesK);
	
  }

} // namespace rpa  
#endif
		   
