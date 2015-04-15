//-*-C++-*-       

#ifndef  FETCH_GEOMETRY_1D_HEADER_H
#define  FETCH_GEOMETRY_1D_HEADER_H

//! \file fetchGeometry1D.h
//! \brief Fetch a 1D Geometry  and its Symmetry

#include <iostream>                                  
#include <iomanip>

#include <sstream>                                  
#include <limits>      
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include "DcaCrystal.h"

namespace rpa {
  using namespace psimag;

  
  //
  template<typename Field, 
	   template<typename> class MatrixTemplate>
  void fetchGeometry1D(int length1d,
		       Crystal<Field,MatrixTemplate>& rCrystal,
		       Crystal<Field,MatrixTemplate>& kCrystal)
  {
    enum{ DIM=1 };

    //======================================== Fetch the real sites cartesian coordinates 
    
    rCrystal.sites.resize(length1d,DIM);
    for (int i=0;i<length1d;i++) {
      for (int j=0;j<DIM;j++) 
	rCrystal.sites(i,j) = i;

    }			
  
    //======================================== Fetch the k sites cartesian coordinates 
    kCrystal.sites.resize(length1d,DIM); 
    for (int i=0;i<length1d;i++) 
      for (int j=0;j<DIM;j++) 
	kCrystal.sites(i,j) = (2.0*M_PI*i/length1d -M_PI);
  
    //======================================== Fetch the real translation difference matricies
    rCrystal.siteDiffToSite.resize(length1d,length1d);
    for (int i=0;i<length1d;i++) {
      for (int j=0;j<length1d;j++) {
	int tmp = i-j;
	rCrystal.siteDiffToSite(i,j)=(tmp<0) ? tmp + length1d : tmp; 
      }
    }

    //======================================== Fetch the reciprocal translation difference matricies

    kCrystal.siteDiffToSite.resize(length1d,length1d);

    for (int i=0;i<length1d;i++) {
      for (int j=0;j<length1d;j++) {
	int tmp = i-j +int(length1d/2.0);
	if (tmp<0) 
	  tmp += length1d;
	if (tmp >= length1d) 
	  tmp-= length1d;
	kCrystal.siteDiffToSite(i,j)= tmp;
      }
    }

    //======================================== Fetch the reciprocal translation addition matricies
    kCrystal.siteSumToSite.resize(length1d,length1d);
    for (int i=0;i<length1d;i++) {
      for (int j=0;j<length1d;j++) {
	int tmp = i+j -int(length1d/2.0);
	if (tmp<0) tmp += length1d;
	if (tmp>=length1d) tmp-= length1d;
	kCrystal.siteSumToSite(i,j)=tmp;
      }
    }
  
    //======================================== Fetch the real group action matricies
    int nGroup=2;
    rCrystal.symmetryGroupAction.resize(nGroup,length1d);
    for (int i=0;i<length1d;i++) {
      for (int j=0;j<nGroup;j++) {
	int tmp;
	if  (j==0) 
	  tmp=i; // identity
	else {// inversion
	  tmp = length1d-i;
	  if (tmp>=length1d) tmp-=length1d;
	}
	rCrystal.symmetryGroupAction(j,i)=tmp;
      }
    }

    //======================================== Fetch the reciprocal group action matricies
    kCrystal.symmetryGroupAction.resize(nGroup,length1d);
    for (int i=0;i<length1d;i++) {
      for (int j=0;j<nGroup;j++) {
	int tmp;
	if  (j==0) 
	  tmp=i; // identity
	else {// inversion
	  tmp = length1d-i;
	  if (tmp>=length1d) tmp-=length1d;
	}
	kCrystal.symmetryGroupAction(j,i)=tmp;
      }
    }
  }

  template<typename Field, template<typename> class MatrixTemplate>
  void fetchMeshGeometry1D(int length1d,
			   int numPoints, 
			   MatrixTemplate<Field>& sitesKmesh)
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

} // rpa namespace

#endif
		   
