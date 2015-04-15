#ifndef BILAYER_H
#define BILAYER_H


#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib> // for atof and atoi

#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "utilities.h"
#include "Range.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class bilayer {
	private:
		typedef MatrixTemplate<Field> 		MatrixType;
		typedef std::complex<Field>			ComplexType;
		typedef MatrixTemplate<ComplexType> ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 						FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		
	public:
		FieldType nbands;


		bilayer(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb)
		{
		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		  	FieldType t,tp,tpp,tperp;

		  	// t = 0.360; tp = 0.3*t; tpp = 0.15*t; tperp = 0.135;
		  	// t = 0.360; tp = 0.3*t; tpp = 0.0*t; tperp = 0.135/4;
		  	// t = 0.180; tp = 0.3*t; tpp = 0.15*t; tperp = 0.026;
		  	// t = 0.180; tp = 0.2*t; tpp = 0.0*t; tperp = 0.026;
		  	t = 1.0; tp = 0.0; tpp = 0.0; tperp = 0.10*t;

			FieldType cx,cy,cz,c2x,c2y;
			cx = cos(k[0]); cy = cos(k[1]); cz = cos(k[2]);
			c2x = cos(2*k[0]); c2y = cos(2*k[1]);
			
			FieldType ek   = -2*t*(cx+cy) + 4*tp*cx*cy - 2*tpp*(c2x+c2y) - param.mu;
			FieldType ekz  = tperp*pow((cx-cy),2)*cz;   
			// FieldType ekz  = -tperp*cz;   

			eigenvals[0] = ek + ekz;

			eigenvects(0,0) =  1.0;
		}

	};

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class orthoIIBilayer {
	private:
		typedef MatrixTemplate<Field> 		MatrixType;
		typedef std::complex<Field>			ComplexType;
		typedef MatrixTemplate<ComplexType> ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 						FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		
	public:
		FieldType nbands;


		orthoIIBilayer(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb)
		{
		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		  	FieldType t,tp,tpp,V,tperp;

		  	// t = 0.558; tp = 0.49*t; tpp=0.5*tp; V = 0.075; tperp=0.2*t;
		  	t = 0.558; tp = 0.49*t; tpp=0.5*tp; V = 0.075; tperp=0.0;

			FieldType cx,cy,cz,c2x,c2y;
			cx = cos(k[0]); cy = cos(k[1]); cz = cos(k[2]);
			c2x = cos(2*k[0]); c2y = cos(2*k[1]);
			
			FieldType ek   = -2*t*cy - 2*tpp*(c2x+c2y) - param.mu;
			// FieldType ekz  = tperp*pow((cx-cy),2)*cz;   
			FieldType ekz  = -tperp*cz;   
			FieldType ekOrtho = sqrt(4.*pow(cx,2)*pow((t-2.*tp*cy),2) + pow(V,2)/4.0);

			eigenvals[0] = ek + ekz - ekOrtho;
			eigenvals[1] = ek + ekz + ekOrtho;

			FieldType r1(1./sqrt(2.));
			eigenvects(0,0) =  r1; eigenvects(0,1) = r1;
			eigenvects(1,0) =  r1; eigenvects(1,1) =-r1;
		}

	};

}

#endif
