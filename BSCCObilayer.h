
#ifndef BSCCOBILAYER_H
#define BSCCOBILAYER_H


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
	class model {
	private:
		typedef MatrixTemplate<Field> 		MatrixType;
		typedef std::complex<Field>			ComplexType;
		typedef MatrixTemplate<ComplexType> ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 						FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		Field t,tp,tpp,tperp,mu;

	public:
		FieldType nbands;
		ComplexMatrixType spinMatrix;
		ComplexMatrixType chargeMatrix;


		model(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb),
			spinMatrix(nbands*nbands,nbands*nbands),
			chargeMatrix(nbands*nbands,nbands*nbands)
		{
			if (nbands != 2) 
				std::cerr<<"Number of orbitals should be 2! Exiting ...\n";

			setupInteractionMatrix();

			if      (param.Case == "BSCCObilayer_UD30") {
			  	t = 0.390; tp = 0.120; tpp=0.045; tperp=0.108; mu = 0.340;
			}
			else if (param.Case == "BSCCObilayer_OD") {
			  	t = 0.360; tp = 0.108; tpp=0.036; tperp=0.108; mu = -param.mu;
			}
		}

		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {

			FieldType cx,cy,c2x,c2y;
			cx = cos(k[0]); cy = cos(k[1]); 
			c2x = cos(2*k[0]); c2y = cos(2*k[1]);

			FieldType ek   = -2*t*(cx + cy) + 4*tp*cx*cy - 2*tpp*(c2x+c2y) + mu;
			FieldType ekz  = tperp/4. * pow((cx-cy),2);

			FieldType r1(1./sqrt(2.));

// #ifdef USE_ABANDONLY
// 			eigenvals[0] = ek + ekz;
// 			eigenvects(0,0) =  1.0; 
// #else			
			eigenvals[0] = ek - ekz;
			eigenvals[1] = ek + ekz;
			eigenvects(0,0) =  r1; eigenvects(0,1) =  r1;
			eigenvects(1,0) =  r1; eigenvects(1,1) = -r1;
// #endif
		}

		void setupInteractionMatrix() {
			size_t nOrb(param.nOrb);
			FieldType U(param.U);

			for (size_t l1 = 0; l1 < nOrb; ++l1) {
				size_t ind1 = l1+l1*nOrb;
				spinMatrix(ind1,ind1)   = U+param.deltaU[l1];
				chargeMatrix(ind1,ind1) = -U-param.deltaU[l1];;
				}	
			}


		std::complex<Field> calcSus(const ComplexMatrixType& sus, const std::string& component = "zz") const {
			std::complex<Field> chiPhys(0.0);
			for (size_t l1 = 0; l1 < param.nOrb; ++l1) {
				for (size_t l2 = 0; l2 < param.nOrb; ++l2){
					size_t ind1(l1+l1*param.nOrb);
					size_t ind2(l2+l2*param.nOrb);
					chiPhys += 0.5*sus(ind1,ind2);
				}
			}
			return chiPhys;
		}

	};

}

#endif

