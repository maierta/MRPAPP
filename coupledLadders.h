// Model file for coupled ladders, 2 bands total
#ifndef COUPLEDLADDERS_H
#define COUPLEDLADDERS_H

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
		typedef std::complex<Field>		    ComplexType;
		typedef MatrixTemplate<ComplexType>	ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 				        FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		
	private:
		size_t msize;
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

			msize = 4;
			setupInteractionMatrix();

		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		    	FieldType t, tl, tp, tlp, cy, cx, c2x, lambdak;
			ComplexType hABk, ekx, emkx;
			const ComplexType ii = ComplexType(0.0,1.0);

			t = 1.0; tp = -0.3; 
			tl  = t  * param.couplingRatio;
			tlp = tp * param.couplingRatio;

			cx = cos(k[0]); cy = cos(k[1]); c2x = cos(2.0*k[0]);
			ekx = exp(ii*k[0]); emkx = conj(ekx);

			hABk     = -t*ekx - tl*emkx - 2.*cy*(tp*ekx + tlp*emkx);
			lambdak  = pow(t,2) + pow(tl,2) + 2.*t*tl*c2x;
			lambdak += 4.*cy*(t*tp + tl*tlp+(t*tlp+tp*tl)*c2x);
			lambdak += 4.*pow(cy,2)*(pow(tp,2)+pow(tlp,2)+2.*tp*tlp*c2x);
			lambdak = sqrt(lambdak);

			// eigen(evals,temp);
			eigenvals[0] = -2*t*cy - lambdak - param.mu; 
			eigenvals[1] = -2*t*cy + lambdak - param.mu; 

			ComplexType v2(conj(hABk)/lambdak);
			if (lambdak!=0.0) {
				eigenvects(0,1) = 1.0; 
				eigenvects(1,1) = v2; 
				eigenvects(0,0) = 1.0; 
				eigenvects(1,0) = -v2; 
				// Normalize
				FieldType normal(sqrt(1.+norm(v2)));
				for (size_t b=0; b<nbands; b++) for (size_t iorb=0; iorb<nbands; iorb++) eigenvects(iorb,b) /= normal;
			} else { // Hamiltonian matrix is diagonal and diagonal elements -2*t*cosky are degenerate
				eigenvects(0,1) = 1.0/sqrt(2.0); 
				eigenvects(1,1) = 1.0/sqrt(2.0); 
				eigenvects(0,0) = 1.0/sqrt(2.0); 
				eigenvects(1,0) = -1.0/sqrt(2.0);
			}
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
