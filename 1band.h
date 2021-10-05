// Model file for  single band model
#ifndef SINGLEBAND_H
#define SINGLEBAND_H

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
		const ComplexType ii;

		ComplexMatrixType spinMatrix;
		ComplexMatrixType chargeMatrix;


		model(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			ii(ComplexType(0.0,1.0)),
			nbands(param.nOrb),
			spinMatrix(nbands*nbands,nbands*nbands),
			chargeMatrix(nbands*nbands,nbands*nbands)
		{
			std::cout << "Single-band Model \n";
			if (nbands != 1) 
				std::cerr<<"Number of orbitals should be 1! Exiting ...\n";

			msize = 1;
			setupInteractionMatrix();

		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
	    	FieldType t, tp, cy, cx;
			t = param.hopping_t;
			tp = param.hopping_tp;
			cx = cos(k[0]); cy = cos(k[1]);
			eigenvals[0] = -2*t*(cx+cy)-4*tp*cx*cy - param.mu;
			eigenvects(0,0) = 1.0;
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

		std::complex<Field> calcSCGap(VectorType& k, size_t band, ComplexMatrixType& Uk) { // dummy function; handled by gaps3g.h directly
			return ComplexType(0.,0.); 
		}

	};
}

#endif
