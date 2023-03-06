// Model file for 1-band model with explicit spin --> 2 bands total
#ifndef SINGLEBAND_WSPIN_H
#define SINGLEBAND_WSPIN_H


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
	class SingleBand_wSpin {
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
		ComplexMatrixType   spinMatrix;
		ComplexMatrixType   chargeMatrix;
		std::vector<int>    spinOfEll;
		std::vector<size_t> orbOfEll;


		SingleBand_wSpin(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb),
			spinMatrix(nbands*nbands,nbands*nbands),
			chargeMatrix(nbands*nbands,nbands*nbands),
			spinOfEll(nbands),
			orbOfEll(nbands)
		{
			if (nbands != 2) 
				std::cerr<<"Number of orbitals should be 2! Exiting ...\n";

			msize = 4;

			// Note that here the "orbitals" denote a combined orbital and spin index
			// The basis is (up;down)
			spinOfEll[0] = +1; orbOfEll[0] = 0;
			spinOfEll[1] = -1; orbOfEll[1] = 0; 

		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		    	FieldType t,tp,cx,cy,cxy;

			t = 1.0; tp = -0.3;
			cx = cos(k[0]); cy = cos(k[1]); cxy = cos(k[0])*cos(k[1]);

			FieldType ek    = -2*t*(cx+cy)-4*tp*cxy - param.mu;
			
			// Write Hamiltonian into eigenvects
			// Basis is (spin up, spin down)
			// Note that Hamiltonian is block-diagonal in 2 pseudospin blocks
			
			eigenvals[0] = ek; // pin up
			eigenvals[1] = ek; // spin down
			eigenvects(0,0) = 1.0; 
			eigenvects(0,1) = 0.0; 
			eigenvects(1,1) = 1.0;
			eigenvects(1,0) = 0.0;

		}

		void setupInteractionMatrix() {
			FieldType U(param.U);
			int s1,s2,s3,s4,o1,o2,o3,o4, l1, l2, l3, l4;

			for (size_t ind1=0; ind1 < msize; ind1++) {
				for (size_t ind2=0; ind2<msize; ind2++) {
					spinMatrix(ind1,ind2) = 0.0;
					l1 = param.indexToOrb(ind1,0); l2 = param.indexToOrb(ind1,1);
					l3 = param.indexToOrb(ind2,0); l4 = param.indexToOrb(ind2,1);
					s1 = spinOfEll[l1]; s2 = spinOfEll[l2];
					s3 = spinOfEll[l3]; s4 = spinOfEll[l4];
					o1 = orbOfEll[l1] ; o2 = orbOfEll[l2];
					o3 = orbOfEll[l3] ; o4 = orbOfEll[l4];

					// U-terms
					if (o1 == o2 && o1 == o3 && o1 == o4) {
						if (s1 == -s2 && s1 == s3 && s1 == -s4)  spinMatrix(ind1,ind2) += U;
						if (s1 == s2 && s1 == -s3 && s1 == -s4)  spinMatrix(ind1,ind2) -= U;
					}
				}
			} 
		
			

			}


		// std::complex<Field> calcSus(const ComplexMatrixType& sus, const std::string& component = "zz") const {
		// 	std::complex<Field> chiPhys(0.0);
		// 	int s1,s2,s3,s4,o1,o2,o3,o4, l1, l2, l3, l4;
		// 	for (size_t ind1=0; ind1 < msize; ind1++) {
		// 		for (size_t ind2=0; ind2<msize; ind2++) {
		// 			l1 = param.indexToOrb(ind1,0); l2 = param.indexToOrb(ind1,1);
		// 			l3 = param.indexToOrb(ind2,0); l4 = param.indexToOrb(ind2,1);
		// 			s1 = spinOfEll[l1]; s2 = spinOfEll[l2];
		// 			s3 = spinOfEll[l3]; s4 = spinOfEll[l4];
		// 			o1 = orbOfEll[l1] ; o2 = orbOfEll[l2];
		// 			o3 = orbOfEll[l3] ; o4 = orbOfEll[l4];

		// 			if (o1 == o2 && o3 == o4) { // always satisfied here for 1-band case!!
		// 				if (component == "zz" && s1 == s2 && s3 == s4) {
		// 					chiPhys += sus(ind1,ind2) * ComplexType(s1,0) * ComplexType(s3,0);
		// 				} else if (component == "+-" && s1 == s3 && s2 == s4 && s1 != s4) // for +- - susc.
		// 					chiPhys += sus(ind1,ind2);
		// 			}
		// 		}
		// 	}
		// 	return 0.25*chiPhys;

		// }

		std::complex<Field> calcSus(const ComplexMatrixType& sus, const std::string& component = "zz") const {
			std::complex<Field> chiPhys(0.0);

			// std::cout << "chiMatrix:" << sus << "\n";

			if (component == "zz") {
				chiPhys = sus(0,0) + sus(3,3) - sus(0,3) - sus(3,0);
			} else if (component == "+-") {
				chiPhys = sus(1,1) + sus(2,2);
			}

			return 0.25*chiPhys;

		}

	};
}

#endif
