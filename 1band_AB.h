// Model file for 1-band model with AB-sublattices --> 2 bands total
#ifndef SINGLEBAND_AB_H
#define SINGLEBAND_AB_H


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
		ComplexMatrixType   spinMatrix;
		ComplexMatrixType   chargeMatrix;
		std::vector<int>    spinOfEll;
		std::vector<size_t> orbOfEll;


		model(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
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
				std::cerr<<"Number of orbitals should be 2 for this model! Exiting ...\n";

			msize = nbands * nbands;

			// Note that here the "orbitals" denote a combined orbital and spin index
			// The basis is (A, up; B, up; A, down; B, down)
			orbOfEll[0] = 0; // A-sublattice
			orbOfEll[1] = 1; // B-sublattice

			setupInteractionMatrix();
		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		  FieldType t,tp,cx,cy,cxy;

			t = 1.0; tp = -0.0;
			cx = cos(k[0]); cy = cos(k[1]); 
      cxy = cos(k[0])*cos(k[1]);

			FieldType ekAB = -2*t*(cx+cy);
      FieldType ekAA = -4*tp*cxy;
			
			// Write Hamiltonian into eigenvects
			// Basis is (A, up; B, up; A, down; B, down)
      // Note that the Hamiltonian is block-diagonal in the spin
			
      for (size_t i=0; i<nbands; i++)
        for (size_t j=0; j<nbands; j++)
          eigenvects(i,j) = 0.0;
          
			eigenvects(0,0) = -param.mu + ekAA; 
			eigenvects(1,1) = -param.mu + ekAA; 

			eigenvects(0,1) = ekAB; 
			eigenvects(1,0) = ekAB; 
    
			eigen(eigenvals, eigenvects);

		}

		void setupInteractionMatrix() {
			FieldType U(param.U);
			int o1,o2,o3,o4, l1, l2, l3, l4;

			for (size_t ind1=0; ind1 < msize; ind1++) {
				for (size_t ind2=0; ind2<msize; ind2++) {
					spinMatrix(ind1,ind2) = 0.0;
					l1 = param.indexToOrb(ind1,0); l2 = param.indexToOrb(ind1,1);
					l3 = param.indexToOrb(ind2,0); l4 = param.indexToOrb(ind2,1);
					o1 = orbOfEll[l1] ; o2 = orbOfEll[l2];
					o3 = orbOfEll[l3] ; o4 = orbOfEll[l4];

					// U-terms
					if (o1 == o2 && o1 == o3 && o1 == o4) {
						spinMatrix(ind1,ind2) = U;
					}
				}
			} 
		
			}


		std::complex<Field> calcSus(const ComplexMatrixType& sus, const std::string& component = "zz") const {
			std::complex<Field> chiPhys(0.0);
			int o1,o2,o3,o4, l1, l2, l3, l4;
			for (size_t ind1=0; ind1 < msize; ind1++) {
				for (size_t ind2=0; ind2<msize; ind2++) {
					l1 = param.indexToOrb(ind1,0); l2 = param.indexToOrb(ind1,1);
					l3 = param.indexToOrb(ind2,0); l4 = param.indexToOrb(ind2,1);
					o1 = orbOfEll[l1] ; o2 = orbOfEll[l2];
					o3 = orbOfEll[l3] ; o4 = orbOfEll[l4];

					if (o1 == o2 && o3 == o4) {
            chiPhys += 0.25 * sus(ind1,ind2);
					}
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
