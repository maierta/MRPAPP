// Model file for Sr2RuO4 with spin-orbit coupling --> 6 bands total
#ifndef SRRUO_SO_H
#define SRRUO_SO_H


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
		std::vector<int> spinOfEll;
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
			if (nbands != 6) 
				std::cerr<<"Number of orbitals should be 6! Exiting ...\n";

			msize = 36;

			// Note that here the "orbitals" denote a combined orbital and spin index
			// The basis is (xz,up;yz,up;xy,down ; xz,down;yz,down;xy,up)
			spinOfEll[0] = +1; orbOfEll[0] = 0;
			spinOfEll[1] = +1; orbOfEll[1] = 1; 
			spinOfEll[2] = -1; orbOfEll[2] = 2; 
			spinOfEll[3] = -1; orbOfEll[3] = 0; 
			spinOfEll[4] = -1; orbOfEll[4] = 1; 
			spinOfEll[5] = +1; orbOfEll[5] = 2; 

			setupInteractionMatrix();

		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		    	FieldType t1,t2,t3,t4,t5,tp,cx,cy,cxy,c2x,c2y;

			t1 = 0.088; t2 = 0.009; t3 = 0.080; t4 = 0.040; t5 = 0.005; 
			tp = 0.0;
			// tp = 0.0044;

			cx = cos(k[0]); cy = cos(k[1]); cxy = cos(k[0])*cos(k[1]);
			c2x = cos(2*k[0]); c2y = cos(2*k[1]);

			FieldType ekXY    = -2*t3*(cx+cy)-4*t4*cxy-2*t5*(c2x+c2y) - param.mu;
			// FieldType ekXY    = -2*t3*(cx+cy)-4*t4*cxy-2*t5*(c2x+c2y) - param.mu - 0.025;
			FieldType ekXZ    = -2*t1*cx-2*t2*cy - param.mu;
			FieldType ekYZ    = -2*t2*cx-2*t1*cy - param.mu;
			FieldType gk      = -4*tp*sin(k[0]) * sin(k[1]);

			FieldType lso     = 0.5*param.lambda_SO;
			
			// Write Hamiltonian into eigenvects
			// Basis is (xz,up;yz,up;xy,down ; xz,down;yz,down;xy,up)
			// Note that Hamiltonian is block-diagonal in 2 pseudospin blocks
			
			const ComplexType ii = ComplexType(0.0,1.0);
			ComplexMatrixType temp(3,3);
			VectorType evals(3);

			for (size_t i=0; i<nbands; i++) 
				for (size_t j=0; j<nbands; j++) eigenvects(i,j) = ComplexType(0.,0.);
					
			for (size_t i=0; i<3; i++) 
				for (size_t j=0; j<3; j++) temp(i,j) = ComplexType(0.,0.);


			temp(0,0) = ekXZ;  
			temp(1,0) = gk + ii*lso; 
			temp(0,1) = gk - ii*lso; 
			temp(0,2) =  ii*lso; 
			temp(2,0) = -ii*lso; 
			temp(1,1) = ekYZ; 
			temp(2,2) = ekXY; 
			temp(1,2) = -lso; 
			temp(2,1) = -lso; 

			eigen(evals,temp);

			for (size_t b=0; b<3; b++) {
					eigenvals[b] = evals[b]; // pseudospin up
					eigenvals[b+3] = evals[b]; // pseudospin down
				for (size_t l=0; l<3; l++) {
					eigenvects(l,b) = temp(l,b); // pseudospin up
				}
				eigenvects(3,b+3) = -conj(eigenvects(0,b)); // pseudospin down xz
				eigenvects(4,b+3) = -conj(eigenvects(1,b)); // pseudospin down yz
				eigenvects(5,b+3) =  conj(eigenvects(2,b)); // pseudospin down xy
			}
			// eigenvects(3,3) = ekXZ;  
			// eigenvects(4,3) = gk - ii*lso; 
			// eigenvects(3,4) = gk + ii*lso; 
			// eigenvects(3,5) =  ii*lso; 
			// eigenvects(5,3) = -ii*lso; 
			// eigenvects(4,4) = ekYZ; 
			// eigenvects(5,5) = ekXY; 
			// eigenvects(4,5) = +lso; 
			// eigenvects(5,4) = +lso; 

			// std::cout << k[0] << "," << k[1] << " , " << eigenvects << "\n";
			// std::cout << k[0] << "," << k[1] << " , " << ekXZ << "," << ekYZ << "," << ekXY << "\n";

			// std::cout << "eigenvects: \n" << eigenvects << "\n"; 

		}

		void setupInteractionMatrix() {
			// size_t nOrb(6);
			FieldType U(param.U);
			FieldType Up(param.Up);
			FieldType J(param.J);
			FieldType Jp(param.Jp);
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
					// U'-terms
					if (o1 != o2 && o1 == o3 && o2 == o4) 
						if (s1 == s3 && s2 == s4)                spinMatrix(ind1,ind2) += Up;
					if (o1 == o2 && o1 != o3 && o3 == o4) 
						if (s1 == s2 && s3 == s4)                spinMatrix(ind1,ind2) -= Up;
					// J-terms
					if (o1 == o2 && o1 != o3 && o3 == o4) 
						if (s1 == s3 && s2 == s4)                spinMatrix(ind1,ind2) += J;
					if (o1 != o2 && o1 == o3 && o2 == o4) 
						if (s1 == s2 && s3 == s4)                spinMatrix(ind1,ind2) -= J;
					// J'-terms
					if (o1 != o2 && o1 == o4 && o2 == o3) {
						if (s1 == s3 && s2 == s4 && s1 != s2)    spinMatrix(ind1,ind2) += Jp;
						if (s1 == s2 && s3 == s4 && s1 != s3)    spinMatrix(ind1,ind2) -= Jp;
					}

				}
			} 
		
			

			}


		std::complex<Field> calcSus(const ComplexMatrixType& sus, const std::string& component = "zz") const {
			std::complex<Field> chiPhys(0.0);
			int s1,s2,s3,s4,o1,o2,o3,o4, l1, l2, l3, l4;
			for (size_t ind1=0; ind1 < msize; ind1++) {
				for (size_t ind2=0; ind2<msize; ind2++) {
					l1 = param.indexToOrb(ind1,0); l2 = param.indexToOrb(ind1,1);
					l3 = param.indexToOrb(ind2,0); l4 = param.indexToOrb(ind2,1);
					s1 = spinOfEll[l1]; s2 = spinOfEll[l2];
					s3 = spinOfEll[l3]; s4 = spinOfEll[l4];
					o1 = orbOfEll[l1] ; o2 = orbOfEll[l2];
					o3 = orbOfEll[l3] ; o4 = orbOfEll[l4];

					if (o1 == o2 && o3 == o4) {
						if (component == "zz" && s1 == s2 && s3 == s4) {
							chiPhys += sus(ind1,ind2) * ComplexType(s1,0) * ComplexType(s3,0);
						} else if (component == "+-" && s1 == s3 && s2 == s4 && s1 != s4) // for +- - susc.
							chiPhys += sus(ind1,ind2);
					}
				}
			}
			return 0.25*chiPhys;

		}

		std::complex<Field> calcSCGap(VectorType& k, size_t band, ComplexMatrixType& Uk) {
			return ComplexType(0.0,0.0); // uses function in gaps3D.h directly for now
		}

	};
}

#endif
