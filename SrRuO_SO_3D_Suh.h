// Model file for Sr2RuO4 with spin-orbit coupling --> 6 bands total; Suh bands
#ifndef SRRUO_SO_3D_SUH_H
#define SRRUO_SO_3D_SUH_H


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
			if (dim != 3) 
				std::cerr<<"Dimension should be 3! Exiting ...\n";

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
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& H0) {

			FieldType ekXZ, ekYZ, ekXY;
			FieldType gxzyz, Txzxy, Tyzxy;
			FieldType cx,cy,cz,cxy,c2x,c2y,sx,sy,cz2,cx2,cy2,sx2,sy2,sz2,s2x,s2y;

			sx = sin(k[0]); sy = sin(k[1]); sx2 = sin(0.5*k[0]); sy2 = sin(0.5*k[1]); sz2 = sin(0.5*k[2]);
			s2x = sin(2*k[0]); s2y = sin(2*k[1]);
			cx = cos(k[0]); cy = cos(k[1]); cz = cos(k[2]); cxy = cos(k[0])*cos(k[1]);
			c2x = cos(2*k[0]); c2y = cos(2*k[1]); 
			cx2 = cos(0.5*k[0]); cy2 = cos(0.5*k[1]); cz2 = cos(0.5*k[2]);


			FieldType t1,t2,t3,t4,t5,tp,tpp,tppp,tint;
			FieldType t6, t7, t8, t9, mu, muxy;
			FieldType t1xy, t2xy, t3xy, t4xy, t5xy, t6xy;

			t1 = 0.3624; t2 = 0.134; 
			t3 = 0.001021; t4 = 0.005727; t5 = 0.04401; 
			t6 = 0.01393; t7 = 0.00752; t8 = -0.002522; t9 = 0.0000228; 
			// t3 = 0; t4 = 0; t5 = 0; 
			// t6 = 0; t7 = 0; t8 = 0; t9 = 0; 

			mu = 0.4385;

			t1xy = 0.262; t2xy = -0.03423; t3xy = 0.04373; 
			t4xy = -0.008069; t5xy = 0.003159; t6xy = -0.001811; 
			// t4xy = 0; t5xy = 0; t6xy = 0; 

			muxy = 0.2186;
			tp = 0.01625; tpp = -0.01995; tppp = 0.00394; tint = -2*0.008304;
			// tp = 0; tpp = 0; tppp = 0; tint = 0;

			ekXZ    = -2*t1*cx-2*t2*cy-2*t3*c2x-2*t4*c2y-4*t5*cx*cy-4*t6*c2x*cy-4*t7*cx*c2y-2*t8*(cz-1)-8*t9*cx2*cy2*cz2 - mu - param.mu;
			ekYZ    = -2*t1*cy-2*t2*cx-2*t3*c2y-2*t4*c2x-4*t5*cy*cx-4*t6*c2y*cx-4*t7*cy*c2x-2*t8*(cz-1)-8*t9*cy2*cx2*cz2 - mu - param.mu;
			ekXY    = -2*t1xy*(cx+cy)-2*t2xy*(c2x+c2y)-4*t3xy*cx*cy -4*t4xy*(c2x*cy+cx*c2y) - 2*t5xy*(cz-1) - 8*t6xy*cx2*cy2*cz2 - muxy - param.mu;

			gxzyz   = -4*tp*sx*sy - 4*tpp*sx2*sy2*cz2 -4*tppp*(s2x*sy+sx*s2y);
			Txzxy   = -4*tint*cx2*sy2*sz2; 
			Tyzxy   = -4*tint*sx2*cy2*sz2; 


			FieldType lso     = 0.5*param.lambda_SO;

			// Write Hamiltonian into eigenvects
			// Basis is (xz,up;yz,up;xy,down ; xz,down;yz,down;xy,up)
			
			const ComplexType ii = ComplexType(0.0,1.0);
			// ComplexMatrixType H0(6,6);

			// for (size_t i=0; i<nbands; i++) 
				// for (size_t j=0; j<nbands; j++) H0(i,j) = ComplexType(0.,0.);
					
			H0(0,0) = ekXZ;
			H0(0,1) = gxzyz - ii*lso;
			H0(0,2) = ii*lso;
			H0(0,3) = 0;
			H0(0,4) = 0;
			H0(0,5) = Txzxy;

			H0(1,0) = gxzyz + ii*lso;
			H0(1,1) = ekYZ;
			H0(1,2) = -lso;
			H0(1,3) = 0;
			H0(1,4) = 0;
			H0(1,5) = Tyzxy;

			H0(2,0) = -ii*lso;
			H0(2,1) = -lso;
			H0(2,2) = ekXY;
			H0(2,3) = Txzxy;
			H0(2,4) = Tyzxy;
			H0(2,5) = 0;

			H0(3,0) = 0;
			H0(3,1) = 0;
			H0(3,2) = Txzxy;
			H0(3,3) = ekXZ;
			H0(3,4) = gxzyz + ii*lso;
			H0(3,5) = ii*lso;

			H0(4,0) = 0;
			H0(4,1) = 0;
			H0(4,2) = Tyzxy;
			H0(4,3) = gxzyz - ii*lso;
			H0(4,4) = ekYZ;
			H0(4,5) = lso;

			H0(5,0) = Txzxy;
			H0(5,1) = Tyzxy;
			H0(5,2) = 0;
			H0(5,3) = -ii*lso;
			H0(5,4) = lso;
			H0(5,5) = ekXY;


			// Optionally add k-SOC terms
			if (param.k_SOC) {
				FieldType t12z, t56z, tdxy, td;

				if (param.Case == "Suh4") {
						t12z = 0; t56z = 0.015; tdxy = 0; td = 0; 			
				} else {
					std::cerr << "Case not implemented \n";
					exit(0);
				}


				FieldType etax   =  8*t12z*cx2*sy2*sz2;
				FieldType etay   = -8*t12z*sx2*cy2*sz2;
				FieldType gammax =  8*t56z*cx2*sy2*sz2;
				FieldType gammay = -8*t56z*sx2*cy2*sz2;
				FieldType alpha  =  4*tdxy*sx*sy;
				FieldType beta   =  2*td*(cx-cy);

				H0(0,2) += alpha - ii*beta;
				H0(0,4) += etax - ii*etay;
				H0(0,5) += -ii*gammay;

				H0(1,2) += -beta - ii*alpha;
				H0(1,3) += -etax + ii*etay;
				H0(1,5) += -ii*gammax;

				H0(2,0) += alpha + ii*beta;
				H0(2,1) += -beta + ii*alpha;
				H0(2,3) += -ii*gammay;
				H0(2,4) += -ii*gammax;

				H0(3,1) += -etax - ii*etay;
				H0(3,2) += ii*gammay;
				H0(3,5) += -alpha - ii*beta;

				H0(4,0) += etax + ii*etay;
				H0(4,2) += ii*gammax;
				H0(4,5) += beta - ii*alpha;

				H0(5,0) += ii*gammay;
				H0(5,1) += ii*gammax;
				H0(5,3) += -alpha + ii*beta;
				H0(5,4) += beta + ii*alpha;
			}

			eigen(eigenvals,H0);

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
						} else if (component == "+-" && s1 == s3 && s2 == s4 && s1 != s4) { // for +- - susc.
							chiPhys += sus(ind1,ind2);
						} else if (component == "xx" && s1 != s2 && s3 != s4) { // for xx - susc.
							chiPhys += sus(ind1,ind2);
						} else if (component == "yy" && s1 != s2 && s3 != s4) { // for yy - susc.
							chiPhys -= sus(ind1,ind2) * ComplexType(s1*s4,0);
						}
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
