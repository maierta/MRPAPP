// Model file for Sr2RuO4 with spin-orbit coupling --> 6 bands total
#ifndef MODELFROMFILESO_H
#define MODELFROMFILESO_H


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
		VectorType dx,dy,dz,ht;
		std::vector<size_t> orb1,orb2;
		int nLines;

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
			readCSVFile();
			msize = int(nbands*nbands);

			// Note that here the "orbitals" denote a combined orbital and spin index
			// The basis is is assumed to be (orb 0, up, orb 0, down, orb 1, up, orb 1, down, ...)
			for (size_t il=0; il<nbands; il++) {
				orbOfEll[il]  = int(il/2);
				spinOfEll[il] = -2 * (il % 2) + 1;
			}

			setupInteractionMatrix();

		}

		void readCSVFile() {
			std::string file = param.tbfile;
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format dx,dy,dz,orb1,orb2,t
			size_t length = 6;
			size_t nLinesTotal(data.size()/length);
			if (conc.rank()==0) std::cout << "tb file contains " << nLinesTotal << " lines\n";
			for (size_t i = 0; i < nLinesTotal; i++) {
				size_t l1(size_t(data[i*length+3]-1));
				size_t l2(size_t(data[i*length+4]-1));
				if (l1<=l2) {
					dx.push_back  (data[i*length]);
					dy.push_back  (data[i*length+1]);
					dz.push_back  (data[i*length+2]);
					orb1.push_back(l1);
					if (orb1[orb1.size()] > nbands)
						std::cerr<<"Number of orbitals exceeds maximum! Exiting ...\n";
					orb2.push_back(l2);
					if (orb2[orb2.size()] > nbands)
						std::cerr<<"Number of orbitals exceeds maximum! Exiting ...\n";
					ht.push_back  (data[i*length+5]);
				}
			}
			nLines = dx.size();
			if (conc.rank()==0) std::cout << nLines <<" entries used from tb input file\n";
		}

		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {

			FieldType exponent(0.);
			int n = eigenvects.n_col();
			std::vector<FieldType> ks(k);
			for (int i=0;i<n;i++) for (int j=0;j<n;j++) eigenvects(i,j)=ComplexType(0.,0.);

			// Now add hopping terms
			for (int i = 0; i < nLines; ++i)
			{
				// if (orb1[i]<=orb2[i]) {
					exponent = (ks[0]*dx[i] + ks[1]*dy[i] + ks[2]*dz[i]);
					ComplexType cs(cos(exponent),sin(exponent));
					eigenvects(orb1[i],orb2[i]) += ht[i] * cs;
				// }
			}

			for (size_t i=0;i<nbands;i++) eigenvects(i,i) -= param.mu;

			// Now add SO coupling terms
			// Note: For this we need to know the basis
			// FieldType lso     = 0.5*param.lambda_SO;
			// const ComplexType ii = ComplexType(0.0,1.0);


			eigen(eigenvals,eigenvects);
	
			return;

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
