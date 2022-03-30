#ifndef TBFROMFILE_H
#define TBFROMFILE_H

// Model file for tight-binding model with hopping parameters read from file

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
		typedef std::complex<Field>		ComplexType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 				FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		VectorType dx,dy,dz,ht;
		std::vector<size_t> orb1,orb2;
		int nLines;
		ComplexMatrixType Lm;
		
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
			std::cout << "In model \n";
			if (param.LS==1) setupLMatrix();
			std::cout << "Reading tb file \n";
			readCSVFile();
			// setupInteractionMatrix();
			// if (param.sublattice==1) fixdr();
			std::cout << "Setting up interaction matrix\n";
			setupInteractionMatrix2();
			std::cout << "Interaction matrix set up\n";
		}
		
		void readCSVFile() {
			std::string file = param.tbfile;
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format dx,dy,dz,orb1,orb2,t
			size_t length = 6;
			// if (param.complexHopping) length=7;
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
					orb2.push_back(l2);
					ht.push_back  (data[i*length+5]);
					// if (param.complexHopping) {
						// hti.push_back  (data[i*length+6]);
					// } else
					// hti.push_back  (0.0);
				}
			}
			nLines = dx.size();
			if (conc.rank()==0) std::cout << nLines <<" entries used from tb input file\n";
		}

		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects, int spin=1) {

			FieldType exponent(0.);
			int n = eigenvects.n_col();
			std::vector<FieldType> ks(k);
			if(param.kTrafo==1) transformK(k,ks);
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

			// Then add approximate spin-orbit coupling term L*S (see Kreisel notes in email from Apr. 4 2013)
			// if (false) {
			//	ComplexType strength(0.0,0.5*param.hybStrength*float(spin));
			//	ComplexType strengthC(conj(strength));
			// //	// std::cout << "Adding hyb term " << strength << "\n";
			//	FieldType kd(ks[0]*param.deltax+ks[1]*param.deltay+ks[2]*param.deltaz);
			//	ComplexType pf(cos(-kd),sin(-kd));

			//	eigenvects(3,4) += 2*strength;
			//	eigenvects(8,9) += 2*strength;
			//	eigenvects(3,9) += 2.*strength*pf;
			//	eigenvects(4,8) += 2.*strengthC*pf;

			//	eigenvects(1,2) += 1*strength;
			//	eigenvects(6,7) += 1*strength;
			//	eigenvects(1,7) += 1.*strength*pf;
			//	eigenvects(2,6) += 1.*strengthC*pf;
			// }

			if (param.LS) {
				// Add full L*S term
				// First expand eigenvects matrix to twice the size
				ComplexMatrixType temp (2*nbands,2*nbands);
				// for (size_t i=0; i<nbands;i++) for (size_t j=0; j<nbands;j++) temp(i,j)=ComplexType(0.0,0.0);
				VectorType tempE(2*nbands,0.0);
				for (size_t i=0; i<nbands;i++) {
					for (size_t j=0; j<nbands;j++) {
						temp(i,j) = eigenvects(i,j);
						temp(i+nbands,j+nbands) = eigenvects(i,j);
					}
				}
				// for (size_t i=0; i<nbands;i++) temp(i,i) -= 0.00001; // small shift for ordering

				// Now add Lm matrix
				for (size_t i=0; i<2*nbands;i++) {
					for (size_t j=0; j<2*nbands;j++) {
						temp(i,j) += Lm(i,j);
					}
				}
				// eigen(eigenvals,eigenvects);
				eigen(tempE,temp);

				// std::cout << tempE << "\n";
				// for (size_t i=0;i<2*nbands;i++) std::cout << temp(i,2);
				// std::cout << "\n";
				// for (size_t i=0;i<2*nbands;i++) std::cout << temp(i,3);
				// std::cout << "\n";
				size_t iband(0);
				for (size_t i=0;i<nbands;i++) eigenvals[i] = tempE[2*i];
				for (size_t band=0; band<nbands;band++) {
					bool b=inspectEigenvector(temp,2*band);
					// if (band==1) std::cout << "b=" << b << "\n";
					// bool b(1);
					// FieldType r1(abs(temp(3,2*band))+abs(temp(3+5,2*band)));
					// FieldType r2(abs(temp(3+nbands,2*band))+abs(temp(3+5+nbands,2*band)));
					(b)?iband=2*band:iband=2*band+1;
					// iband=2*band;
					for (size_t j=0; j<nbands;j++) eigenvects(j,band) = temp(j,iband);
					normalize(eigenvects,band);

					// if (band==1) {
					//	std::cout << "Normalized EV: ";
					//	for (size_t i=0;i<nbands;i++) std::cout << eigenvects(i,band) << "\n";
					// }
				}
				return;
			}

			// if (param.hyb) ComplexMatrixType HkOrg(eigenvects);

			if (param.hyb) {
				// FieldType fofk(0.5-0.25*(cos(ks[0])+cos(ks[1])));
				// FieldType fofk(1.0);
				// FieldType fofk((sin(ks[1])>0) - (sin(ks[1])<0));
				FieldType fofk((ks[1]>0) - (ks[1]<0));

				ComplexType hyb(0.0,fofk*param.hybStrength);
				// ComplexType hyb(fofk*param.hybStrength,0.0);
				FieldType kd(ks[0]*param.deltax+ks[1]*param.deltay+ks[2]*param.deltaz);
				ComplexType pf(cos(-kd),sin(-kd));

				eigenvects(1,2) += hyb;
				eigenvects(6,7) += hyb;
				eigenvects(1,7) += hyb*pf;
				eigenvects(2,6) += hyb*pf;

				// eigenvects(1,6) += hyb*pf;
				// eigenvects(2,7) += hyb*pf;

				eigen(eigenvals,eigenvects);
				return;
			}

			eigen(eigenvals,eigenvects);
			phaseFactor(k,eigenvects);
			return;

		}

		inline void phaseFactor(const VectorType& kOrg, ComplexMatrixType& eigenvects) {
			// Note: Here we assume that in the TB input data, all orbitals sit on the same site
			if (param.sublattice==1) {
				std::vector<FieldType> ks(kOrg);
				if (param.kTrafo==1) transformK(kOrg,ks);
				FieldType exponent(-(ks[0]*param.deltax + ks[1]*param.deltay + ks[2]*param.deltaz));
				ComplexType cs(cos(exponent),sin(exponent));
				// ComplexType cs(cos(param.pi_f/3.),sin(param.pi_f/3.));
				// ComplexType cs(1.0,0.0);
				for (size_t orb=size_t(nbands/2); orb<nbands; orb++) for (size_t band=0; band<nbands; band++) {
					// eigenvects(orb,band) *= ComplexType(cos(exponent),sin(exponent));
					eigenvects(orb,band) *= cs;
				}
			}
		}

		void setupLMatrix() {
			// Assume the basis is (z^2,yz,xz,xy,x^2-y^2)
			// Setup the full 20 x 20 matrix, but only the upper triangle
			// const FieldType& lambda(param.hybStrength);
			if (conc.rank()==0) std::cout << "Setting up L matrix with lambda = " << param.hybStrength << "\n";
			for (size_t i=0;i<Lm.n_row();i++) for (size_t j=0;j<Lm.n_col();j++) Lm(i,j) = ComplexType(0.0,0.0);
			// 0:z^2; 1:yz; 2:xz; 3:xy; 4:x^2-y^2
			const ComplexType I (0.0,1.0);
			ComplexMatrixType Lplus(10,10);
			ComplexMatrixType Lz(10,10);

			// Lz for up
			Lz(1,2) =  I  ; // (yz,xz)
			Lz(3,4) =  2*I; // (xy,x2-y2)
			// L+
			Lplus(0,1) = I*sqrt(3.);
			Lplus(0,2) =  -sqrt(3.);
			Lplus(1,0) = -I*sqrt(3.);
			Lplus(1,3) = -1.0         ;
			Lplus(1,4) = -I       ;
			Lplus(2,0) =  sqrt(3.);
			Lplus(2,3) =  I       ;
			Lplus(2,4) = -1.0        ;
			Lplus(3,1) =  1.0      ;
			Lplus(3,2) =  -I        ;
			Lplus(4,1) =  I        ;
			Lplus(4,2) =  1.0      ;

			// Now extend to second Fe
			for (size_t i =0; i<5; i++) {
				for (size_t j =0; j<5; j++) {
					Lz(i+5,j+5) = Lz(i,j);
					Lplus(i+5,j+5) = Lplus(i,j);
				}
			}

			// Now build full Lmatrix Lm(0:20,0:20)
			for (size_t i =0; i<10; i++) {
				for (size_t j =0; j<10; j++) {
					Lm(i,j) =  Lz(i,j);       // spin up
					Lm(i+10,j+10) = -Lz(i,j); // spin down

					Lm(i,j+10) = Lplus(i,j); // L+
				}
			}

			// Finally muliply with 0.5*hybStrength
			for (size_t i =0; i<20; i++) for (size_t j =0; j<20; j++) Lm(i,j) *= 0.5*param.hybStrength;

			if (conc.rank()==0) std::cout << "Done setting up L matrix \n";

		}

		void normalize(ComplexMatrixType& matrix, size_t colIndex) {
			std::vector<ComplexType> vec(matrix.n_row(),0.0);
			matrix.getCol(colIndex,vec);
			FieldType r1(0.0);
			for (size_t i=0;i<vec.size();i++) {
				r1 += pow(abs(vec[i]),2);
			}
			for (size_t i=0;i<vec.size();i++) matrix(i,colIndex) /= sqrt(r1);
			// std::cout << "In Normalize: r1=" << r1 << "\n";
			return;
		}

		bool inspectEigenvector(ComplexMatrixType& matrix, size_t colIndex) {
			std::vector<ComplexType> vec(matrix.n_row(),0.0);
			matrix.getCol(colIndex,vec);
			FieldType r1(0.0);
			FieldType r2(0.0);
			for (size_t i=0;i<vec.size()/2;i++) {
				r1 += pow(abs(vec[i]),2);
				r2 += pow(abs(vec[i+vec.size()/2]),2);
			}
			bool r1gr2(0);
			(r1>r2)?r1gr2=1:r1gr2=0;
			return r1gr2;
		}


		// void setupInteractionMatrix() {
		// 	size_t nOrb(param.nOrb);
		// 	size_t msize(size_t(nOrb*nOrb));
		// 	size_t limit(nOrb);
		// 	FieldType U(param.U);
		// 	FieldType Up(param.Up);
		// 	FieldType J(param.J);
		// 	FieldType Jp(param.Jp);

		// 	if (param.sublattice==1) limit=nOrb<10?nOrb/2:5; //Note: This only works for Fe-type problems with 2 Fe in unit cell where each Fe has 5 orbitals

			
		// 	ComplexMatrixType spinSubMatrix(limit*limit,limit*limit);
		// 	ComplexMatrixType chargeSubMatrix(limit*limit,limit*limit);
				
		// 		// First the diagonal terms
		// 		for (size_t l1 = 0; l1 < limit; ++l1) {
		// 				for (size_t l2 = 0; l2 < limit; ++l2) {
		// 					size_t ind1 = l2+l1*limit;
		// 					if (l1==l2) {
		// 						spinSubMatrix(ind1,ind1)   = U+param.deltaU[l1];
		// 						chargeSubMatrix(ind1,ind1) = -U-param.deltaU[l1];;
		// 						} else {
		// 						spinSubMatrix(ind1,ind1)   = Up;
		// 						chargeSubMatrix(ind1,ind1) = Up-2*J;
		// 						}
		// 				}
		// 			}	
		// 		// Off-diagonal terms
		// 		for (size_t l1=0; l1 < limit; l1++) {
		// 			size_t ind1 = l1+l1*limit;
		// 			for (size_t l2=0; l2 < limit; l2++) {
		// 				size_t ind2 = l2+l2*limit;
		// 				if (l2!=l1) {
		// 					spinSubMatrix(ind1,ind2)   = J;
		// 					chargeSubMatrix(ind1,ind2) = -2.*Up+J;
		// 				}
		// 			}
		// 		}
		// 		// Finally the pair hopping terms
		// 		for (size_t l1=0; l1 < limit; l1++) {
		// 			for (size_t l2=0; l2 < limit; l2++) {
		// 				size_t ind1 = l2+l1*limit;
		// 				size_t ind2 = l1+l2*limit;
		// 				if (l2!=l1) {
		// 					spinSubMatrix(ind1,ind2)   = Jp;
		// 					chargeSubMatrix(ind1,ind2) = -Jp;
		// 				}
		// 			}
		// 		}
		// 	if (param.sublattice==0) {
		// 		for (size_t i=0; i<msize; i++) for (size_t j=0; j<msize; j++) {
		// 			spinMatrix(i,j) = spinSubMatrix(i,j);
		// 			chargeMatrix(i,j) = chargeSubMatrix(i,j);
		// 		}
		// 		} else {
		// 			for(size_t l1=0; l1<limit; l1++) for (size_t l2=0; l2<limit; l2++) {
		// 			for(size_t l3=0; l3<limit; l3++) for (size_t l4=0; l4<limit; l4++) {
						
		// 				size_t ind1=l2+l1*limit;
		// 				size_t ind2=l4+l3*limit;
	
		// 				size_t ind3=l2+l1*nOrb;
		// 				size_t ind4=l4+l3*nOrb;

		// 				spinMatrix(ind3,ind4) = spinSubMatrix(ind1,ind2);
		// 				chargeMatrix(ind3,ind4) = chargeSubMatrix(ind1,ind2);

		// 				ind3=l2+limit+(l1+limit)*nOrb; // position of 2. Fe d-orbitals is shifted by limit wrt 1. Fe d-orbs.
		// 				ind4=l4+limit+(l3+limit)*nOrb; // position of 2. Fe d-orbitals is shifted by limit wrt 1. Fe d-orbs.

		// 				spinMatrix(ind3,ind4) = spinSubMatrix(ind1,ind2);
		// 				chargeMatrix(ind3,ind4) = chargeSubMatrix(ind1,ind2);
		// 			}
		// 			}	
		// 		}					
		// 	}

		void setupInteractionMatrix2() {
			size_t nOrb(param.nOrb);
			FieldType U(param.U);
			FieldType Up(param.Up);
			FieldType J(param.J);
			FieldType Jp(param.Jp);
			// if (conc.rank()==0) std::cout << "U:"<<U<<" U'"<<Up<<" J:"<<J<<" Jp:"<<Jp<<"\n";

			for (size_t i=0; i<spinMatrix.n_row(); ++i) for (size_t j=0; j<spinMatrix.n_col(); ++j) {
				spinMatrix(i,j) = 0.0;
				chargeMatrix(i,j) = 0.0;
			}

			// First the diagonal terms (U and U')
			for (size_t l1 = 0; l1 < nOrb; ++l1) {
				for (size_t l2 = 0; l2 < nOrb; ++l2) {
					if (param.orbToSite[l1] != param.orbToSite[l2]) continue; // orbital l1 and l2 belong to different sites
					// size_t ind1 = l2+l1*nOrb;
					size_t ind1 = param.OrbsToIndex(l1,l2);
					if (l1==l2) {
						spinMatrix(ind1,ind1)   = U;
						chargeMatrix(ind1,ind1) = -U;
						} else {
						spinMatrix(ind1,ind1)   = Up;
						chargeMatrix(ind1,ind1) = Up-2*J;
						}
					}
				}	
			// The the off-diagonal terms
			for (size_t l1=0; l1 < nOrb; l1++) {
				for (size_t l2=0; l2 < nOrb; l2++) {
					if (param.orbToSite[l1] != param.orbToSite[l2]) continue; // orbital l1 and l2 belong to different sites
					if (l2!=l1) {
					// size_t ind1 = l1+l1*nOrb;
					size_t ind1 = param.OrbsToIndex(l1,l1);
					// size_t ind2 = l2+l2*nOrb;
					size_t ind2 = param.OrbsToIndex(l2,l2);
						spinMatrix(ind1,ind2)   = J;
						chargeMatrix(ind1,ind2) = -2.*Up+J;
					}
				}
			}
			// Finally the pair hopping terms
			for (size_t l1=0; l1 < nOrb; l1++) {
				for (size_t l2=0; l2 < nOrb; l2++) {
					if (param.orbToSite[l1] != param.orbToSite[l2]) continue; // orbital l1 and l2 belong to different sites
					if (l2!=l1) {
					// size_t ind1 = l2+l1*nOrb;
					size_t ind1 = param.OrbsToIndex(l1,l2);
					// size_t ind2 = l1+l2*nOrb;
					size_t ind2 = param.OrbsToIndex(l2,l1);
						spinMatrix(ind1,ind2)   = Jp;
						chargeMatrix(ind1,ind2) = -Jp;
					}
				}
			}

			// std::cout << "U matrix: " << spinMatrix << "\n";
			}


		void fixdr() { // This is to shift position of all orbitals to the same site
			for (int i = 0; i < nLines; ++i)
			{
				if (orb1[i]<nbands/2 && orb2[i]>=nbands/2) {
					dx[i]=dx[i]+param.deltax;
					dy[i]=dy[i]+param.deltay;
					dy[i]=dy[i]+param.deltaz;
				}
				else if (orb1[i]>=nbands/2 && orb2[i]<nbands/2) {
					dx[i]=dx[i]-param.deltax;
					dy[i]=dy[i]-param.deltay;
					dy[i]=dy[i]-param.deltaz;
				}
			}
		}

		void transformK(const std::vector<FieldType>& kIn, std::vector<FieldType>& kOut) {
			kOut[0] = kIn[0]*param.WTrafo(0,0)+kIn[1]*param.WTrafo(1,0)+kIn[2]*param.WTrafo(2,0);
			kOut[1] = kIn[0]*param.WTrafo(0,1)+kIn[1]*param.WTrafo(1,1)+kIn[2]*param.WTrafo(2,1);
			kOut[2] = kIn[0]*param.WTrafo(0,2)+kIn[1]*param.WTrafo(1,2)+kIn[2]*param.WTrafo(2,2);
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
			Field factor(1.0);
			// if (param.sublattice==1) factor=param.nSitesPerUnitCell;
			return chiPhys/factor;

		}

		std::complex<Field> calcSCGap(VectorType& k, size_t band, ComplexMatrixType& Uk) { // dummy function; handled by gaps3g.h directly
			return ComplexType(0.,0.); 
		}
	};
}

#endif
