#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H


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
#include "SrRuO.h"
#include "BaFeAs_5orb.h"
#include "KFe2Se2.h"
#include "FourOrbital.h"
#include "bilayer.h"

namespace rpa {
	extern "C" void
#ifdef glyph
	zheev_
#else
	zheev
#endif
	(char *,char *,int *,std::complex<double> *,int *,
						  double *,std::complex<double> *,
						  int *,double *,int *);
	inline void GEEV(char jobvl,char jobvr,int n,psimag::Matrix<std::complex<double> > &A, int lda,
					 std::vector<double>  &w,std::vector<std::complex<double> > &work,
					 int lwork,std::vector<double> &rwork,int *info)
					 {
#ifdef glyph
						zheev_
#else
						zheev
#endif
						(&jobvl,&jobvr,&n,&(A(0,0)),&lda,
						      &(w[0]),&(work[0]),&lwork,&(rwork[0]),info);
						}

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class bandstructure {

	private:
		typedef MatrixTemplate<Field>	MatrixType;
		typedef std::complex<Field>		ComplexType;
		typedef MatrixTemplate<ComplexType>	ComplexMatrixType;
		typedef std::vector<Field>      VectorType;
		typedef Field					FieldType;
		typedef momentumDomain<Field,MatrixTemplate,ConcurrencyType> kDomain;
		typedef PsimagLite::Range<ConcurrencyType> RangeType;

		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t nbands;
		int nLines;
		VectorType dx,dy,dz,ht;
		// VectorType hti;
		std::vector<size_t> orb1,orb2;
		kDomain kmesh_;
		const kDomain& kmesh;
		bool caching_;
		std::vector<bool> cachedK;
		std::vector<VectorType> ev;
		std::vector<ComplexMatrixType> ak;
		ComplexMatrixType Lm;

	public:

		bandstructure(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
				  ConcurrencyType& concurrency, const kDomain& kmeshIn,
				  bool caching):
			param(parameters),
			conc(concurrency),
			nbands(param.nOrb),
			kmesh_(param,conc,0,2), // not being used
			kmesh(kmeshIn),
			// caching_((kmesh.nktot<=2e6)?caching:false),
			caching_(false),
			cachedK(kmesh.nktot,false),
			ev(caching_?kmesh.nktot:0,VectorType(nbands)),
			ak(caching_?kmesh.nktot:0,ComplexMatrixType(nbands,nbands)),
			Lm(2*nbands,2*nbands)
		{
			// if (kmesh.nktot>=16384) caching_=false;
			if (param.tbfile!="") readCSVFile();
			if (conc.rank()==0) std::cout << "Caching=" << caching_ << "\n";
			if (param.LS==1) setupLMatrix();
			// if (param.sublattice==1) fixdr();
		}

		// Constructor without kmesh input
		bandstructure(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
				  ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			nbands(param.nOrb),
			kmesh_(param,conc,0,2),
			kmesh(kmesh_),
			caching_(false),
			cachedK(0),
			ev(0),
			ak(0,ComplexMatrixType(0,0)),
			Lm(2*nbands,2*nbands)
		{
			// if (kmesh.nktot>=16384) caching_=false;
			if (param.tbfile!="") readCSVFile();
			std::cout << "Caching=" << caching_ << "\n";
			if (param.LS==1) setupLMatrix();
			// if (param.sublattice==1) fixdr();
		}


		void getEkAndAk(VectorType& kvec, VectorType& eigenvals, ComplexMatrixType& eigenvects, int spin=1)  {
			if(!caching_) {
				getBands(kvec,eigenvals,eigenvects,spin);
				phaseFactor(kvec,eigenvects);
				return;
			} else {
				// NOTE: This scheme only works when tb input file has all atoms within a unit cell on the same site, i.e.
				// when the phase-factors associated with positions with a unit cell are taken into account here. If the tb input
				// file has these shifts already in it, the phase factors will already be taken care of and are included in the
				// eigenvectors. But in this case, since the caching maps the input k-vector to the 1. BZ, the phase factors
				// won't be correct. So for now, we will switch off mapping to 1. BZ
				// Also, when dim=2 and q has a finite qz, there seems to be a problem. Therefore caching is turned off!
				size_t ik; FieldType residue;
				VectorType kOrg(kvec);
				kmesh.kToik(kvec,ik,residue,false); // false = mapping to 1.BZ switched off
				if (fabs(residue)<=1.0e-5) { // kvec falls on k[ik] in kmesh
					if (cachedK[ik]) {
						eigenvals = ev[ik];
						eigenvects = ak[ik];
						phaseFactor(kOrg,eigenvects);
						return;
					} else {
					    VectorType kik(3); kmesh.momenta.getRow(ik,kik);
						getBands(kik,eigenvals,eigenvects,spin);
						ev[ik] = eigenvals;
						ak[ik] = eigenvects;
						cachedK[ik] = true;
						phaseFactor(kOrg,eigenvects);
						return;
					}
				} else {
					getBands(kOrg,eigenvals,eigenvects,spin);
					phaseFactor(kOrg,eigenvects);
					return;
				}
			}
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

		void transformK(const std::vector<FieldType>& kIn, std::vector<FieldType>& kOut) {
			kOut[0] = kIn[0]*param.WTrafo(0,0)+kIn[1]*param.WTrafo(1,0)+kIn[2]*param.WTrafo(2,0);
			kOut[1] = kIn[0]*param.WTrafo(0,1)+kIn[1]*param.WTrafo(1,1)+kIn[2]*param.WTrafo(2,1);
			kOut[2] = kIn[0]*param.WTrafo(0,2)+kIn[1]*param.WTrafo(1,2)+kIn[2]*param.WTrafo(2,2);
		}

		inline void getBands(const VectorType& k,
								   VectorType& eigenvals,
								   ComplexMatrixType& eigenvects,
								   int spin=1)  {

			// eigenvects=ComplexMatrixType(nbands,nbands,nbands,0.0);
#ifdef USE_SRRUO
			SrRuO<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvals,eigenvects);
#elif USE_BILAYER
			orthoIIBilayer<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			// bilayer<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvals,eigenvects);
#elif USE_BSCCObilayer
			BSCCObilayer<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvals,eigenvects);
#elif USE_BILAYER_FESC
			bilayerFESC<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvals,eigenvects);
#elif USE_BAFEAS
			BaFeAs<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvects);
			for (size_t i=0;i<nbands;i++) eigenvects(i,i) -= param.mu;
			eigen(eigenvals,eigenvects);
			return;
#elif USE_KFE2SE2
			KFe2Se2<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvects);
			for (size_t i=0;i<nbands;i++) eigenvects(i,i) -= param.mu;
			eigen(eigenvals,eigenvects);
			return;
#elif USE_FOURORBITAL
			FourOrbital<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
			s.getBands(k,eigenvects);
			for (size_t i=0;i<nbands;i++) eigenvects(i,i) -= param.mu;
			eigen(eigenvals,eigenvects);
			return;
#else
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
			return;
#endif
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


		void calcBandStructure(std::string file,bool printOccupation) {
			size_t nktot(kmesh.nktot);
			RangeType range(0,nktot,conc);
			std::vector<std::vector<FieldType> > ek(nktot,VectorType(nbands,0));
			std::vector<MatrixType> weights(nktot,MatrixType(nbands,nbands));
			ComplexMatrixType ak(nbands,nbands);

			VectorType occupation(nktot,0);

			for (;!range.end();range.next()) {
				size_t ik = range.index();
				std::vector<FieldType> k(3);
				kmesh.momenta.getRow(ik,k);
				getEkAndAk(k,ek[ik],ak);
				if (printOccupation) 
					for (size_t i=0;i<nbands;i++) occupation[ik] += fermi(ek[ik][i],1./param.temperature);

				for (size_t iband=0; iband<nbands; iband++) 
					for (size_t iorb=0; iorb<nbands; iorb++) weights[ik](iorb,iband)=abs(ak(iorb,iband));
			}

			for (size_t ik=0;ik<nktot;ik++) {
				conc.reduce(ek[ik]);
				conc.reduce(weights[ik]);
			}
			if (printOccupation) {
				conc.reduce(occupation);
				FieldType occ(0.0);
				for (size_t ik=0;ik<nktot;ik++) occ += occupation[ik];
				occ /= FieldType(nktot);
				if (conc.rank()==0) std::cout << "\n\t\tFilling = " << 2 * occ  << "\n\n";
			}

			if (conc.rank()==0) {
				const char *filename = file.c_str();
				std::ofstream os(filename);
				int precision=5;
				os.precision(precision);
				for (size_t ik=0; ik<nktot; ik++) {
					std::vector<FieldType> k(3);
					kmesh.momenta.getRow(ik,k);
					os << std::setw(15);
					os << k[0]/param.pi_f;
					os << std::setw(15);
					os << k[1]/param.pi_f;
					os << std::setw(15);
					os << k[2]/param.pi_f;
					os << std::setw(15);
					for (size_t i = 0; i < nbands; ++i) os << ek[ik][i] << std::setw(15);

					for (size_t iorb=0; iorb < nbands; iorb++) 
						for (size_t iband=0; iband < nbands; iband++) os << weights[ik](iorb,iband) << std::setw(15);
					
					os << "\n";
				}
				os.close();
			}

		}

		inline FieldType fermi(const FieldType& energy, const FieldType& invT) {
			FieldType xx=energy*invT;
			FieldType f;
			if (xx <= 0.) f = 1./(1.+exp(xx));
			else if (xx > 0.) f = exp(-xx)/(1.+exp(-xx));
			else f = 0.5;
			return f;
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

		void eigen(VectorType& eigenvals, ComplexMatrixType& matrix) const {
			int n = matrix.n_row()
;			int lwork = 2*n-1;
			std::vector<ComplexType> work(lwork);
			std::vector<FieldType> rwork(3*n-2);
			int info;

			GEEV('V','U',n,matrix,n,eigenvals,work,lwork,rwork,&info);
			if (info!=0) {
				throw std::runtime_error("zheev: failed\n");
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
			// Lplus(0,1) = I*sqrt(3.);
			// Lplus(0,2) =  -sqrt(3.);
			// Lplus(1,0) =   sqrt(3.);
			// Lplus(1,3) = I         ;
			// Lplus(1,4) = 1.0       ;
			// Lplus(2,0) = I*sqrt(3.);
			// Lplus(2,3) = 1.0       ;
			// Lplus(2,4) = -I        ;
			// Lplus(3,1) =  1.0      ;
			// Lplus(3,2) =  I        ;
			// Lplus(4,1) = -I        ;
			// Lplus(4,2) =  1.0      ;

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

		void diag2x2(const FieldType& eplus, const FieldType& eminus, const FieldType& exy,
				 VectorType& w, ComplexMatrixType& v) {

			FieldType sign((exy >=0) - (exy < 0));
			FieldType squareRoot(sqrt(pow(eminus,2)+pow(exy,2)));
			w[0] = eplus - squareRoot;
			w[1] = eplus + squareRoot;

			FieldType uk,vk;
			coherence_factors(eminus,exy,uk,vk);

			v(0,0) = -vk;
			v(1,0) = sign*uk;
			v(0,1) = sign*uk;
			v(1,1) = vk;

			return;
		}

		void coherence_factors(const FieldType& e1,
							   const FieldType& e2,
							   FieldType& uk,
							   FieldType& vk) {

			FieldType E(sqrt(pow(e1,2)+pow(e2,2)));
			if (E==0) {
				uk = 1.0; vk = 1.0;
				return;
			} else {
				uk = sqrt(0.5*(1.0+e1/E));
				vk = sqrt(0.5*(1.0-e1/E));
				return;
			}

		}



	};

}

#endif
