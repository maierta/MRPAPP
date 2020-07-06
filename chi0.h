//-*-C++-*-

#ifndef CHI0_H
#define CHI0_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "Fermi.h"
#include "bandstructure.h"
#include "susInt.h"
#include "gaps2D.h"
#include "gaps3D.h"
#include "sepBasis.h"
#include "utilities.h"
#include "model.h"

namespace rpa {

	// class to hold the susceptibility matrix and some functions
	template<typename Field,template<typename> class MatrixTemplate,
	         typename ConcurrencyType>
	class susc:
		public MatrixTemplate<std::complex<Field> >
	{
	private:
		rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;

	public:
 		typedef psimag::Matrix<std::complex<Field> > BaseType;
		size_t nOrb,msize;

	private:
		susc() {}

	public:
		susc(const susc& other):
			BaseType (other),
			param(other.param),
			conc(other.conc),
			nOrb(other.nOrb),
			msize(other.msize)
		{}

		susc& operator=(const susc& other) {
			BaseType::operator=(other);
			nOrb = other.nOrb;
			msize = other.msize;
			return (*this);
		}

		susc(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			 ConcurrencyType& concurrency):
		BaseType (parameters.nOrb*parameters.nOrb,parameters.nOrb*parameters.nOrb),
		param(parameters),
		conc(concurrency),
		nOrb(param.nOrb),
		msize(nOrb*nOrb)
		{

		}

		void setLowerTriangle() {
			for (size_t i = 0; i < msize; ++i) for (size_t j = 0; j < i; ++j) {
				(*this)(i,j) = conj((*this)(j,i));
			}
		}

		std::complex<Field> calcSus() const {
			std::complex<Field> chiPhys(0.0);
			for (size_t l1 = 0; l1 < nOrb; ++l1) {
				for (size_t l2 = 0; l2 < nOrb; ++l2){
					size_t ind1(l1+l1*nOrb);
					size_t ind2(l2+l2*nOrb);
					chiPhys += 0.5*(*this)(ind1,ind2);
				}
			}
			Field factor(1.0);
			// if (param.sublattice==1) factor=param.nSitesPerUnitCell;
			return chiPhys/factor;

		}

		void allReduce() {
			conc.allReduce((*this));
		}

		void print() {
			for (size_t ind1 = 0; ind1 < msize; ++ind1) {
				for (size_t ind2 = 0; ind2 < msize; ++ind2){
					std::cout << ind1<<","<<ind2<<" , " << (*this)(ind1,ind2);
				}
			}

		}

	};

	// class to calculate the susceptibility matrix
	template<typename Field, typename SuscType, typename BandsType,
			 typename GapType,
	         template<typename> class MatrixTemplate,
	         typename ConcurrencyType>
	class calcChi0Matrix {

	private:
		typedef Field 							FieldType;
		typedef std::complex<Field>				ComplexType;
		typedef MatrixTemplate<Field> 			MatrixType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      		VectorType;
		typedef std::vector<ComplexType>      	ComplexVectorType;
		typedef PsimagLite::Range<ConcurrencyType> RangeType;

		rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		const momentumDomain<Field,MatrixTemplate,ConcurrencyType>& kmesh;
		size_t nktot;
		BandsType& bands;

		ConcurrencyType& conc;
		size_t nOrb;
		size_t msize;
		FieldType invT;
		bool kMap_;
		bool calcOnlyDiagonal_;

		VectorType chi0k;

	public:



		void printChi0k(const VectorType& q) {
			std::ofstream os("chi0k.txt");
			// std::ofstream os2("chi0kpq.txt");
			int width(10);
			os.precision(width);
			// os2.precision(width);
			os << std::fixed;
			// os2 << std::fixed;
			for (size_t ik = 0; ik < kmesh.nktot; ++ik)	{
				std::vector<FieldType> k(3);
				std::vector<FieldType> kq(3);
				kmesh.momenta.getRow(ik,k);
				for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
				kmesh.mapTo1BZ(kq);
				os << k[0]/param.pi_f << " , " << k[1]/param.pi_f << " , " << k[2]/param.pi_f << " , " << chi0k[ik] << "\n";
				// os2 << kq[0]/param.pi_f << " , " << kq[1]/param.pi_f << " , " << kq[2]/param.pi_f << " , " << chi0k[ik] << "\n";
			}
		}



		// Constructor for w=0, Delta=0 calculation
		calcChi0Matrix(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			const momentumDomain<Field,psimag::Matrix,ConcurrencyType>& kmeshIn,
			const VectorType& q,
			BandsType& bandsIn,
			ConcurrencyType& concurrency,
			SuscType& chi0matrix,
			const bool& cacheB=0):

		param(parameters),
		kmesh(kmeshIn),
		nktot(kmesh.nktot),
		bands(bandsIn),
		conc(concurrency),
		nOrb(param.nOrb),
		msize(nOrb*nOrb),
		invT(1./param.temperature)

		{
			for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++)
					chi0matrix(i,j) = ComplexType(0.0,0.0);

			ComplexMatrixType Sm(nOrb,nOrb);
			VectorType *ekq(&bands.ekq[0]);
			VectorType *ek(&bands.ek[0]);
			ComplexMatrixType *Mkq(&bands.Mkq[0]);
			ComplexMatrixType *Mk(&bands.Mk[0]);

			for (size_t ik = 0; ik < nktot; ++ik)	{
				if (cacheB) {
					ekq = &bands.ekq[ik]; ek = &bands.ek[ik];
					Mkq = &bands.Mkq[ik]; Mk = &bands.Mk[ik];
				} else { // Need to calculate band tensors first
					std::vector<FieldType> k(3);
					std::vector<FieldType> kq(3);
					kmesh.momenta.getRow(ik,k);
					for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
					bands.calculateBandTensors(k, bands.ek[0], bands.ak[0], bands.Mk[0], bands.gapk[0], bands.MkFF[0], param.scState, 0);
					bands.calculateBandTensors(kq, bands.ekq[0], bands.akq[0], bands.Mkq[0], bands.gapkq[0], bands.MkqFF[0], param.scState, 1);
				}
				// Now build S-matrix
				for (size_t band1 = 0; band1 < nOrb; ++band1){
					for (size_t band2 = 0; band2 < nOrb; ++band2){
						Sm(band1,band2) = susInt((*ekq)[band1],(*ek)[band2],invT);
					}
				}

				// Now build Mkq(l1,l3,b1) * Sm(b1,b2) * Mk(l4,l2,b2)
				ComplexMatrixType c0(Mkq->n_row(),Sm.n_col());
				ComplexMatrixType c1(c0.n_row(),Mk->n_col());
				matMul((*Mkq),Sm,c0);
				matMul(c0,(*Mk),c1);

				for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
					size_t l1 = param.indexToOrb(i,0); size_t l3 = param.indexToOrb(i,1);
					size_t l4 = param.indexToOrb(j,0); size_t l2 = param.indexToOrb(j,1);
					size_t ind1(param.OrbsToIndex(l1,l2));
					size_t ind2(param.OrbsToIndex(l3,l4));
					chi0matrix(ind1,ind2) += c1(i,j);
				}

			}
			for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
					chi0matrix(i,j) /= ComplexType(kmesh.nktot,0.0);
					// std::cout << chi0matrix(i,j) << "\n";
				}
		}

		// Constructor for finite w calculation, pre- band diagonalization
		calcChi0Matrix(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			const momentumDomain<Field,psimag::Matrix,ConcurrencyType>& kmeshIn,
			const VectorType& q,
			BandsType& bandsIn,
			ConcurrencyType& concurrency,
			SuscType& chi0matrix,
			const FieldType& omega,
			const bool& cacheB=0):

		param(parameters),
		kmesh(kmeshIn),
		nktot(kmesh.nktot),
		bands(bandsIn),
		conc(concurrency),
		nOrb(param.nOrb),
		msize(nOrb*nOrb),
		invT(1./param.temperature)

		{
			for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++)
					chi0matrix(i,j) = ComplexType(0.0,0.0);

			ComplexMatrixType SmGG(nOrb,nOrb);
			ComplexMatrixType SmFF(nOrb,nOrb);
			VectorType *ekq, *ek;
			ComplexVectorType *gapkq, *gapk;
			ComplexMatrixType *Mkq, *Mk, *MkqFF, *MkFF;

			ekq = &bands.ekq[0]; ek = &bands.ek[0];
			gapkq = &bands.gapkq[0]; gapk = &bands.gapk[0];
			Mkq = &bands.Mkq[0]; Mk = &bands.Mk[0];
			MkqFF = &bands.MkqFF[0]; MkFF = &bands.MkFF[0];

			for (size_t ik = 0; ik < nktot; ++ik)	{

				if (cacheB) { // Band tensors are already pre-calculated and we just need to point to it
					ekq = &bands.ekq[ik]; ek = &bands.ek[ik];
					gapkq = &bands.gapkq[ik]; gapk = &bands.gapk[ik];
					Mkq = &bands.Mkq[ik]; Mk = &bands.Mk[ik];
					MkqFF = &bands.MkqFF[ik]; MkFF = &bands.MkFF[ik];
				} else { // Need to calculate band tensors first
					std::vector<FieldType> k(3);
					std::vector<FieldType> kq(3);
					kmesh.momenta.getRow(ik,k);
					for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
					bands.calculateBandTensors(k, bands.ek[0], bands.ak[0], bands.Mk[0], bands.gapk[0], bands.MkFF[0], param.scState, 0);
					bands.calculateBandTensors(kq, bands.ekq[0], bands.akq[0], bands.Mkq[0], bands.gapkq[0], bands.MkqFF[0], param.scState, 1);
				}

				// First build S-matrix
				for (size_t band1 = 0; band1 < nOrb; ++band1){
					for (size_t band2 = 0; band2 < nOrb; ++band2){
						if (!param.scState) {
							SmGG(band1,band2) = susInt((*ekq)[band1],(*ek)[band2],invT,omega,param.damp);
						} else {
							SmGG(band1,band2) = susIntBCSGG((*ekq)[band1],(*ek)[band2],(*gapkq)[band1],(*gapk)[band2],invT,omega,param.damp);
							SmFF(band1,band2) = susIntBCSFF((*ekq)[band1],(*ek)[band2],(*gapkq)[band1],(*gapk)[band2],invT,omega,param.damp,param.parity);
							// SmGG(band1,band2) = susIntBCS((*ekq)[band1],(*ek)[band2],(*gapkq)[band1],(*gapk)[band2],invT,omega,param.damp,param.signF);
						}
					}
				}

				// Now build Mkq(l1,l3,b1) * SmGG(b1,b2) * Mk(l4,l2,b2)
				ComplexMatrixType c0(Mkq->n_row(),SmGG.n_col());
				ComplexMatrixType c1(c0.n_row(),Mk->n_col());
				ComplexMatrixType c2(c0.n_row(),Mk->n_col());
				matMul((*Mkq),SmGG,c0);
				matMul(c0,(*Mk),c1);

				if (param.scState) {
					// Now build MkqFF(l1,l3,b1) * SmFF(b1,b2) * MkFF(l4,l2,b2)
					matMul((*MkqFF),SmFF,c0);
					matMul(c0,(*MkFF),c2);
					// matMul((*Mkq),SmFF,c0);
					// matMul(c0,(*Mk),c2);
				}

				for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
					size_t l1 = param.indexToOrb(i,0); size_t l3 = param.indexToOrb(i,1);
					size_t l4 = param.indexToOrb(j,0); size_t l2 = param.indexToOrb(j,1);
					size_t ind1 = param.OrbsToIndex(l1,l2);
					size_t ind2 = param.OrbsToIndex(l3,l4);
					chi0matrix(ind1,ind2) += c1(i,j);
					// if (param.scState) chi0matrix(ind1,ind2) += c2(i,j);

					if (param.scState) { // Add FF term
						l4 = param.indexToOrb(i,0); l1 = param.indexToOrb(i,1);
						l2 = param.indexToOrb(j,0); l3 = param.indexToOrb(j,1);
						size_t ind1 = param.OrbsToIndex(l1,l2);
						size_t ind2 = param.OrbsToIndex(l3,l4);
						chi0matrix(ind1,ind2) += c2(i,j);
					}
				}

			}
			for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
					chi0matrix(i,j) /= ComplexType(kmesh.nktot,0.0);
					// std::cout << chi0matrix(i,j) << "\n";
				}

		}

		// Constructor for w=0, Delta=0 calculation, band diagonalization on the fly
		// calcChi0Matrix(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
		// 	const momentumDomain<Field,psimag::Matrix,ConcurrencyType>& kmeshIn,
		// 	BandsType& bandsIn,
		// 	const std::vector<Field>& q,
		// 	ConcurrencyType& concurrency,
		// 	SuscType& chi0matrix,
		// 	const bool kMap=0,
		// 	const bool calcOnlyDiagonal=0):
                //
		// param(parameters),
		// kmesh(kmeshIn),
		// bands(bandsIn),
		// conc(concurrency),
		// nOrb(param.nOrb),
		// msize(nOrb*nOrb),
		// invT(1./param.temperature),
		// kMap_(kMap),
		// calcOnlyDiagonal_(calcOnlyDiagonal),
		// chi0k(kMap?kmeshIn.nktot:0,0.0)
                //
		// {
		// 	VectorType k(3),kq(3);
		// 	VectorType ek(nOrb,0),ekq(nOrb,0);
		// 	ComplexMatrixType ak(nOrb,nOrb), akq(nOrb,nOrb), susTerm(nOrb,nOrb);
                //
		// 	for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++)
		// 			chi0matrix(i,j) = ComplexType(0.0,0.0);
                //
		// 	for (size_t ik = 0; ik < kmesh.nktot; ++ik)	{
		// 		kmesh.momenta.getRow(ik,k);
		// 		for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
		// 		bands.getEkAndAk(k,ek,ak);
		// 		bands.getEkAndAk(kq,ekq,akq);
		// 		for (size_t band1 = 0; band1 < nOrb; ++band1){
		// 			for (size_t band2 = 0; band2 < nOrb; ++band2){
		// 				ComplexType r1 = ComplexType(susInt(ekq[band1],ek[band2],invT),0);
		// 				// r1 = susInt(ekq[band1],ek[band2],invT,FieldType(0.0),param.damp);
                //
		// 				for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++) {
		// 					size_t l1 = param.indexToOrb(i,0); size_t l2 = param.indexToOrb(i,1);
		// 					size_t l3 = param.indexToOrb(j,0); size_t l4 = param.indexToOrb(j,1);
                //
		// 					// Check if l1,l2 and l3,l4 are on the same site
		// 					if ((param.orbToSite[l1]!=param.orbToSite[l2]) || (param.orbToSite[l3]!=param.orbToSite[l4])) continue;
		// 					if (calcOnlyDiagonal && ((l1!=l2) || (l3!=l4))) continue;
		// 					// ComplexType c1 = ak(l4,band2)  * conj(ak(l2,band2))
		// 			                       // * akq(l1,band1) * conj(akq(l3,band1));
                //
		// 					ComplexType c1 = computeM(l1,l2,l3,l4,band1,band2,ak,akq);
		// 					// if (l1==0 && l2==0 && l3==0 && l4==0 && band1==band2) std::cout << "band1,band2,c1: " << k[0]<<","<<k[1]<<","<<band1<<","<<band2<<","<<c1<<","<<r1 << ","<<ekq[band1]<<","<<ek[band2]<<"\n";
		// 					// if (l1==l2 && l3==l4) std::cout << "l1,l2,l3,l4,band1,band2,c1: " << l1 << " , " << l2 << " , " << l3 << " , " << l4 << " , " << k[0]<<","<<k[1]<<","<<band1<<","<<band2<<"," << r1 << " , " <<c1<<"," <<ekq[band1]<<","<<ek[band2]<<"\n";
                //
		// 					chi0matrix(i,j) += r1*c1 ;
		// 					if (kMap_ && l1==l2 && l3==l4) chi0k[ik] += real(c1*r1);
		// 					// if (i==0 && j==0) std::cout << "chi0matrix " << chi0matrix(i,j) << "\n";
		// 				}
		// 			}
		// 		}
		// 	}
		// 	for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++) {
		// 			chi0matrix(i,j) /= ComplexType(kmesh.nktot,0.0);
		// 			// std::cout << i <<","<<j<<","<<chi0matrix(i,j);
		// 		}
		// 	chi0matrix.setLowerTriangle();
		// 	if (kMap_) printChi0k(q);
                //
		// }


		// Constructor for fininte w, finite Delta calculation, band diagonalization on the fly ek,ak 
		// calcChi0Matrix(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
		// 	const momentumDomain<Field,psimag::Matrix,ConcurrencyType>& kmeshIn,
		// 	BandsType& bandsIn,
		// 	const std::vector<Field>& q,
		// 	ConcurrencyType& concurrency,
		// 	SuscType& chi0matrix,
		// 	GapType& DeltaIn,
		// 	const FieldType& omega=FieldType(0.0),
		// 	const bool kMap=0):
		// param(parameters),
		// kmesh(kmeshIn),
		// bands(bandsIn),
		// conc(concurrency),
		// nOrb(param.nOrb),
		// msize(nOrb*nOrb),
		// invT(1./param.temperature),
		// kMap_(kMap),
		// chi0k(kMap?kmeshIn.nktot:0,0.0)
                //
		// {
		// 	GapType& Delta(DeltaIn);
		// 	VectorType k(3),kq(3);
		// 	VectorType ek(nOrb,0),ekq(nOrb,0);
		// 	ComplexMatrixType ak(nOrb,nOrb), akq(nOrb,nOrb), susTerm(nOrb,nOrb);
                //
                //
		// 	for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++)
		// 			chi0matrix(i,j) = ComplexType(0.0,0.0);
                //
		// 	for (size_t ik = 0; ik < kmesh.nktot; ++ik)	{
		// 		kmesh.momenta.getRow(ik,k);
		// 		for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
		// 		bands.getEkAndAk(k,ek,ak);
		// 		bands.getEkAndAk(kq,ekq,akq);
		// 		kmesh.mapTo1BZ(kq); // make sure kq is in 1. BZ because gap may not be periodic in kz
		// 		for (size_t band1 = 0; band1 < nOrb; ++band1){
		// 			ComplexType gap1 = Delta(kq,band1,akq);
		// 			// gap1 *= pow(param.Omega0,2)/(pow(ekq[band1],2)+pow(param.Omega0,2)); // Lorentzian cut-off
		// 			gap1 *= exp(-std::pow(ekq[band1],2) / std::pow(param.Omega0,2)); // Gaussian cut-off
		// 			for (size_t band2 = 0; band2 < nOrb; ++band2){
		// 				ComplexType r1(0.0);
		// 				ComplexType gap2 = Delta(k,band2,ak);
		// 				// gap2 *= pow(param.Omega0,2)/(pow(ek[band2],2)+pow(param.Omega0,2)); // Lorentzian cut-off
		// 				gap2 *= exp(-std::pow(ek[band2],2) / std::pow(param.Omega0,2)); // Gaussian cut-off
		// 				r1 = susIntBCS(ekq[band1],ek[band2],gap1,gap2,invT,omega,param.damp,param.signF);
                //
                //
		// 				for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
		// 					size_t l1 = param.indexToOrb(i,0); size_t l2 = param.indexToOrb(i,1);
		// 					size_t l3 = param.indexToOrb(j,0); size_t l4 = param.indexToOrb(j,1);
                //
		// 					// ComplexType c1 = ak(l4,band2)  * conj(ak(l2,band2))
		// 					// 	           * akq(l1,band1) * conj(akq(l3,band1));
                //
		// 					// Check if l1,l2 and l3,l4 are on the same site
		// 					if ((param.orbToSite[l1]!=param.orbToSite[l2]) || (param.orbToSite[l3]!=param.orbToSite[l4])) continue;
                //
		// 					ComplexType c1 = computeM(l1,l2,l3,l4,band1,band2,ak,akq);
                //
		// 					// std::cout << "k,r1===== " << k << "," <<r1 <<  "\n";
                //
                //
		// 					chi0matrix(i,j) += c1*r1 ;
                //
		// 					if (kMap_ && l1==l2 && l3==l4) chi0k[ik] += imag(c1*r1);
		// 				}
		// 			}
		// 		}
		// 	}
		// 	for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++)
		// 			chi0matrix(i,j) /= ComplexType(kmesh.nktot,0.0);
		// 	// chi0matrix.setLowerTriangle();
                //
		// 	if (kMap_) printChi0k(q);
		// }

		// // Constructor for finite w, Delta=0 calculation, band diagonalization on the fly
		// calcChi0Matrix(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
		// 	const momentumDomain<Field,psimag::Matrix,ConcurrencyType>& kmeshIn,
		// 	BandsType& bandsIn,
		// 	const std::vector<Field>& q,
		// 	ConcurrencyType& concurrency,
		// 	SuscType& chi0matrix,
		// 	const FieldType& omega,
		// 	const bool kMap=0):
		// param(parameters),
		// kmesh(kmeshIn),
		// bands(bandsIn),
		// conc(concurrency),
		// nOrb(param.nOrb),
		// msize(nOrb*nOrb),
		// invT(1./param.temperature),
		// kMap_(kMap),
		// chi0k(kMap?kmeshIn.nktot:0,0.0)
                //
		// {
		// 	VectorType k(3),kq(3);
		// 	VectorType ek(nOrb,0),ekq(nOrb,0);
		// 	ComplexMatrixType ak(nOrb,nOrb), akq(nOrb,nOrb), susTerm(nOrb,nOrb);
                //
                //
		// 	for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++)
		// 			chi0matrix(i,j) = ComplexType(0.0,0.0);
                //
		// 	for (size_t ik = 0; ik < kmesh.nktot; ++ik)	{
		// 		kmesh.momenta.getRow(ik,k);
		// 		for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
		// 		bands.getEkAndAk(k,ek,ak);
		// 		bands.getEkAndAk(kq,ekq,akq);
		// 		for (size_t band1 = 0; band1 < nOrb; ++band1){
		// 			for (size_t band2 = 0; band2 < nOrb; ++band2){
		// 				ComplexType r1(0.0);
		// 				r1 = susInt(ekq[band1],ek[band2],invT,omega,param.damp);
                //
		// 				// if (imag(r1)>=1.e-5) std::cout << band1 << "," << band2 << "," << ekq[band1]<<","<<ek[band2]<<","<<imag(r1) << "\n";
                //
		// 				// for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++) {
		// 				for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
		// 					size_t l1 = param.indexToOrb(i,0); size_t l2 = param.indexToOrb(i,1);
		// 					size_t l3 = param.indexToOrb(j,0); size_t l4 = param.indexToOrb(j,1);
                //
		// 					// ComplexType c1 = ak(l4,band2)  * conj(ak(l2,band2))
		// 					// 	           * akq(l1,band1) * conj(akq(l3,band1));
                //
		// 					// Check if l1,l2 and l3,l4 are on the same site
		// 					if ((param.orbToSite[l1]!=param.orbToSite[l2]) || (param.orbToSite[l3]!=param.orbToSite[l4])) continue;
                //
		// 					ComplexType c1 = computeM(l1,l2,l3,l4,band1,band2,ak,akq);
                //
		// 					chi0matrix(i,j) += c1*r1 ;
		// 					// chi0matrix(i,j) += r1 ;
		// 					if (kMap_ && l1==l2 && l3==l4) chi0k[ik] += imag(c1*r1);
                //
		// 				}
		// 			}
		// 		}
		// 	}
		// 	// for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++) {
		// 	for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) {
		// 			chi0matrix(i,j) /= ComplexType(kmesh.nktot,0.0);
		// 			// std::cout << chi0matrix(i,j) << "\n";
		// 		}
		// 	// chi0matrix.print();
		// 	// chi0matrix.setLowerTriangle(); // For finite omega there is a damping term w + ii*damp which breaks this symmetry
		// 	if (kMap_) printChi0k(q);
                //
		// }




		// Emery (3-band) d-p-p model


		// Constructor for w=0  calculation for Emery problem with k,k',q-dependent interactions
		// band diagonalization on the fly
		calcChi0Matrix(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			const momentumDomain<Field,psimag::Matrix,ConcurrencyType>& kmeshIn,
			BandsType& bandsIn,
			const std::vector<Field>& q,
			ConcurrencyType& concurrency,
			ComplexMatrixType& chi0,
			ComplexMatrixType& chi0_g,
			ComplexMatrixType& chi0_gg):
		param(parameters),
		kmesh(kmeshIn),
		bands(bandsIn),
		conc(concurrency),
		nOrb(param.nOrb),
		msize(nOrb*nOrb),
		invT(1./param.temperature)

		{
			VectorType k(3),kq(3);
			VectorType ek(nOrb,0),ekq(nOrb,0);
			ComplexMatrixType ak(nOrb,nOrb), akq(nOrb,nOrb), susTerm(nOrb,nOrb);

			for (size_t i=0; i<chi0.n_row(); i++) for (size_t j=0; j<chi0.n_col(); j++)
				chi0(i,j) = ComplexType(0.0);
			for (size_t i=0; i<chi0_gg.n_row(); i++) for (size_t j=0; j<chi0_gg.n_col(); j++)
				chi0_gg(i,j) = ComplexType(0.0);
			for (size_t i=0; i<chi0_g.n_row(); i++) for (size_t j=0; j<chi0_g.n_col(); j++)
				chi0_g(i,j) = ComplexType(0.0);

			for (size_t ik=0; ik<kmesh.nktot; ik++) { // loop over k-points
				kmesh.momenta.getRow(ik,k);
				sepBasis<FieldType,psimag::Matrix,ConcurrencyType> basis(param,conc,k);
			for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
				bands.getEkAndAk(k,ek,ak);
				bands.getEkAndAk(kq,ekq,akq);
				for (size_t band1=0; band1<nOrb; ++band1) for (size_t band2=0; band2<nOrb; ++band2)  { // loop over band indices
					ComplexType r1 = ComplexType(susInt(ekq[band1],ek[band2],invT),0); // f(e_k+q-f_ek/(ek+q-ek))
					for (size_t i=0; i<msize; i++) for (size_t j=0; j<msize; j++)  { // loop over orbital indices
						size_t l1 = param.indexToOrb(i,0); size_t l2 = param.indexToOrb(i,1);
						size_t l3 = param.indexToOrb(j,0); size_t l4 = param.indexToOrb(j,1);
						ComplexType me = computeM(l1,l2,l3,l4,band1,band2,ak,akq);
						ComplexType integrand = r1*me;

						if ((l1 == l2) && (l3 == l4)) chi0(l1,l3) += integrand; // chi0_l1,l2 only diag. elements
						for (size_t iB=0; iB<19; iB++) for (size_t jB=0; jB<19; jB++)  // loop over basis set
							chi0_gg(iB,jB) += integrand*basis(iB,l4,l3)*basis(jB,l2,l1);
							// chi0_gg(iB,jB) += integrand*basis(jB,l4,l3)*basis(iB,l2,l1);
							// chi0_gg(iB,jB) += integrand*basis(iB,l1,l2)*basis(jB,l3,l4);
						if (l3==l4) for (size_t iB=0; iB<19; iB++) {
							chi0_g(iB,l3) -= integrand*basis(iB,l2,l1);
							// chi0_g(iB,l3) -= integrand*basis(iB,l1,l2);
							// if (iB==3 && l3 == 6) std::cout << "integrand: " << integrand << "," << basis(iB,l2,l1) << "\n";
						}
					}
				}
			}

			for (size_t i=0;i<chi0.n_row();i++) for (size_t j=0;j<chi0.n_col();j++)
				chi0(i,j) /= ComplexType(kmesh.nktot,0.0);
			for (size_t i=0;i<chi0_g.n_row();i++) for (size_t j=0;j<chi0_g.n_col();j++)
				chi0_g(i,j) /= ComplexType(kmesh.nktot,0.0);
			for (size_t i=0;i<chi0_gg.n_row();i++) for (size_t j=0;j<chi0_gg.n_col();j++)
				chi0_gg(i,j) /= ComplexType(kmesh.nktot,0.0);
		}



	inline ComplexType computeM(size_t l1,size_t l2,size_t l3,size_t l4, size_t band1,size_t band2,
								ComplexMatrixType& ak, ComplexMatrixType& akq) {

		// ComplexType c1 = ak(l4,band2)  * conj(ak(l2,band2))
		// 	           * akq(l1,band1) * conj(akq(l3,band1));  see Kreisel et al., PRB 88, 094522 (2013), Eqs. (4) and (5)

		const FieldType& ar(real(ak(l4,band2)));
		const FieldType& ai(imag(ak(l4,band2)));
		const FieldType& br(real(ak(l2,band2)));
		const FieldType& bi(imag(ak(l2,band2)));
		const FieldType& cr(real(akq(l1,band1)));
		const FieldType& ci(imag(akq(l1,band1)));
		const FieldType& dr(real(akq(l3,band1)));
		const FieldType& di(imag(akq(l3,band1)));

		FieldType t1r(ar*br+ai*bi);
		FieldType t1i(ai*br-ar*bi);
		FieldType t2r(cr*dr+ci*di);
		FieldType t2i(ci*dr-cr*di);

		ComplexType c1(t1r*t2r-t1i*t2i,t1r*t2i+t1i*t2r);

		return c1;
	}

	};

	// Class to calculate chi0 matrix on a q-mesh
	template<typename Field, typename SuscType, typename BandsType,
	         typename GapType,
	         template<typename> class MatrixTemplate,
	         typename ConcurrencyType>
	class chi0q:
		public std::vector<SuscType> {
	private:
		typedef Field 							FieldType;
		typedef std::complex<Field>				ComplexType;
		typedef MatrixTemplate<Field> 			MatrixType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      		VectorType;
		typedef PsimagLite::Range<ConcurrencyType> RangeType;

		rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		momentumDomain<Field,psimag::Matrix,ConcurrencyType>& qMesh;
		std::string file;
		size_t nOrb, msize;

		model<FieldType, MatrixTemplate, ConcurrencyType> model;

	public:
 		typedef std::vector<SuscType> BaseType;

		chi0q(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			momentumDomain<Field,psimag::Matrix,ConcurrencyType>& qMeshIn,
			ConcurrencyType& concurrency):
			BaseType (qMeshIn.nktot,SuscType(parameters,concurrency)),
			param(parameters),
			conc(concurrency),
			qMesh(qMeshIn),
			file("none"),
			nOrb(param.nOrb),
			msize(nOrb*nOrb),
			model(param,conc)
		{
			calcChi0q();
		}

		chi0q(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			momentumDomain<Field,psimag::Matrix,ConcurrencyType>& qMeshIn,
			std::string fileIn,
			ConcurrencyType& concurrency):
			BaseType (qMeshIn.nktot,SuscType(parameters,concurrency)),
			param(parameters),
			conc(concurrency),
			qMesh(qMeshIn),
			file(fileIn),
			nOrb(param.nOrb),
			msize(nOrb*nOrb),
			model(param,conc)
		{
			if (qMesh.nktot > 0) {
				// if (file=="none") calcChi0q();
				if (param.readChi==0) calcChi0q();
				else {
					if(conc.rank()==0) readChiqTxt(file,qMesh);
						// conc.broadcast(qMesh.nk);
						// conc.broadcast(qMesh.nkz);
						// conc.broadcast(qMesh.nktot);
						conc.broadcast(qMesh.momenta);
						for (size_t iq=0;iq<qMesh.nktot;iq++) conc.broadcast((*this)[iq]);
				}
			}
		}

		void calcChi0q() {
			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			kmesh.set_momenta(false);
			BandsType bands(param,model,conc,kmesh,true);
			RangeType range(0,qMesh.nktot,conc);
			for (;!range.end();range.next()) {

				size_t iq = range.index();

				std::vector<FieldType> q(3);
				qMesh.momenta.getRow(iq,q);

				calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
				             calcChi0(param,kmesh,q,bands,conc,(*this)[iq],param.cacheBands);

				// for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) (*this)[iq](i,j) = calcChi0(i,j);
				// (*this)[iq].setLowerTriangle();

				if (conc.rank()==0) {
					std::cout.precision(7);
					std::cout << "iq = " << iq << " q= " << q << " of " << qMesh.nktot
	                          << " total. ChiPhys=" << (*this)[iq].calcSus()
	                          << "\n";
                }
			}
			// for (size_t iq=0;iq<qMesh.nktot;iq++) conc.allReduce((*this)[iq]);
			for (size_t iq=0;iq<qMesh.nktot;iq++) (*this)[iq].allReduce();
			if (conc.rank()==0) {
				std::cout << "Now printing out chiq \n";
				writeChiqTxt();
			}
		}

		void writeChiqTxt() {
			std::ofstream os("susOfQFull.txt");
			std::ofstream os2("susRPA.txt");
			int width(10);
			os.precision(width);
			os2.precision(width);
			os << std::fixed;
			os2 << std::fixed;
			os << qMesh.nk << " , " << qMesh.nkz << "\n";
			std::vector<FieldType> q(3);
			for (size_t iq=0;iq<qMesh.nktot;iq++) {
				qMesh.momenta.getRow(iq,q);
				os << q[0] << " , " << q[1] << " , " << q[2] << " , ";
				for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
					os << real((*this)[iq](l1,l2))<< " , ";
				}
				for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
					os << imag((*this)[iq](l1,l2))<< " , ";
				}
     			ComplexType sus0((*this)[iq].calcSus());
				os << real(sus0) << " , " << imag(sus0);
				os << "\n";
				SuscType chiRPA(param,conc);
     			calcRPAResult((*this)[iq],model.spinMatrix,chiRPA,q);
     			ComplexType susR(chiRPA.calcSus());
     			os2 << q[0] << " , " << q[1] << " , " << q[2] << " , ";
     			os2 << real(susR) << ","  << imag(susR) << "\n";
			}

		}

		void readChiqTxt(const std::string& file,momentumDomain<Field,psimag::Matrix,ConcurrencyType>& qMesh) {
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format qx,qy,qz,chi_{l1,l2,l3,l4}(q) chi_phys(q)
			size_t step=3+msize*msize*2+2; // 3 qx,qy,qz, 2*msize*msize for real and imag. part of chi, 2 for chiPhys
			// qMesh.nk = size_t(data[0]);
			// qMesh.nkz = size_t(data[1]);
			// qMesh.nktot = qMesh.nk*qMesh.nk*qMesh.nkz;
			size_t nqTot(qMesh.nktot);
			std::cout << "nqTot,data.size()="<<qMesh.nk << ","<< qMesh.nkz << " , "<< data.size() << "\n";
			for (size_t iq = 0; iq < nqTot; iq++)
			{
				qMesh.momenta(iq,0) = data[2+iq*step];
				qMesh.momenta(iq,1) = data[2+iq*step+1];
				qMesh.momenta(iq,2) = data[2+iq*step+2];
				for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++)
					(*this)[iq](l1,l2) = data[2+iq*step+3+l2+l1*msize]; // real part
				for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++)
					(*this)[iq](l1,l2) = ComplexType(real((*this)[iq](l1,l2)),data[2+iq*step+3+msize*msize+l2+l1*msize]); // imag. part
				std::vector<FieldType> q(3);
				qMesh.momenta.getRow(iq,q);
				std::cout << "q=" << qMesh.momenta(iq,0) << ","<<qMesh.momenta(iq,1) << "," << qMesh.momenta(iq,2) << "; chiq=" << (*this)[iq].calcSus() << "\n";
			}
			std::cout << "nqTot,nq,nqz=" << qMesh.nktot <<","<<qMesh.nk<<","<<qMesh.nkz<< "\n";
		}

	};

	template<typename FieldType, typename SuscType>
	void writeChiqTxt2(const std::vector<std::vector<FieldType> >& qField,
					   const std::vector<SuscType>& chiField,
					   const std::string& filename) {
		std::ofstream os(filename.c_str());
		int width(10);
		os.precision(width);
		os << std::fixed;
		size_t nq(qField.size());
		size_t msize(chiField[0].msize);
		for (size_t iq=0;iq<nq;iq++) {
			os << qField[iq][0] << " , " << qField[iq][1] << " , " << qField[iq][2] << " , ";
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				os << real(chiField[iq](l1,l2))<< " , ";
			}
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				os << imag(chiField[iq](l1,l2))<< " , ";
			}
			os << real(chiField[iq].calcSus());
			os << "\n";
		}
	}

	template<typename FieldType, typename SuscType>
	void readChiqTxt2(std::vector<std::vector<FieldType> >& qField,
				      std::vector<SuscType>& chiField,
				      const std::string& filename) {

		std::vector<FieldType> data;
		typedef std::complex<FieldType> ComplexType;
		loadVector(data,filename);
		// We assume that each line has the format qx,qy,qz,chi_{l1,l2,l3,l4}(q) chi_phys(q)
		size_t msize(chiField[0].msize);
		size_t step=3+msize*msize*2+1; // 3 qx,qy,qz, 2*msize*msize for real and imag. part of chi, 1 for chiPhys
		size_t nq(data.size()/step);
		std::cout << "We have " << nq << " points in the stored file " << filename.c_str() << "\n";
		std::cout << "data.size()" << data.size() << "\n";
		if (qField.size()!=nq) {std::cerr << "Number of q-points in file not correct! \n"; exit(1);}
		for (size_t iq=0; iq<nq; iq++) {
			qField[iq][0] = data[0 + iq*step];
			qField[iq][1] = data[1 + iq*step];
			qField[iq][2] = data[2 + iq*step];
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				FieldType r1 = data[l2+l1*msize + 3 + iq*step];
				FieldType i1 = data[msize*msize + l2+l1*msize + 3 + iq*step];
				chiField[iq](l1,l2) = ComplexType(r1,i1);
			}

			// std::cout << "iq,q,chi.calcSus,file.calcSus: " << iq << ", " << qField[iq] << ", "
			          // << chiField[iq].calcSus() << ", " << data[2*msize*msize + 3 + iq*step] << "\n";
		}
	}


}

#endif
