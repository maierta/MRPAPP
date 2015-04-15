#ifndef BANDS_H
#define BANDS_H


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
	class Bands {
	private:
		typedef MatrixTemplate<Field> 	MatrixType;
		typedef std::complex<Field>		ComplexType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      VectorType;
		typedef Field 					FieldType;
		const rpa::parameters<Field,MatrixTemplate>& param;
		ConcurrencyType& conc;
		size_t nk;
		size_t dim;
		int nLines;
		std::vector<Field> dx,dy,dz,ht;
		std::vector<size_t> orb1,orb2;


		
	public:
		FieldType nbands;
		MatrixType bandsEk;
		std::vector<ComplexMatrixType> bandsAk;

		Bands(const rpa::parameters<Field,MatrixTemplate>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			nk(param.nkInt),
			dim(param.dimension),
			nbands(param.nOrb),
			bandsEk(0,0),
			bandsAk(0)
		
		{
			if (param.tbfile!="") readCSVFile();
			if (param.sublattice==1) fixdr();
		}
		
		void calcBandsKMesh(const size_t nk, const size_t nkz, const size_t dim,
			                const FieldType& kxmin, const FieldType& kxmax,
			                const FieldType& kymin, const FieldType& kymax,
			                const FieldType& kzmin, const FieldType& kzmax,
			                const size_t saveBands=0) {
			
			std::cout << "Calculating band structure \n";

			std::vector<FieldType> ek(param.nOrb);
			psimag::Matrix<std::complex<FieldType> > ak(param.nOrb,param.nOrb);
			momentumDomain<FieldType,psimag::Matrix> momentumDomain1(param,nk,nkz,dim);
			if (dim==3) momentumDomain1.set_momenta(kxmin,kxmax,kymin,kymax,kzmin,kzmax);
			else if(dim==2) momentumDomain1.set_momenta(kxmin,kxmax,kymin,kymax);
			if (conc.rank()==0) std::cout << "momentumDomain.nktot="<<momentumDomain1.nktot <<"\n";
			size_t nktot(momentumDomain1.nktot);
			VectorType occupation(nktot,0);
			// FieldType occupation(0.);

			if (saveBands==1) {
				bandsEk.resize(nktot,param.nOrb);
				bandsAk.resize(nktot);
				for (size_t ik=0;ik<nktot;ik++) bandsAk[ik]=ComplexMatrixType(param.nOrb,param.nOrb);
			}


			MatrixType ekFull(nktot,param.nOrb);

			typedef PsimagLite::Range<ConcurrencyType> RangeType;
			RangeType range(0,nktot,conc);
			if (conc.rank()==0) std::cout << "range: "<< range.end() << "\n";
			for (;!range.end();range.next()) {
				size_t ik = range.index();
				// if (conc.rank()==0) std::cout << ik << " of " << momentumDomain1.nktot << " total momenta done \n";
				std::vector<FieldType> k(3);
				momentumDomain1.momenta.getRow(ik,k);
				calcBandsK(k,ek,ak);
				// if (conc.rank()==0 && k[2]==0.0) std::cout << "k=(" << k[0] << " , " << k[1] << " , " << k[2] << "; ek2 = " << ek[1] << "\n"; 
				for (size_t i = 0; i < param.nOrb; ++i) {
					ekFull(ik,i) = ek[i];
					bandsEk(ik,i) = ek[i];
					for (size_t j=0;j<param.nOrb;j++) bandsAk[ik](i,j) = ak(i,j);
				}
				
				for (int i=0;i<nbands;i++) occupation[i] += fermi(ek[i],1./param.temperature);
			}
			conc.reduce(ekFull);
			conc.reduce(occupation);

			if (conc.rank()==0) {
				std::ofstream os("ek.dat");
				// os.flags(std::ios::left);
				// os.fill(' ');
				int precision=5;
				// int spacer=precision+10;
				os.precision(precision);
				// os.width(spacer);
				for (size_t ik=0; ik<nktot; ik++) {
					std::vector<FieldType> k(3);
					momentumDomain1.momenta.getRow(ik,k);
					os << std::setw(15);
					os << k[0]/param.pi_f;
					// os.width(spacer);
					os << std::setw(15);
					os << k[1]/param.pi_f;
					// os.width(spacer);
					os << std::setw(15);
					os << k[2]/param.pi_f;
					// os.width(spacer);
					os << std::setw(15);
					for (size_t i = 0; i < param.nOrb; ++i) os << ekFull(ik,i) << std::setw(15);
					os << "\n";
				}
				os.close();
			}

			FieldType occ(0.0);
			for (size_t i=0;i<nktot;i++) occ += occupation[i];
			occ /= FieldType(nktot*param.sublattice*2);
			if (conc.rank()==0) std::cout << "Filling = " << 2.*occ << "\n";

			if (conc.rank()==0) std::cout << "Done calculating band structure \n";

		}

		inline void calcBandsK1(const VectorType& k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
			FieldType t=1.0;
			FieldType mu=0.0;
			eigenvals[0] = -2.*t*(cos(k[0])+cos(k[1])) - mu;
			eigenvects(0,0) = 1.0;
		}

		inline void calcBandsK(const VectorType& k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
			// const ComplexType ii = ComplexType(0.0,1.0);
			if (param.Case=="") {
				FieldType exponent(0.);
				int n = eigenvects.n_col();
				std::vector<FieldType> ks(3,0);
				if(param.kTrafo==1) {
					ks[0] = k[0]*param.WTrafo(0,0)+k[1]*param.WTrafo(1,0)+k[2]*param.WTrafo(2,0);
					ks[1] = k[0]*param.WTrafo(0,1)+k[1]*param.WTrafo(1,1)+k[2]*param.WTrafo(2,1);
					ks[2] = k[0]*param.WTrafo(0,2)+k[1]*param.WTrafo(1,2)+k[2]*param.WTrafo(2,2);
				} else {
					ks[0] = k[0]; ks[1] = k[1]; ks[2] = k[2];
				}
				for (int i=0;i<n;i++) for (int j=0;j<n;j++) eigenvects(i,j)=0.;

				for (int i = 0; i < nLines; ++i)
				{
					if (orb1[i]<=orb2[i]) {
						exponent = (ks[0]*dx[i] + ks[1]*dy[i] + ks[2]*dz[i]); 
						eigenvects(orb1[i],orb2[i]) += ht[i] * ComplexType(cos(exponent),sin(exponent));
					}
				}

				// if (k[0]==-0.5*param.pi_f && k[1]==param.pi_f && k[2]==0.0) {
				// 	std::cout << "H: " << eigenvects << "\n";
				// }
				diag(eigenvals,eigenvects);
				for (size_t i=0;i<nbands;i++) eigenvals[i] -= param.mu;
			} else if (param.Case=="SrRuO") {
				SrRuO<FieldType,MatrixTemplate,ConcurrencyType> s(param,conc);
				s.getBands(k,eigenvals,eigenvects);
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

		FieldType fermiVelocity2D(VectorType& k, size_t band) {
			
			std::vector<FieldType> ekpx(param.nOrb);
			std::vector<FieldType> ekpy(param.nOrb);
			std::vector<FieldType> ekmx(param.nOrb);
			std::vector<FieldType> ekmy(param.nOrb);
			psimag::Matrix<std::complex<FieldType> > ak(param.nOrb,param.nOrb);


			FieldType dk(0.0001);
			VectorType kpx(3,0);
			VectorType kpy(3,0);
			VectorType kmx(3,0);
			VectorType kmy(3,0);

		    kpx[0]=k[0]+dk ; kpx[1]=k[1]; kpx[2]=k[2];
		    kmx[0]=k[0]-dk ; kmx[1]=k[1]; kmx[2]=k[2];
		    kpy[0]=k[0]    ; kpy[1]=k[1]+dk; kpy[2]=k[2];
		    kmy[0]=k[0]    ; kmy[1]=k[1]-dk; kmy[2]=k[2];

			calcBandsK(kpx,ekpx,ak);
			calcBandsK(kpy,ekpy,ak);
			calcBandsK(kmx,ekmx,ak);
			calcBandsK(kmy,ekmy,ak);
	           
		    FieldType rx=(ekpx[band]-ekmx[band])/(2.*dk);
		    FieldType ry=(ekpy[band]-ekmy[band])/(2.*dk);
    
		    FieldType vel = sqrt(pow(rx,2)+pow(ry,2));

		    // std::cout << "vel in fermiVelocity2D:"<< vel << ","<<rx<<","<<ry<<"\n";

		    return vel;
		}

		void readCSVFile() {
			std::string file = param.tbfile;
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format dx,dy,dz,orb1,orb2,t
			nLines=data.size()/6;
			for (int i = 0; i < nLines; i++)
			{
				dx.push_back  (data[i*6]);
				dy.push_back  (data[i*6+1]);
				dz.push_back  (data[i*6+2]);
				orb1.push_back(size_t(data[i*6+3])-1);
				orb2.push_back(size_t(data[i*6+4])-1);
				ht.push_back  (data[i*6+5]);
 			}
		}

		void fixdr() {
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

		void diag(VectorType& eigenvals, ComplexMatrixType& matrix)
		{
			int n = matrix.n_row();
			int lwork = 2*n-1;
			std::vector<ComplexType> work(lwork);
			std::vector<FieldType> rwork(3*n-2);
			int info;

			GEEV('V','U',n,matrix,n,eigenvals,work,lwork,rwork,&info);
			if (info!=0) {
				throw std::runtime_error("zheev: failed\n");
			}
		}

		
	};

}

#endif





















