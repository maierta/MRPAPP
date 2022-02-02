// Model file for Sr2RuO4 with spin-orbit coupling --> 6 bands total
#ifndef SRRUO_SO_3D_H
#define SRRUO_SO_3D_H


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
		typedef std::complex<Field>	        ComplexType;
		typedef MatrixTemplate<ComplexType>	ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 				FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		
	private:
		size_t msize;
		// std::vector<std::vector<Field> > kField;
		// std::vector<size_t> bandField;
		// std::vector<std::vector<Field> > gapField(3, 0);
	public:
		FieldType nbands;
		ComplexMatrixType   spinMatrix;
		ComplexMatrixType   chargeMatrix;
		std::vector<int>    spinOfEll;
		std::vector<size_t> orbOfEll;

		std::vector<FieldType> kxGap;
		std::vector<FieldType> kyGap;
		std::vector<FieldType> kzGap;
		std::vector<FieldType> DeltaGap1;
		std::vector<FieldType> DeltaGap2;
		std::vector<FieldType> DeltaGap3;

		size_t nTotal;

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

			// if (param.options.find("calcSus")!=std::string::npos && param.scState) {
			// 	if (conc.rank()==0) {
			// 		std::cout << "Now reading in GapFile\n";
			// 		readGapFromFile();
			// 		std::cout << "GapFile was read in \n";
			// 	}

			// 	conc.broadcast(nTotal);
			// 	if (conc.rank()!=0) {
			// 		kxGap.resize(nTotal);
			// 		kyGap.resize(nTotal);
			// 		kzGap.resize(nTotal);
			// 		DeltaGap1.resize(nTotal);
			// 		DeltaGap2.resize(nTotal);
			// 		DeltaGap3.resize(nTotal);
			// 	}
			// 	conc.broadcast(kxGap);
			// 	conc.broadcast(kyGap);
			// 	conc.broadcast(kzGap);
			// 	conc.broadcast(DeltaGap1);
			// 	conc.broadcast(DeltaGap2);
			// 	conc.broadcast(DeltaGap3);
			// 	if (conc.rank()==0) std::cout << "... and broadcast\n";
			// }


		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& H0) {

			FieldType ekXZ, ekYZ, ekXY;
			FieldType gxzyz, Txzxy, Tyzxy;
			FieldType cx,cy,cxy,c2x,c2y,sx,sy,cz2,cx2,cy2,sx2,sy2,sz2;

			sx = sin(k[0]); sy = sin(k[1]); sx2 = sin(0.5*k[0]); sy2 = sin(0.5*k[1]); sz2 = sin(0.5*k[2]);
			cx = cos(k[0]); cy = cos(k[1]); cxy = cos(k[0])*cos(k[1]);
			c2x = cos(2*k[0]); c2y = cos(2*k[1]); 
			cx2 = cos(0.5*k[0]); cy2 = cos(0.5*k[1]); cz2 = cos(0.5*k[2]);


			FieldType t1,t2,t3,t4,t5,t11,t12,tint;

			t1  = 0.088; t2 = 0.009; t3 = 0.080; t4 = 0.040; t5 = 0.005; 
			t11 = 0.003; t12 = -0.005; tint = -0.004;
			// t11 = 0.0; t12 = 0.0; tint = 0.0;


			gxzyz   = -4*t11*sx*sy - 4*t12*sx2*sy2*cz2;
			Txzxy   = -4*tint*cx2*sy2*sz2; 
			Tyzxy   = -4*tint*sx2*cy2*sz2; 

			ekXZ    = -2*t1*cx-2*t2*cy - param.mu;
			ekYZ    = -2*t2*cx-2*t1*cy - param.mu;
			ekXY    = -2*t3*(cx+cy)-4*t4*cxy-2*t5*(c2x+c2y) - param.mu;

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

				if (param.Case == "Cobo0") {
						t12z = 0; t56z = 0; tdxy = 0; td = 0; 			
				} else if (param.Case == "Cobo1") {
						t12z = 0.005; t56z = 0.003; tdxy = 0; td = 0; 			
				} else if (param.Case == "Cobo2") {
						t12z = 0.005; t56z = 0.004; tdxy = 0; td = 0; 			
				} else if (param.Case == "Cobo3") {
						t12z = 0.005; t56z = 0.005; tdxy = 0; td = 0; 
				} else if (param.Case == "Cobo4") {
						t12z = 0.005; t56z = 0.003; tdxy = 0.002; td = 0.002; 
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

		// template<typename FieldType>
		// void readSCGap(std::vector<std::vector<FieldType> > & kField,
		// 			      std::vector<size_t>& bandField,
		// 			      std::vector<FieldType>& gapField,
		// 			      const std::string& filename) {
                //
		// 	std::vector<FieldType> data;
		// 	loadVector(data,filename);
		// 	// We assume that each line has the format qx,qy,qz,chi_{l1,l2,l3,l4}(q) chi_phys(q)
		// 	size_t step = 5; //  kx,ky,kz, band index, gap
		// 	size_t nq(data.size()/step);
		// 	std::cout << "We have " << nq << " points in the stored gap file " << filename.c_str() << "\n";
		// 	std::cout << "data.size()" << data.size() << "\n";
		// 	// if (kField.size()!=nq) {std::cerr << "Number of k-points in file not correct! \n"; exit(1);}
		// 	for (size_t iq=0; iq<nq; iq++) {
		// 		std::vector<FieldType> k(3, 0);
		// 		k[0] = data[0 + iq*step];
		// 		k[1] = data[1 + iq*step];
		// 		k[2] = data[2 + iq*step];
		// 		kField.push_back(k);
		// 		bandField.push_back(data[3 + iq*step]);
		// 		gapField.push_back(data[4 + iq*step]);
		// 		}
                //
		// 		// std::cout << "iq,q,chi.calcSus,file.calcSus: " << iq << ", " << qField[iq] << ", "
		// 			  // << chiField[iq].calcSus() << ", " << data[2*msize*msize + 3 + iq*step] << "\n";
		// 	}

		std::complex<Field> calcSCGapReadIn(VectorType& k, size_t band, ComplexMatrixType& Uk) {
			// return ComplexType(0.0,0.0); // uses function in gaps3D.h directly for now
			// std::cout << "kx:"<<k[0]<<", ky:"<<k[1]<<", kz:"<< k[2]<<"\n";
			FieldType delta=1.0e-5;
			for (size_t ik=0; ik < kxGap.size(); ik++) {
				bool x = (fabs(k[0]-kxGap[ik]) < delta); 
				bool y = (fabs(k[1]-kyGap[ik]) < delta); 
				bool z = (fabs(k[2]-kzGap[ik]) < delta); 
				if (x & y & z) {
					if (band==0 || band==1) {
						return ComplexType(DeltaGap1[ik],0);
					} else if (band==2 || band==3) {
						return ComplexType(DeltaGap2[ik]);
					} else if (band==4 || band==5) {
						return ComplexType(DeltaGap3[ik]);
					} 
				}
			}
			std::cout << "k-point not found in gap file !!! \n";
			exit(0);
		}

		std::complex<Field> calcSCGap(VectorType& k, size_t band, ComplexMatrixType& Uk) {
			ComplexType Delta;
			const ComplexType ii = ComplexType(0.0,1.0);
			if (param.gAmpl == "") {
				ComplexType g3, g4;
				FieldType sx, sy, s2x, s2y, d1, d2;
				sx = sin(k[0]);
				sy = sin(k[1]);
				s2x = sin(2*k[0]);
				s2y = sin(2*k[1]);
				d1 = sin(k[0])*cos(k[1]);
				d2 = sin(k[1])*cos(k[0]);

				if (band==0 || band==1) {
					g3 = -0.2410*sx - 0.1925*sy - 0.090*s2x - 0.072*s2y - 0.1270*d1 - 0.1013*d2;
					g4 = -0.2410*sy + 0.1925*sx - 0.090*s2y + 0.072*s2x - 0.1270*d2 + 0.1013*d1;
				} else if (band==2 || band==3) {
					g3 = -0.1114*sx - 0.089*sy - 0.0073*s2x - 0.0057*s2y - 0.1370*d1 - 0.1095*d2;
					g4 = -0.1114*sy + 0.089*sx - 0.0073*s2y + 0.0057*s2x - 0.1370*d2 + 0.1095*d1;
				} else if (band==4 || band==5) {
					g3 = 0.0425*sx + 0.0340*sy - 0.0227*s2x - 0.0182*s2y + 0.0170*d1 + 0.0136*d2;
					g4 = 0.0425*sy - 0.0240*sx - 0.0227*s2y + 0.0182*s2x + 0.0170*d2 - 0.0136*d1;
				} 
				Delta = param.Delta0*(g3 + ii*g4) * sin(k[2]/2);

			} else if (param.gAmpl == "SrRuO_A1g") {
				FieldType cxs, cxy;
				cxs = cos(k[0]) + cos(k[1]);
				cxy = cos(k[0]) * cos(k[1]);

				if (band==0 || band==1) {
					Delta = 0.3700 + 0.2454*cxs - 0.0564*cxy;
				} else if (band==2 || band==3) {
					Delta = 0.5483 + 0.6702*cxs + 0.6870*cxy;
				} else if (band==4 || band==5) {
					Delta = 0.3023 + 0.6283*cxs + 0.9289*cxy;
				} 
			} else if (param.gAmpl == "SrRuO_B1g") {
				FieldType cd, cd2;
				cd  = cos(k[0]) - cos(k[1]);
				cd2 = cos(2*k[0]) - cos(2*k[1]);

				if (band==0 || band==1) {
					Delta = -0.8236*cd - 0.3012*cd2;
				} else if (band==2 || band==3) {
					Delta = -0.0045*cd - 0.0685*cd2;
				} else if (band==4 || band==5) {
					Delta = 0.0145*cd - 0.0743*cd2;
				} 
			}

			return Delta * param.Delta0;
		}

		void readGapFromFile() {
			std::string file = param.scGapfile;
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format dx,dy,dz,orb1,orb2,t
			size_t length = 6;
			// if (param.complexHopping) length=7;
			size_t nLinesTotal(data.size()/length);
			nTotal = nLinesTotal;
			if (conc.rank()==0) std::cout << "gapfile contains " << nLinesTotal << " lines\n";
			for (size_t i = 0; i < nLinesTotal; i++) {
					kxGap.push_back    (data[i*length]);
					kyGap.push_back    (data[i*length+1]);
					kzGap.push_back    (data[i*length+2]);
					DeltaGap1.push_back(data[i*length+3]);
					DeltaGap2.push_back(data[i*length+4]);
					DeltaGap3.push_back(data[i*length+5]);
			}
			size_t nLines = kxGap.size();
			if (conc.rank()==0) std::cout << nLines <<" entries used from gap input file\n";
		}


	};
}

#endif
