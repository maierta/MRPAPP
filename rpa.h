//-*-C++-*-

#ifndef RPA_H
#define RPA_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"


namespace rpa {

	extern "C" void 
#ifdef glyph
	zgemm_
#else
	zgemm
#endif
	(char *,char *,int *,int *,int *,std::complex<double> *,
						  std::complex<double> *,int *,std::complex<double> *,
						  int *,std::complex<double> *,std::complex<double> *,int *);
	
	extern "C" void 
#ifdef glyph
	zgetrf_
#else
	zgetrf
#endif
	(int *,int *,std::complex<double> *,
						   int *,int *,int *);
	
	extern "C" void 
#ifdef glyph
	zgetri_
#else
	zgetri
#endif
	(int *,std::complex<double> *, 
						   int *, int *, 
						   std::complex<double> *, int *, int *);

	inline void GEMM(char transa,char transb,int m,int n,int k,
					 std::complex<double> &alpha,
				     psimag::Matrix<std::complex<double> > &a,int lda,
				     psimag::Matrix<std::complex<double> > &b,int ldb,
				     std::complex<double> &beta,
				     psimag::Matrix<std::complex<double> > &c,int ldc) 
					 {
#ifdef glyph						
						zgemm_
#else
						zgemm
#endif
						(&transa,&transb,&m,&n,&k,&alpha,&(a(0,0)),
							  &lda,&(b(0,0)),&ldb,&beta,&(c(0,0)),&ldc);
						}
	inline void GETRF(int m, int n, 
	 			      psimag::Matrix<std::complex<double> > &a,int lda, 
	 			      std::vector<int> &ipiv, int *info) 
					 {
#ifdef glyph						
						zgetrf_
#else
						zgetrf
#endif
						(&m,&n,&(a(0,0)),&lda,&(ipiv[0]),info);
						}
	inline void GETRI(int m, 
				      psimag::Matrix<std::complex<double> > &a, int lda, 
				      std::vector<int> &ipiv, 
					  psimag::Matrix<std::complex<double> > &work, int lwork, int *info) 
					 {
#ifdef glyph						
						zgetri_
#else
						zgetri
#endif
						(&m,&(a(0,0)),&lda,&(ipiv[0]),&(work(0,0)),&lwork,info);
						}


	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class interaction {

		private:
			typedef std::complex<Field>				ComplexType;
			typedef MatrixTemplate<Field> 			MatrixType;
			typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
			const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
			Field U,Up,J,Jp;
			size_t nOrb,msize;
			// std::vector<std::vector<Field> >  orbPos;
			ComplexType ii;

		public:
			ComplexMatrixType spinMatrix;
			ComplexMatrixType chargeMatrix;


		interaction(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters):
			param(parameters),
			U(param.U),
			Up(param.Up),
			J(param.J),
			Jp(param.Jp),
			nOrb(param.nOrb),
			msize(size_t(nOrb*nOrb)),
			// orbPos(nOrb,std::vector<Field>(3,0.0)),
			ii(1.0,0.0),
			spinMatrix(msize,msize),
			chargeMatrix(msize,msize)
		{
			setupInteractionMatrix();
		}


		void setupInteractionMatrix() {
			if (param.Case== "YBCO_orthoII_perpStripes" 
				|| param.Case == "1band" 
				|| param.Case == "bilayer"
				|| param.Case == "trilayer"
				|| param.Case == "Checkerboard"
				|| param.Case == "bilayer_Harr_Seb"
				|| param.Case == "BSCCObilayer_OD"
						) {
				// Only diagonal terms
				for (size_t l1 = 0; l1 < nOrb; ++l1) {
					size_t ind1(l1+l1*nOrb);
					spinMatrix  (ind1,ind1)   = U;
					chargeMatrix(ind1,ind1)   = -U;
				}	

			} else if (param.Case=="EmeryOnlyUd") {
				// Only diagonal 11 terms
				for (size_t l1 = 0; l1 < 1; ++l1) {
					size_t ind1(l1+l1*nOrb);
					spinMatrix  (ind1,ind1)   = U;
					chargeMatrix(ind1,ind1)   = -U;
				}	
			} else {  // general multi-orbital model
				size_t limit(nOrb);
				if (param.sublattice==1) limit=nOrb<10?nOrb/2:5;
				std::cout << "In rpa.h: limit=" << limit << "\n";

				
				ComplexMatrixType spinSubMatrix(limit*limit,limit*limit);
				ComplexMatrixType chargeSubMatrix(limit*limit,limit*limit);
					
					// First the diagonal terms
					for (size_t l1 = 0; l1 < limit; ++l1) {
							for (size_t l2 = 0; l2 < limit; ++l2) {
								size_t ind1 = l2+l1*limit;
								if (l1==l2) {
									spinSubMatrix(ind1,ind1)   = U+param.deltaU[l1];
									chargeSubMatrix(ind1,ind1) = -U-param.deltaU[l1];;
									} else {
									spinSubMatrix(ind1,ind1)   = Up;
									chargeSubMatrix(ind1,ind1) = Up-2*J;
									}
							}
						}	
					// Off-diagonal terms
					for (size_t l1=0; l1 < limit; l1++) {
						size_t ind1 = l1+l1*limit;
						for (size_t l2=0; l2 < limit; l2++) {
							size_t ind2 = l2+l2*limit;
							if (l2!=l1) {
								spinSubMatrix(ind1,ind2)   = J;
								chargeSubMatrix(ind1,ind2) = -2.*Up+J;
							}
						}
					}
					// Finally the pair hopping terms
					for (size_t l1=0; l1 < limit; l1++) {
						for (size_t l2=0; l2 < limit; l2++) {
							size_t ind1 = l2+l1*limit;
							size_t ind2 = l1+l2*limit;
							if (l2!=l1) {
								spinSubMatrix(ind1,ind2)   = Jp;
								chargeSubMatrix(ind1,ind2) = -Jp;
							}
						}
					}
				if (param.sublattice==0) {
					for (size_t i=0; i<msize; i++) for (size_t j=0; j<msize; j++) {
						spinMatrix(i,j) = spinSubMatrix(i,j);
						chargeMatrix(i,j) = chargeSubMatrix(i,j);
					}
					} else {
						for(size_t l1=0; l1<limit; l1++) for (size_t l2=0; l2<limit; l2++) {
						for(size_t l3=0; l3<limit; l3++) for (size_t l4=0; l4<limit; l4++) {
							
							size_t ind1=l2+l1*limit;
							size_t ind2=l4+l3*limit;
		
							size_t ind3=l2+l1*nOrb;
							size_t ind4=l4+l3*nOrb;

							spinMatrix(ind3,ind4) = spinSubMatrix(ind1,ind2);
							chargeMatrix(ind3,ind4) = chargeSubMatrix(ind1,ind2);

							ind3=l2+limit+(l1+limit)*nOrb; // position of 2. Fe d-orbitals is shifted by limit wrt 1. Fe d-orbs.
							ind4=l4+limit+(l3+limit)*nOrb; // position of 2. Fe d-orbitals is shifted by limit wrt 1. Fe d-orbs.

							spinMatrix(ind3,ind4) = spinSubMatrix(ind1,ind2);
							chargeMatrix(ind3,ind4) = chargeSubMatrix(ind1,ind2);
						}
						}	
					}					
				}
		}

		void calcRPAResult(ComplexMatrixType& matrix0, 
						   ComplexMatrixType& interactionMatrix, 
						   ComplexMatrixType& matrix1, 
						   std::vector<Field> q=std::vector<Field>(3,0.0)) {
			int n = interactionMatrix.n_row();
			// int m = matrix0.n_col();
			// int k = spinMatrix.n_col();
			ComplexType alpha(1.0);
			ComplexType beta(0.0);
			std::vector<int> ipiv(n);
			int info;
			ComplexMatrixType work(n,n);
			ComplexMatrixType c(n,n);
			int lwork(n);

			// std::cout << "spinMatrix: " << "\n" << spinMatrix;

			// GEMM('N','N',m,n,k,alpha,spinMatrix,n,matrix0,m,beta,c,n);
			matMul(interactionMatrix,matrix0,c);
			// Result of matrix multiplication is in c
			// std::cout << "matrix0: " << "\n" << matrix0;
			for (size_t i = 0; i < c.n_row(); ++i) for (size_t j = 0; j < c.n_col(); ++j) c(i,j) = -c(i,j);
			for (size_t i = 0; i < c.n_row(); ++i) c(i,i) = 1.+c(i,i);
			// Now invert
			GETRF(n,n,c,n,ipiv,&info);
			if (info!=0) throw std::runtime_error("GETRF: failed\n");
			GETRI(n,c,n,ipiv,work,lwork,&info);
			if (info!=0) throw std::runtime_error("GETRI: failed\n");
			// Now multiply result with matrix0
			// std::cout << "inv(c): " << "\n" << c;
			// GEMM('N','N',m,n,k,alpha,matrix0,m,c,n,beta,matrix1,n);
			matMul(matrix0,c,matrix1);
			// std::cout << "matrix1: " << "\n" << matrix1;

		}

	};








}

#endif
