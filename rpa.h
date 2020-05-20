//-*-C++-*-

#ifndef RPA_H
#define RPA_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"
// #include "tbFromFile.h"
// #include "SrRuO_SO.h"
// #include "1band_wSpin.h"
// #include "coupledLadders.h"
#include "model.h"


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
			typedef Field					FieldType;
			typedef std::complex<Field>			ComplexType;
			typedef MatrixTemplate<Field> 			MatrixType;
			typedef MatrixTemplate<ComplexType> 		ComplexMatrixType;
			const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
			ConcurrencyType& conc;
			Field U,Up,J,Jp;
			size_t nOrb,msize;
			// std::vector<std::vector<Field> >  orbPos;
			ComplexType ii;

		public:
			model<FieldType, MatrixTemplate, ConcurrencyType> model;

// 			// Model specific needed for charge and spin matrices
// #ifdef USE_SRRUO
// 			SrRuO_SO<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_1BANDWSPIN
// 			SingleBand_wSpin<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_BILAYER
// 			orthoIIBilayer<FieldType,MatrixTemplate,ConcurrencyType> model;
// 			// bilayer<FieldType,MatrixTemplate,ConcurrencyType> s;
// #elif USE_BILAYER_1BAND
// 			bilayer<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_BSCCObilayer
// 			BSCCObilayer<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_BILAYER_FESC
// 			bilayerFESC<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_BAFEAS
// 			BaFeAs<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_KFE2SE2
// 			KFe2Se2<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_FOURORBITAL
// 			FourOrbital<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_TBFILE
// 			tbFromFile<FieldType,MatrixTemplate,ConcurrencyType> model;
// #elif USE_COUPLEDLADDERS
// 			coupledLadders<FieldType,MatrixTemplate,ConcurrencyType> model;
// #endif
			// ComplexMatrixType spinMatrix;
			// ComplexMatrixType chargeMatrix;


		interaction(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			U(param.U),
			Up(param.Up),
			J(param.J),
			Jp(param.Jp),
			nOrb(param.nOrb),
			msize(size_t(nOrb*nOrb)),
			// orbPos(nOrb,std::vector<Field>(3,0.0)),
			ii(1.0,0.0),
			model(param,conc)
		{
			model.setupInteractionMatrix();
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
