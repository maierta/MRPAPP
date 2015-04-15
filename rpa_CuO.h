//-*-C++-*-

#ifndef RPA_CUO_H
#define RPA_CUO_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"


namespace rpa {


	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class interactionEmery { 

		private:
			typedef std::complex<Field>				ComplexType;
			typedef MatrixTemplate<Field> 			MatrixType;
			typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
			typedef std::vector<Field> 				VectorType;
			const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
			Field U_d_c,U_d_s,U_p_c,U_p_s,U_pd_c,U_pd_s,U_pp_c,U_pp_s,U_d_coupl,U_p_coupl,U_pd_coupl,U_pp_coupl;
			size_t nOrb,msize;

		public:
			// ComplexMatrixType interactionMatrix;
			ComplexMatrixType V_X_c;
			ComplexMatrixType V_X_s;
			ComplexMatrixType V_D;
			ComplexMatrixType V_X_coupl;
			ComplexMatrixType V_D_coupl;
			ComplexMatrixType V_Charge;
			ComplexMatrixType V_Spin;
			ComplexMatrixType V_Charge_coupl;
			ComplexMatrixType V_Spin_coupl;


		interactionEmery(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters):
			param(parameters),
			U_d_c(param.U_d_c),
			U_d_s(param.U_d_s),
			U_p_c(param.U_p_c),
			U_p_s(param.U_p_s),
			U_pd_c(param.U_pd_c),
			U_pd_s(param.U_pd_s),
			U_pp_c(param.U_pp_c),
			U_pp_s(param.U_pp_s),
			U_d_coupl(param.U_d_coupl),
			U_p_coupl(param.U_p_coupl),
			U_pd_coupl(param.U_pd_coupl),
			U_pp_coupl(param.U_pp_coupl),
			nOrb(param.nOrb),
			msize(size_t(nOrb*nOrb)),
			// interactionMatrix(19,19),
			V_X_c(19,19),
			V_X_s(19,19),
			V_D(19,19),
			V_X_coupl(19,19),
			V_D_coupl(19,19),
			V_Charge(19,19),
			V_Spin(19,19),
			V_Charge_coupl(19,19),
			V_Spin_coupl(19,19)

		{
			setupExchangeMatrix();
		}



		void setupExchangeMatrix() {
				// Now set up U^s and U^c matrices (0=O px; 1=O py; 2=Cu d)
				// V_X is diagonal
				for (size_t i=0;i<V_X_c.n_row();++i) for (size_t j=0;j<V_X_c.n_row();++j) {
					V_X_c(i,j) = 0.0; V_X_s(i,j) = 0.0;  V_X_coupl(i,j) = 0.0; 
					V_D(i,j) = 0.0; V_D_coupl(i,j) = 0.0;
				}
				for (int i = 0; i < 4; ++i)   {V_X_c(i,i) = 2.*U_pd_c; V_X_s(i,i) = 2.*U_pd_s; V_X_coupl(i,i) = 2.*U_pd_coupl;}
				for (int i = 11; i < 15; ++i) {V_X_c(i,i) = 2.*U_pd_c; V_X_s(i,i) = 2.*U_pd_s; V_X_coupl(i,i) = 2.*U_pd_coupl;}
				for (int i = 4; i < 8; ++i)   {V_X_c(i,i) = 2.*U_pp_c; V_X_s(i,i) = 2.*U_pp_s; V_X_coupl(i,i) = 2.*U_pp_coupl;}
				for (int i = 15; i < 19; ++i) {V_X_c(i,i) = 2.*U_pp_c; V_X_s(i,i) = 2.*U_pp_s; V_X_coupl(i,i) = 2.*U_pp_coupl;}
				for (int i = 9; i < 11; ++i)  {V_X_c(i,i) =     U_p_c; V_X_s(i,i) =     U_p_s; V_X_coupl(i,i) =     U_p_coupl;}
				V_X_c(8,8) = U_d_c; V_X_s(8,8) = U_d_s; V_X_coupl(8,8) = U_d_coupl; 
		}

		void setupDirectMatrix(std::vector<Field>& q) {
				Field cx(cos(0.5*q[0]));
				Field cy(cos(0.5*q[1]));

				V_D(8,8)   = U_d_c;
				V_D(9,9)   = U_p_c;
				V_D(10,10) = U_p_c;
				V_D(8,9)   = 2.*U_pd_c*cx;
				V_D(9,8)   = 2.*U_pd_c*cx;
				V_D(8,10)  = 2.*U_pd_c*cy;
				V_D(10,8)  = 2.*U_pd_c*cy;
				V_D(9,10)  = 4.*U_pp_c*cx*cy;
				V_D(10,9)  = 4.*U_pp_c*cx*cy;

				V_D_coupl(8,8)   = U_d_coupl;
				V_D_coupl(9,9)   = U_p_coupl;
				V_D_coupl(10,10) = U_p_coupl;
				V_D_coupl(8,9)   = 2.*U_pd_coupl*cx;
				V_D_coupl(9,8)   = 2.*U_pd_coupl*cx;
				V_D_coupl(8,10)  = 2.*U_pd_coupl*cy;
				V_D_coupl(10,8)  = 2.*U_pd_coupl*cy;
				V_D_coupl(9,10)  = 4.*U_pp_coupl*cx*cy;
				V_D_coupl(10,9)  = 4.*U_pp_coupl*cx*cy;
			}

		void setupVBare(std::vector<Field>& q) {
			
			setupDirectMatrix(q); // sets V_D(q)

			for (int i = 0; i < 19; ++i) for (int j = 0; j < 19; ++j) {
					V_Charge(i,j)     = -V_X_c(i,j) + 2.*V_D(i,j);
					V_Spin  (i,j)     = -V_X_s(i,j);
					
					V_Charge_coupl(i,j)      = -V_X_coupl(i,j) + 2.*V_D_coupl(i,j);
					V_Spin_coupl(i,j)        = -V_X_coupl(i,j);
			}
		}

		void calcRPAResult(ComplexMatrixType& chi0, 
						   size_t interactionType, 
						   ComplexMatrixType& fullInteraction, 
						   ComplexMatrixType& interactionMatrix, 
						   std::vector<Field> q=std::vector<Field>(3,0.0)) {

			// ComplexMatrixType interactionMatrix(19,19);
			setupDirectMatrix(q);	
			setupVBare(q);
			if (interactionType == 0) { // charge
				for (int i = 0; i < 19; ++i) for (int j = 0; j < 19; ++j) {
					// interactionMatrix(i,j) = -V_X(i,j) + 2.*V_D(i,j);
					interactionMatrix(i,j) = V_Charge(i,j);
				}
			} else if (interactionType==1) { // spin
				for (int i = 0; i < 19; ++i) for (int j = 0; j < 19; ++j) {
					// interactionMatrix(i,j) = -V_X(i,j);
					interactionMatrix(i,j) = V_Spin(i,j);
				}
			}


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
			matMul(interactionMatrix,chi0,c);
			for (size_t i = 0; i < c.n_row(); ++i) c(i,i) = 1.+c(i,i);
			// Now invert
			GETRF(n,n,c,n,ipiv,&info);
			if (info!=0) throw std::runtime_error("GETRF: failed\n");
			GETRI(n,c,n,ipiv,work,lwork,&info);
			if (info!=0) throw std::runtime_error("GETRI: failed\n");
			matMul(c,interactionMatrix,fullInteraction);



		}

		void calcChiRPAFromGammaRPA(ComplexMatrixType& Gamma, size_t interactionType, std::vector<Field>& q,
			                        ComplexMatrixType& ChiRPA) {
			int n = Gamma.n_row();
			// int m = matrix0.n_col();
			// int k = spinMatrix.n_col();
			ComplexType alpha(1.0);
			ComplexType beta(0.0);
			std::vector<int> ipiv(n);
			int info;
			ComplexMatrixType work(n,n);
			ComplexMatrixType c(n,n);
			int lwork(n);

			setupVBare(q); 
			ComplexMatrixType V_int(19,19);
			for (size_t i=0;i<19;i++) for (size_t j=0;j<19;j++) {
				if (interactionType==1) V_int(i,j) = V_Spin(i,j);
				else                    V_int(i,j) = V_Charge(i,j);
			}

			for (size_t i=0;i<19;i++) for (size_t j=0;j<19;j++) ChiRPA(i,j) = Gamma(i,j)- V_int(i,j);
			// Now invert V_int
			GETRF(n,n,V_int,n,ipiv,&info);
			if (info!=0) throw std::runtime_error("GETRF: failed\n");
			GETRI(n,V_int,n,ipiv,work,lwork,&info);
			if (info!=0) throw std::runtime_error("GETRI: failed\n");
			// Now multiply (Gamma-V_int) from both sides with 1/V_int
			matMul(ChiRPA,V_int,c);
			matMul(V_int,c,ChiRPA); // ChiRPA now contains the RPA chi matrix
			for (size_t i=0;i<19;i++) for (size_t j=0;j<19;j++) ChiRPA(i,j) *= -1;
		}

		void calcGammaFromChiRPA(ComplexMatrixType& V_bare, ComplexMatrixType& V_coupl, 
								 ComplexMatrixType& chiRPA, 
			                     ComplexMatrixType& Gamma) {
			// Gamma = V_Bare - V_Bare*Chi_RPA*V_Bare
			// setupVBare(q); 
			ComplexMatrixType c(19,19);
			matMul(V_coupl,chiRPA,c);
			matMul(c,V_coupl,Gamma);
			for (size_t i=0;i<19;i++) for (size_t j=0;j<19;j++) Gamma(i,j) = V_bare(i,j)-Gamma(i,j);
		}

	};

}

#endif
