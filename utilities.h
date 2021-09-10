//-*-C++-*-
#ifndef UTILITIES_H
#define UTILITIES_H

#include "Matrix.h"
#include <vector>



extern "C" void 
#ifdef glyph
	zheev_
#else
	zheev
#endif
	(char *,char *,int *,std::complex<double> *,int *,double *,std::complex<double> *,int *,double *,int *);
	
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
	dgetrf_
#else
	dgetrf
#endif
	(int *,int *,double *,int *,int *,int *);
	
	extern "C" void 
#ifdef glyph
	dgetri_
#else
	dgetri
#endif
	(int *, double *, int *, int *, double *, int *, int *);


inline void GETRF(int m, int n, 
 			      psimag::Matrix<double> &a,int lda, 
 			      std::vector<int> &ipiv, int *info) 
				 {
#ifdef glyph						
				psimag::LAPACK::dgetrf_
#else
				dgetrf
#endif
					(&m,&n,&(a(0,0)),&lda,&(ipiv[0]),info);
				}
inline void GETRI(int m, 
			      psimag::Matrix<double> &a, int lda, 
			      std::vector<int> &ipiv, 
				  psimag::Matrix<double> &work, int lwork, int *info) 
				 {
#ifdef glyph						
				psimag::LAPACK::dgetri_
#else
				dgetri
#endif
					(&m,&(a(0,0)),&lda,&(ipiv[0]),&(work(0,0)),&lwork,info);
				}

extern "C" void 
#ifdef glyph
	zgetrf_
#else
	zgetrf
#endif
	(int *,int *,std::complex<double> *,int *,int *,int *);
	
	extern "C" void 
#ifdef glyph
	zgetri_
#else
	zgetri
#endif
	(int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

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



typedef psimag::Matrix<std::complex<double> > ComplexMatrixType;
typedef psimag::Matrix<double>  MatrixType;


inline void matMul(ComplexMatrixType& matrix0, ComplexMatrixType& matrix1, ComplexMatrixType& matrix2) {

	int m = matrix0.n_row();
	int n = matrix1.n_col();
	int k = matrix0.n_col();
	std::complex<double> alpha(1.0);
	std::complex<double> beta(0.0);
	char transa('N');
	char transb('N');

#ifdef glyph
	psimag::BLAS::zgemm_
#else
	zgemm
#endif
		(&transa,&transb,&m,&n,&k,&alpha,&(matrix0(0,0)),
	  	 &m,&(matrix1(0,0)),&k,&beta,&(matrix2(0,0)),&m);
}

template<typename MType>
inline void calcInverse(MType& matrix) {

		int n = matrix.n_row();
		std::vector<int> ipiv(n);
		int info;
		MType work(n,n);
		int lwork(n);
		GETRF(n,n,matrix,n,ipiv,&info);
		if (info!=0) throw std::runtime_error("GETRF: failed\n");
		GETRI(n,matrix,n,ipiv,work,lwork,&info);
		if (info!=0) throw std::runtime_error("GETRI: failed\n");
}

template<typename FieldType>
std::ostream& operator<<(std::ostream& os,std::vector<FieldType>& v)
{
	for (size_t i=0;i<v.size()-1;i++) os <<v[i]<< " "; os << v[v.size()-1];
	return os;
}

// template<typename FieldType, template<typename> class MatrixTemplate>
// FieldType calcChiPhys(const rpa::parameters<FieldType,MatrixTemplate>& param, const MatrixType& chi) {
// 	FieldType chiPhys(0.0);
// 	// diagonal terms
// 	for (size_t l1 = 0; l1 < param.nOrb; ++l1)
// 	{
// 		for (size_t l2 = 0; l2 < param.nOrb; ++l2)
// 		{
// 			size_t ind1(l1+l1*param.nOrb);
// 			size_t ind2(l2+l2*param.nOrb);
// 			chiPhys += 0.5*real(chi(ind1,ind2)) ; 
// 		}
// 	}
// 	FieldType factor(1.0);
// 	if (param.sublattice==1) factor=2.0;

// 	return chiPhys/factor;
// }

// Needed for reading in CSV files (provided by G. Alvarez)
bool hasNonBlank(const std::string& b)
{
	for (size_t i=0;i<b.length();i++)
	if (b[i]!=' ' || b[i]!='\n' || b[i]!='\t') return true;
	return false;
}

template<typename FieldType>
void loadVector(std::vector<FieldType>& v,const std::string& myfile)
{
	std::ifstream fin(myfile.c_str());
	std::string buffer("");
	while(!fin.eof()) {
		char c;
		fin.get(c);
		if (c==',' || c=='\n') {
			if (hasNonBlank(buffer)) 
			v.push_back(atof(buffer.c_str())); // use atoi for ints
			buffer="";
		} else {
			buffer = buffer + c;
		}
	}
	fin.close();
}


template<typename FieldType, template<typename> class MatrixTemplate>
std::ostream& operator<<(std::ostream& os, std::vector<MatrixTemplate<std::complex<FieldType> > > & v)
{

	size_t nel(v.size());
	size_t nrow(v[0].n_row());
	size_t ncol(v[0].n_col());
	
	for (size_t i=0;i<nel;i++) {
		for (size_t j=0;j<nrow;j++) {
			for (size_t k=0;k<ncol;k++) {
				os << v[i](j,k);
			}
		}
		os << "\n";
	}
	return os;
}

template<typename FieldType, template<typename> class MatrixTemplate>
std::istream& operator>>(std::istream& is, std::vector<MatrixTemplate<std::complex<FieldType> > > & v)
{

	size_t nel(v.size());
	size_t nrow(v[0].n_row());
	size_t ncol(v[0].n_col());
	
	for (size_t i=0;i<nel;i++) {
		for (size_t j=0;j<nrow;j++) {
			for (size_t k=0;k<ncol;k++) {
				is >> v[i](j,k);
			}
		}
	}
	return is;
}


template<typename FieldType, template<typename> class MatrixTemplate>
void eigen(std::vector<FieldType>& eigenvals, MatrixTemplate<std::complex<FieldType> >& matrix) {
	int n = matrix.n_row();
	int lwork = 2*n-1;
	std::vector<std::complex<FieldType> > work(lwork);
	std::vector<FieldType> rwork(3*n-2);
	int info;

	GEEV('V','U',n,matrix,n,eigenvals,work,lwork,rwork,&info);
	if (info!=0) {
		throw std::runtime_error("zheev: failed\n");
	}
}

template<typename FieldType, template<typename> class MatrixTemplate>
void calcRPAResult(MatrixTemplate<std::complex<FieldType> >& matrix0, 
				   MatrixTemplate<std::complex<FieldType> >& interactionMatrix, 
				   MatrixTemplate<std::complex<FieldType> >& matrix1, 
				   std::vector<FieldType> q=std::vector<FieldType>(3,0.0)) {
	int n = interactionMatrix.n_row();
	// int m = matrix0.n_col();
	// int k = spinMatrix.n_col();
	// std::complex<FieldType> alpha(1.0);
	// std::complex<FieldType> beta(0.0);
	std::vector<int> ipiv(n);
	int info;
	MatrixTemplate<std::complex<FieldType> > work(n,n);
	MatrixTemplate<std::complex<FieldType> > c(n,n);
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

template<typename FieldType, typename SuscType, template<typename> class MatrixTemplate>
void calcRPAResult(SuscType& matrix0, 
		   MatrixTemplate<std::complex<FieldType> >& interactionMatrix, 
		   SuscType& matrix1, 
		   std::vector<FieldType> q=std::vector<FieldType>(3,0.0)) {
	int n = interactionMatrix.n_row();
	std::vector<int> ipiv(n);
	int info;
	MatrixTemplate<std::complex<FieldType> > work(n,n);
	MatrixTemplate<std::complex<FieldType> > c(n,n);
	int lwork(n);


	matMul(interactionMatrix,matrix0,c);
	// Result of matrix multiplication is in c
	for (size_t i = 0; i < c.n_row(); ++i) for (size_t j = 0; j < c.n_col(); ++j) c(i,j) = -c(i,j);
	for (size_t i = 0; i < c.n_row(); ++i) c(i,i) = 1.+c(i,i);
	// Now invert
	GETRF(n,n,c,n,ipiv,&info);
	if (info!=0) throw std::runtime_error("GETRF: failed\n");
	GETRI(n,c,n,ipiv,work,lwork,&info);
	if (info!=0) throw std::runtime_error("GETRI: failed\n");
	// Now multiply result with matrix0
	matMul(matrix0,c,matrix1);

}
#endif
