//-*-C++-*-

#ifndef SEPBASIS_H
#define SEPBASIS_H

#include <string>
#include <vector>
#include <fstream>
#include "parameters.h"


namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class sepBasis {

		private:
			typedef std::complex<Field>				ComplexType;
			typedef MatrixTemplate<Field> 			MatrixType;
			typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
			typedef std::vector<Field> 				VectorType;
			const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
			ConcurrencyType& conc;
			VectorType& k;

		public:
			MatrixType gMatrix;


		sepBasis(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, 
				 ConcurrencyType& concurrency,
				 VectorType& kIn):
			param(parameters),
			conc(concurrency),
			k(kIn),
			gMatrix(19,9)
		{
			setupGMatrix();
		}

		Field operator()(const size_t iB,const size_t l1,const size_t l2) {
			size_t ind(l2+l1*param.nOrb);
			return gMatrix(iB,ind);
		}

		void setupGMatrix() {
			for (int iB=0;iB<19;iB++) for (int l1=0;l1<3;l1++) for (int l2=0;l2<3;l2++){
				size_t ind(l2+l1*param.nOrb);
				gMatrix(iB,ind) = getGFor(iB,l1,l2);
			}
		}

		Field getGFor(const size_t iB,const size_t l1,const size_t l2) {
			// indexing for Emery model: orbital 0 = d; 1 = px; 2 = py

			switch (iB)
			{
				case 0:
					if (l1==0 && l2==1) return cos(0.5*k[0]); else return 0.0;
					break;
				case 1:
					if (l1==0 && l2==1) return sin(0.5*k[0]); else return 0.0;
					break;
				case 2:
					if (l1==0 && l2==2) return cos(0.5*k[1]); else return 0.0;
					break;
				case 3:
					if (l1==0 && l2==2) return sin(0.5*k[1]); else return 0.0;
					break;
				case 4:
					if (l1==1 && l2==2) return cos(0.5*k[0])*cos(0.5*k[1]); else return 0.0;
					break;
				case 5:
					if (l1==1 && l2==2) return cos(0.5*k[0])*sin(0.5*k[1]); else return 0.0;
					break;
				case 6:
					if (l1==1 && l2==2) return sin(0.5*k[0])*cos(0.5*k[1]); else return 0.0;
					break;
				case 7:
					if (l1==1 && l2==2) return sin(0.5*k[0])*sin(0.5*k[1]); else return 0.0;
					break;
				case 8:
					if (l1==0 && l2==0) return 1.0; else return 0.0;
					break;
				case 9:
					if (l1==1 && l2==1) return 1.0; else return 0.0;
					break;
				case 10:
					if (l1==2 && l2==2) return 1.0; else return 0.0;
					break;
				case 11:
					if (l1==1 && l2==0) return cos(0.5*k[0]); else return 0.0;
					break;
				case 12:
					if (l1==1 && l2==0) return sin(0.5*k[0]); else return 0.0;
					break;
				case 13:
					if (l1==2 && l2==0) return cos(0.5*k[1]); else return 0.0;
					break;
				case 14:
					if (l1==2 && l2==0) return sin(0.5*k[1]); else return 0.0;
					break;
				case 15:
					if (l1==2 && l2==1) return cos(0.5*k[0])*cos(0.5*k[1]); else return 0.0;
					break;
				case 16:
					if (l1==2 && l2==1) return cos(0.5*k[0])*sin(0.5*k[1]); else return 0.0;
					break;
				case 17:
					if (l1==2 && l2==1) return sin(0.5*k[0])*cos(0.5*k[1]); else return 0.0;
					break;
				case 18:
					if (l1==2 && l2==1) return sin(0.5*k[0])*sin(0.5*k[1]); else return 0.0;
					break;
				default:
					return 0;
			}
		}
	};
}
#endif
