//-*-C++-*-

#ifndef GAP_H
#define GAP_H


#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib> // for atof and atoi

#include "Matrix.h"
#include "parameters.h"
#include "CrystalHarmonics.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, 
			 typename ConcurrencyType>
	class gap {

	private:
		typedef MatrixTemplate<Field> 	MatrixType;
		typedef std::complex<Field>		ComplexType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      VectorType;
		typedef Field 					FieldType;
		typedef PsimagLite::Range<ConcurrencyType> RangeType;

		const rpa::parameters<Field,MatrixTemplate>& param;
		ConcurrencyType& conc;
		size_t nbands;
		VectorType a,b,c;

	public:


		gap(const rpa::parameters<Field,MatrixTemplate>& parameters,
			ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			nbands(param.nOrb),
			a(nbands,0),
			b(nbands,0),
			c(nbands,0)
		{			
			setAmplitudes();
		}

		ComplexType operator()(VectorType& k, size_t band) {
			functorSwave<FieldType>    swave;
			functorSwaveRot<FieldType> swaveRot;
			functorDwave<FieldType>    dwave;
			functorDwaveRot<FieldType> dwaveRot;

			FieldType Delta(0.0);

			if (param.gapSym=="s")         Delta = swave(k,a[band],b[band],c[band]);
			else if (param.gapSym=="d")    Delta = dwave(k,a[band],b[band],c[band]);
			else if (param.gapSym=="sRot") Delta = swaveRot(k,a[band],b[band],c[band]);
			else if (param.gapSym=="dRot") Delta = dwaveRot(k,a[band],b[band],c[band]);
			return ComplexType(Delta,0.0);
		}


// This is the place where different gaps are hardcoded 
// in terms of amplitudes a,b,c for the s- or d-wave crystal harmonics

		void setAmplitudes() {
			if (param.gAmpl == "LaOFeAs_s_1") { // simple s+- for 5-orbital 1111 model 2D
				a[1] =  param.Delta0;
				a[2] =  param.Delta0;
				a[3] = -param.Delta0;
			} else if (param.gAmpl == "KFeSe_s_1") { // Mazin s+- for 10-orbital KFe2Se2 model 2D
				a[7] =   1.0 * param.Delta0;
				a[6] =  -0.5 * param.Delta0;
				c[6] =  +0.7 * param.Delta0;
			} else if (param.gAmpl == "KFeSe_s_2") { // s++ for 10-orbital KFe2Se2 model 2D
				a[7] =   1.0 * param.Delta0;
				a[6] =  0.65 * param.Delta0;
				c[6] =  -0.5 * param.Delta0;
			} else if (param.gAmpl == "KFeSe_xs_1") { // xs for 10-orbital KFe2Se2 model 2D
				a[7] =   -0.3 * param.Delta0;
				b[7] =   -5.5 * param.Delta0;
				a[6] =  3.65 * param.Delta0; // constant
				b[6] =  4.8  * param.Delta0; // cos(kx)+cos(ky)
				c[6] =  5.37 * param.Delta0; // cos(kx)*cos(ky)
			} else if (param.gAmpl == "KFeSe_d_1") { // d-wave for 10-orbital KFe2Se2 model 2D
				a[7] = -0.58 * param.Delta0;
				a[6] =  0.60 * param.Delta0;
			} else {
				std::cout << "Gap not implemented! Bailing out.\n";
				exit(0);
			}
		}

	};

}


#endif
