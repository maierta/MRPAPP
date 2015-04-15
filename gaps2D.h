//-*-C++-*-

#ifndef GAPS2D_H
#define GAPS2D_H


#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib> // for atof and atoi

#include "Matrix.h"
#include "parameters.h"
#include "CrystalHarmonics2D.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, 
			 typename ConcurrencyType>
	class gap2D {    // simple s+- for 5-orbital 1111 model 2D

	private:
		typedef Field 				    FieldType;
		typedef std::complex<Field>		ComplexType;
		typedef std::vector<Field>      VectorType;
		typedef PsimagLite::Range<ConcurrencyType> RangeType;

		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t nbands;
		VectorType a,b,c;
		FieldType (*crystHarm)(const std::vector<FieldType>&, 
							   const FieldType&, 
						       const FieldType&, 
						       const FieldType&);

	public:

		gap2D() {}
				
		gap2D(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			  ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			nbands(param.nOrb),
			a(nbands,0),
			b(nbands,0),
			c(nbands,0)
		{			
			setParams2D();
		}

		ComplexType operator()(VectorType& k, size_t band) {
			FieldType Delta(0.0);
			Delta = crystHarm(k,a[band],b[band],c[band])*param.Delta0;
			return ComplexType(Delta,0.0);
		}
	
// This is the place where different gaps are hardcoded 
// in terms of amplitudes a,b,c for the s- or d-wave crystal harmonics

		void setParams2D() {
			if (param.gAmpl == "LaOFeAs_s_1") { // simple s+- for 5-orbital 1111 model 2D
				crystHarm = &swave;
				a[1] =  param.Delta0;
				a[2] =  param.Delta0;
				a[3] = -param.Delta0;
			} else if (param.gAmpl == "KFeSe_s_1") { // Mazin s+- for 10-orbital KFe2Se2 model 2D
				crystHarm = &swaveRot;
				a[7] =   1.0 * param.Delta0;
				a[6] =  -0.5 * param.Delta0;
				c[6] =  +0.7 * param.Delta0;
			} else if (param.gAmpl == "KFeSe_s_2") { // s++ for 10-orbital KFe2Se2 model 2D
				crystHarm = &swaveRot;
				a[7] =   1.0 * param.Delta0;
				a[6] =  0.65 * param.Delta0;
				c[6] =  -0.5 * param.Delta0;
			} else if (param.gAmpl == "KFeSe_xs_1") { // xs for 10-orbital KFe2Se2 model 2D
				crystHarm = &swaveRot;
				a[7] =   -0.3 * param.Delta0;
				b[7] =   -5.5 * param.Delta0;
				a[6] =  3.65 * param.Delta0; // constant
				b[6] =  4.8  * param.Delta0; // cos(kx)+cos(ky)
				c[6] =  5.37 * param.Delta0; // cos(kx)*cos(ky)
			} else if (param.gAmpl == "KFeSe_d_1") { // d-wave for 10-orbital KFe2Se2 model 2D
				crystHarm = &dwaveRot;
				a[7] = -0.58 * param.Delta0;
				a[6] =  0.60 * param.Delta0;
			}  else if (param.gAmpl == "KFeSe_d_2") { // nodal d-wave for hybridized 10-orbital KFe2Se2 model 2D
				crystHarm = &dwaveRot;
				a[7] = -0.9 * param.Delta0;
				b[7] =  0.6 * param.Delta0;
				c[7] = -1.2 * param.Delta0;
				a[6] =  0.988 * param.Delta0;
				b[6] = -0.78 * param.Delta0;
				c[6] =  2.08 * param.Delta0;
			} else if (param.gAmpl == "dwave") {
                a[0] = 1.0;
                b[0] = 0.0;
                c[0] = 0.0;
                crystHarm = &dwave;

            }else {
				std::cout << "Gap not implemented! Bailing out.\n";
				exit(0);
			}
		}

	};

}


#endif
