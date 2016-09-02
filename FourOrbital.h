#ifndef FOURORBITAL_H
#define FOURORBITAL_H


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
	class FourOrbital {
	private:
		typedef MatrixTemplate<Field> 		MatrixType;
		typedef std::complex<Field>			ComplexType;
		typedef MatrixTemplate<ComplexType> ComplexMatrixType;
		typedef std::vector<Field>      	VectorType;
		typedef Field 						FieldType;
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		size_t dim;
		
	public:
		FieldType nbands;


		FourOrbital(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb)
		{
			if (nbands!=4) {
				std::cout << "Number of orbitals should be 4 !!! \n";
				exit(0);
			}
		}
		
		inline void getBands(const VectorType k, ComplexMatrixType& eigenvects) {

			FieldType cx,cy,sx,sy;
			cx = cos(k[0]); 
			cy = cos(k[1]); 
			sx  = sin(k[0]); 
			sy  = sin(k[1]); 

			const FieldType t1(0.5);
			const FieldType t2(0.15);
			const FieldType t3(-0.175);
			const FieldType t4(-0.2);
			const FieldType t5(0.8);
			const FieldType t6(-0.450);
			const FieldType t7(-0.460);
			const FieldType t8(0.005);
			const FieldType t9(-0.8);
			const FieldType t10(-0.4);
			const FieldType t17(0.9);
			const FieldType Deltaxy(-0.6);
			const FieldType Deltax2y2(-2.0);


			const ComplexType ii = ComplexType(0.0,1.0);

		    eigenvects(0,0) = -2.*t2*cx-2.*t1*cy-4.*t3*cx*cy;
		    eigenvects(1,1) = -2.*t1*cx-2.*t2*cy-4.*t3*cx*cy;
		    eigenvects(0,1) = -4.*t4*sx*sy;
		    eigenvects(2,2) = -2.*t5*(-cx-cy)-4.*t6*cx*cy+Deltaxy;
		    eigenvects(0,2) = -4.*ii*t7*sx+8.*ii*t8*sx*cy;
		    eigenvects(1,2) = -4.*ii*t7*sy+8.*ii*t8*sy*cx;
		    eigenvects(3,3) = -2.*t17*(-cx-cy)-4.*t9*cx*cy+Deltax2y2;
		    eigenvects(0,3) = -4.*ii*t10*sy;
		    eigenvects(1,3) =  4.*ii*t10*sx;
		    eigenvects(2,3) =  0.0;

		}

	};
}

#endif
