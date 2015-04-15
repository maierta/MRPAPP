#ifndef SRRUO_H
#define SRRUO_H


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
	class SrRuO {
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


		SrRuO(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb)
		{
		}
		
		inline void getBands(const VectorType k, VectorType& eigenvals, ComplexMatrixType& eigenvects) {
		    FieldType tXYx,tXYy,tXZx,tXZy,tYZx,tYZy,tXYp,tXZp,tYZp,muXY,muXZ,muYZ,tperp;
		    tXYx = 0.44; tXZx = 0.31 ; tYZx = 0.045;
			tXYy = 0.44; tXZy = 0.045; tYZy = 0.31;
			tXYp =-0.14; tXZp = 0.01 ; tYZp = 0.01;
			muXY = 0.50; muXZ = 0.24 ; muYZ = 0.24;
			tperp = 0.1;
			FieldType cx,cy,cxy;
			cx = cos(k[0]); cy = cos(k[1]); cxy = cos(k[0])*cos(k[1]);
			FieldType ekXY    = -2*tXYx*cx-2*tXYy*cy+4*tXYp*cxy-muXY;
			FieldType ekXZ    = -2*tXZx*cx-2*tXZy*cy+4*tXZp*cxy-muXZ;
			FieldType ekYZ    = -2*tYZx*cx-2*tYZy*cy+4*tYZp*cxy-muYZ;
			FieldType ekPlus  = 0.5*(ekXZ+ekYZ);
			FieldType ekMinus = 0.5*(ekXZ-ekYZ);
			FieldType denom = sqrt(ekMinus*ekMinus+tperp*tperp);

			eigenvals[0] = ekPlus - denom;
			eigenvals[1] = ekPlus + denom;
			eigenvals[2] = ekXY;

			FieldType uk = sqrt(0.5*(1.0+ekMinus/denom));
			FieldType vk = sqrt(0.5*(1.0-ekMinus/denom));

			eigenvects(0,2) = 0.0; eigenvects(1,2) = 0.0; eigenvects(2,2) = 1.0;
			eigenvects(0,0) = -vk; eigenvects(1,0) = -uk; eigenvects(2,0) = 0.0;
			eigenvects(0,1) = -uk; eigenvects(1,1) =  vk; eigenvects(2,1) = 0.0;

			if (fabs(uk*uk+vk*vk-1.0)>=1.0e-5) std::cout << "Problem: EV not normalized!!" << uk*uk+vk*vk << " \n";

		}

	};
}

#endif
