#ifndef BAFEAS_5ORB_H
#define BAFEAS_5ORB_H


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
	class BaFeAs {
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


		BaFeAs(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb)
		{
		}
		
		inline void getBands(const VectorType k, ComplexMatrixType& eigenvects) {
		    
		    MatrixType tx(5,5),ty(5,5),txx(5,5),txy(5,5),txxy(5,5),txyy(5,5),
		    	       txxyy(5,5),tz(5,5),txz(5,5),txxz(5,5),txyz(5,5),txxyz(5,5);

		    size_t m,n;

		    m=0;
		    tx(m,m)		=	-0.0604;
		    ty(m,m)		=	-0.3005;
		    txx(m,m)	=	 0.0253;
		    txy(m,m)	=	 0.2388;
		    txxy(m,m)	=	-0.0414;
		    txyy(m,m)	=	-0.0237;
		    txxyy(m,m)	=	 0.0158;
		    tz(m,m)		=	 0.0;
		    txz(m,m)	=	-0.0101;
		    txxz(m,m)	=	 0.0126;
		    txyz(m,m)	=	 0.0;

		    m=2;
		    tx(m,m)		=	 0.3378;
		    ty(m,m)		=	 0.0;
		    txx(m,m)	=	 0.0011;
		    txy(m,m)	=	-0.0947;
		    txxy(m,m)	=	 0.0;
		    txyy(m,m)	=	 0.0;
		    txxyy(m,m)	=	 0.0;
		    tz(m,m)		=	 0.0;
		    txz(m,m)	=	 0.0;
		    txxz(m,m)	=	 0.0;
		    txyz(m,m)	=	 0.0;

		    m=3;
		    tx(m,m)		=	 0.1965;
		    ty(m,m)		=	 0.0;
		    txx(m,m)	=	-0.0528;
		    txy(m,m)	=	 0.1259;
		    txxy(m,m)	=	-0.032;
		    txyy(m,m)	=	 0.0;
		    txxyy(m,m)	=	 0.0045;
		    tz(m,m)		=	 0.1001;
		    txz(m,m)	=	 0.0662;
		    txxz(m,m)	=	 0.0;
		    txyz(m,m)	=	 0.0421;

		    m=4;
		    tx(m,m)		=	-0.0656;
		    ty(m,m)		=	 0.0;
		    txx(m,m)	=	 0.0001;
		    txy(m,m)	=	 0.0;
		    txxy(m,m)	=	 0.01;
		    txyy(m,m)	=	 0.0;
		    txxyy(m,m)	=	 0.0047;
		    tz(m,m)		=	 0.0563;
		    txz(m,m)	=	-0.0036;
		    txxz(m,m)	=	 0.0;
		    txyz(m,m)	=	 0.0;

		    m=0;n=1;
		    tx(m,n)		=	0.0;
		    txy(m,n)	=	0.1934;
		    txxy(m,n)	=  -0.0325;
		    txxyy(m,n)	=   0.0158;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0;
		    txyz(m,n)	=  -0.0168;
		    txxyz(m,n)	=   0.0;

		    m=0;n=2;
		    tx(m,n)		=  -0.4224;
		    txy(m,n)	=	0.0589;
		    txxy(m,n)	=   0.0005;
		    txxyy(m,n)	=   0.0;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0;
		    txyz(m,n)	=   0.0;
		    txxyz(m,n)	=   0.0;

		    m=0;n=3;
		    tx(m,n)		=   0.1549;
		    txy(m,n)	=  -0.007;
		    txxy(m,n)	=  -0.0055;
		    txxyy(m,n)	=   0.0;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0524;
		    txyz(m,n)	=   0.0349;
		    txxyz(m,n)	=   0.0018;

		    m=0;n=4;
		    tx(m,n)		=  -0.0526;
		    txy(m,n)	=  -0.0862;
		    txxy(m,n)	=   0.0;
		    txxyy(m,n)	=   0.0;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0;
		    txyz(m,n)	=  -0.0203;
		    txxyz(m,n)	=   0.0;

		    m=1;n=3;
		    tx(m,n)		=   0.0;
		    txy(m,n)	=   0.0;
		    txxy(m,n)	=   0.0;
		    txxyy(m,n)	=   0.0;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0566;
		    txyz(m,n)	=   0.0;
		    txxyz(m,n)	=   0.0283;

		    m=2;n=3;
		    tx(m,n)		=   0.0;
		    txy(m,n)	=   0.0;
		    txxy(m,n)	=  -0.0108;
		    txxyy(m,n)	=   0.0;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0;
		    txyz(m,n)	=   0.0;
		    txxyz(m,n)	=   0.0;

		    m=2;n=4;
		    tx(m,n)		=  -0.2845;
		    txy(m,n)	=   0.0;
		    txxy(m,n)	=   0.0046;
		    txxyy(m,n)	=   0.0;
		    tz(m,n)		=   0.0;
		    txz(m,n)	=   0.0;
		    txyz(m,n)	=   0.0;
		    txxyz(m,n)	=   0.0;

		    m=3;n=4;
		    tx(m,n)		=   0.0;
		    txy(m,n)	=  -0.0475;
		    txxy(m,n)	=   0.0;
		    txxyy(m,n)	=   0.0004;
		    tz(m,n)		=  -0.019;
		    txz(m,n)	=  -0.0023;
		    txyz(m,n)	=   0.0;
		    txxyz(m,n)	=   0.0;


			FieldType cx,cy,c2x,c2y,sx,sy,s2x,s2y,cz,sz;
			const ComplexType ii = ComplexType(0.0,1.0);
			cx = cos(k[0]); cy = cos(k[1]); cz = cos(k[2]);
			sx = sin(k[0]); sy = sin(k[1]); sz = sin(k[2]);
			c2x = cos(2.0*k[0]); c2y = cos(2.0*k[1]);
			s2x = sin(2.0*k[0]); s2y = sin(2.0*k[1]);

			eigenvects(0,0) = 2.0*tx(0,0)*cx + 2.0*ty(0,0)*cy +
			                  4.0*txy(0,0)*cx*cy + 2.0*txx(0,0)*(c2x-c2y) +
			                  4.0*txxy(0,0)*c2x*cy + 4.0*txyy(0,0)*c2y*cx + 
			                  4.0*txxyy(0,0)*c2x*c2y + 4.0*txz(0,0)*(cx+cy)*cz + 
			                  4.0*txxz(0,0)*(c2x-c2y)*cz
			                  + 0.0987;

			eigenvects(1,1) = 2.0*ty(0,0)*cx + 2.0*tx(0,0)*cy +
			                  4.0*txy(0,0)*cx*cy - 2.0*txx(0,0)*(c2x-c2y) +
			                  4.0*txyy(0,0)*c2x*cy + 4.0*txxy(0,0)*c2y*cx + 
			                  4.0*txxyy(0,0)*c2x*c2y + 4.0*txz(0,0)*(cx+cy)*cz -
			                  4.0*txxz(0,0)*(c2x-c2y)*cz
			                  + 0.0987;
			
			eigenvects(2,2) = 2.0*tx(2,2)*(cx+cy) + 4.0*txy(2,2)*cx*cy +
							  2.0*txx(2,2)*(c2x+c2y)
							  - 0.3595;

			eigenvects(3,3) = 2.0*tx(3,3)*(cx+cy) + 4.0*txy(3,3)*cx*cy + 
							  2.0*txx(3,3)*(c2x+c2y) + 4.0*txxy(3,3)*(c2x*cy+c2y*cx) + 
							  4.0*txxyy(3,3)*c2x*c2y + 2.0*tz(3,3)*cz +
							  4.0*txz(3,3)*(cx+cy)*cz + 8.0*txyz(3,3)*cx*cy*cz
							  + 0.2078;

			eigenvects(4,4) = 2.0*tx(4,4)*(cx+cy) + 2.0*txx(4,4)*(c2x+c2y) + 
							  4.0*txxy(4,4)*(c2x*cy+c2y*cx) + 
							  4.0*txxyy(4,4)*c2x*c2y + 2.0*tz(4,4)*cz +
							  4.0*txz(4,4)*(cx+cy)*cz
							  - 0.7516;

			eigenvects(0,1) = 4.0*txy(0,1)*sx*sy + 4.0*txxy(0,1)*(s2x*sy+s2y*sx) + 
							  4.0*txxyy(0,1)*s2x*s2y + 8.0*txyz(0,1)*sx*sy*cz;

			eigenvects(0,2) = 2.0*ii*tx(0,2)*sy + 4.0*ii*txy(0,2)*sy*cx - 
							  4.0*ii*txxy(0,2)*(s2y*cx-c2x*sy);

			eigenvects(1,2) = 2.0*ii*tx(0,2)*sx + 4.0*ii*txy(0,2)*sx*cy - 
							  4.0*ii*txxy(0,2)*(s2x*cy-c2y*sx);

			eigenvects(0,3) = 2.0*ii*tx(0,3)*sx + 4.0*ii*txy(0,3)*cy*sx +
			                  4.0*ii*txxy(0,3)*s2x*cy + 4.0*ii*txz(0,3)*sx*cz -
			                  4.0*txz(1,3)*sx*sz + 8.0*ii*txyz(0,3)*cy*sx*cz +
			                  8.0*ii*txxyz(0,3)*s2x*cy*cz -
			                  8.0*txxyz(1,3)*s2x*cy*sz;

			eigenvects(1,3) = -2.0*ii*tx(0,3)*sy - 4.0*ii*txy(0,3)*cx*sy -
			                   4.0*ii*txxy(0,3)*s2y*cx - 4.0*ii*txz(0,3)*sy*cz -
			                   4.0*txz(1,3)*sy*sz - 8.0*ii*txyz(0,3)*cx*sy*cz -
			                   8.0*ii*txxyz(0,3)*s2y*cx*cz -
			                   8.0*txxyz(1,3)*s2y*cx*sz;

			eigenvects(0,4) = +2.0*ii*tx(0,4)*sy - 4.0*ii*txy(0,4)*sy*cx -
							   8.0*ii*txyz(0,4)*sy*cx*cz;

			eigenvects(1,4) = -2.0*ii*tx(0,4)*sx + 4.0*ii*txy(0,4)*sx*cy +
							   8.0*ii*txyz(0,4)*sx*cy*cz;

			eigenvects(2,3) = 4.0*txxy(2,3)*(s2y*sx-s2x*sy);

			eigenvects(2,4) = 2.0*tx(2,4)*(cx-cy) + 4.0*txxy(2,4)*(c2x*cy-c2y*cx);

			eigenvects(3,4) = 4.0*txy(3,4)*sx*sy + 4.0*txxyy(3,4)*s2x*s2y + 2.0*ii*tz(3,4)*sz +
							  4.0*ii*txz(3,4)*(cx+cy)*sz;

		}

	};
}

#endif
