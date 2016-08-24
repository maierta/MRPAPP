#ifndef KFE2SE2_H
#define KFE2SE2_H


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
	class KFe2Se2 {
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


		KFe2Se2(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, ConcurrencyType& concurrency):
			param(parameters),
			conc(concurrency),
			dim(param.dimension),
			nbands(param.nOrb)
		{
			if (nbands!=5) {
				std::cout << "Number of bands should be 5 !!! \n";
				exit(0);
			}
		}
		
		inline void getBands(const VectorType k, ComplexMatrixType& eigenvects) {

			MatrixType tx(5,5);
			MatrixType ty(5,5);
			MatrixType txx(5,5);
			MatrixType tyy(5,5);
			MatrixType txy(5,5);
			MatrixType txxy(5,5);
			MatrixType txyy(5,5);
			MatrixType txxyy(5,5);
			MatrixType tz(5,5);
			MatrixType txz(5,5);
			MatrixType txxz(5,5);
			MatrixType txyz(5,5);
			MatrixType txxyz(5,5);

			size_t n,m;

			m=0;
		    tx(m,m)    = -0.00719007;
		    ty(m,m)    = -0.33281;
		    txx(m,m)   = -0.0571593;
		    tyy(m,m)   =  0.0231101;
		    txy(m,m)   = 0.206179;
		    txxy(m,m)  = -0.0257147;
		    txyy(m,m)  = -0.0242853;
		    txxyy(m,m) = 0.0304743;
		    tz(m,m)    = 0.0;
		    txz(m,m)   = 0.000150236;
		    txxz(m,m)  = -0.005;
		    txyz(m,m)  = 0.0;

		    m=2;
		
		    tx(m,m)    = 0.351144;
		    ty(m,m)    = 0.0;
		    txx(m,m)   = -0.00379318;
		    txy(m,m)   = -0.0384905;
		    txxy(m,m)  = 0.0;
		    txyy(m,m)  = 0.0;
		    txxyy(m,m) = 0.0;
		    tz(m,m)    = 0.0;
		    txz(m,m)   = 0.0;
		    txxz(m,m)  = 0.0;
		    txyz(m,m)  = 0.0;
		
		    m=3;

		    tx(m,m)    = 0.1;
		    ty(m,m)    = 0.0;
		    txx(m,m)   = -0.00727457;
		    txy(m,m)   = 0.0609889;
		    txxy(m,m)  = -0.0166046;
		    txyy(m,m)  = 0.0;
		    txxyy(m,m) = -0.00358076;
		    tz(m,m)    = 0.0552816;
		    txz(m,m)   = 0.0292682;
		    txxz(m,m)  = 0.0;
		    txyz(m,m)  = 0.0147694;
		
		    m=4;
		
		    tx(m,m)    = 0.03;
		    ty(m,m)    = 0.0;
		    txy(m,m)   = -0.0456289;
		    txx(m,m)   = 0.00625709;
		    txy(m,m)   = -0.0456289;
		    txxy(m,m)  = 0.00666156;
		    txyy(m,m)  = 0.0;
		    txxyy(m,m) = -0.0192934;
		    tz(m,m)    = -0.0165274;
		    txz(m,m)   = 0.002;
		    txxz(m,m)  = 0.0;
		    txyz(m,m)  = 0.0;

		    m=0;n=1;
		
		    tx(m,n)    = 0.0;
		    txy(m,n)   = 0.0941787;
		    txxy(m,n)  = -0.025;
		    txxyy(m,n) = 0.0304743;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.0;
		    txyz(m,n)  = 0.00135765;
		    txxyz(m,n) = 0.0;
		
		    m=0;n=2;
		
		    tx(m,n)    = -0.29371;
		    txy(m,n)   = 0.0574321;
		    txxy(m,n)  = 0.012815;
		    txxyy(m,n) = 0.0;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.0;
		    txyz(m,n)  = 0.0;
		    txxyz(m,n) = 0.0;
		
		
		    m=0;n=3;
		
		    tx(m,n)    = 0.240416;
		    txy(m,n)   = 0.0721671;
		    txxy(m,n)  = 0.0169706;
		    txxyy(m,n) = 0.0;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.00848528;
		    txyz(m,n)  = 0.0;
		    txxyz(m,n) = -0.00707107;
		
		
		    m=0;n=4;
		
		    tx(m,n)    = -0.063064;
		    txy(m,n)   = -0.104846;
		    txxy(m,n)  = 0.0;
		    txxyy(m,n) = 0.0;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.0;
		    txyz(m,n)  = -0.00370312;
		    txxyz(m,n) = 0.0;
		
		    m=1;n=3;
		
		    tx(m,n)    = 0.0;
		    txy(m,n)   = 0.0;
		    txxy(m,n)  = 0.0;
		    txxyy(m,n) = 0.0;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.0;
		    txyz(m,n)  = 0.0;
		    txxyz(m,n) = 0.0;
		
		    m=2;n=3;
		
		    tx(m,n)    = 0.0;
		    txy(m,n)   = 0.0;
		    txxy(m,n)  = -0.0204324;
		    txxyy(m,n) = 0.0;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.0;
		    txyz(m,n)  = 0.0;
		    txxyz(m,n) = 0.0;
		
		    m=2;n=4;
		
		    tx(m,n)    = -0.303532;
		    txy(m,n)   = 0.0;
		    txxy(m,n)  = -0.00497091;
		    txxyy(m,n) = 0.0;
		    tz(m,n)    = 0.0;
		    txz(m,n)   = 0.0;
		    txyz(m,n)  = 0.0;
		    txxyz(m,n) = 0.0;
		
		    m=3;n=4;
		
		    tx(m,n)    = 0.0;
		    txy(m,n)   = -0.0649461;
		    txxy(m,n)  = 0.0;
		    txxyy(m,n) = 0.01251;
		    tz(m,n)    = -0.0361762;
		    txz(m,n)   = 0.00841571;
		    txyz(m,n)  = 0.0;
		    txxyz(m,n) = 0.0;


			FieldType cx,cy,cz,c2x,c2y,sx,sy,sz,s2x,s2y;
			cx = cos(k[0]); 
			cy = cos(k[1]); 
			cz = cos(0);
			c2x = cos(2.*k[0]);  
			c2y = cos(2.*k[1]);  

			sx  = sin(k[0]); 
			sy  = sin(k[1]); 
			sz  = sin(0);
			s2x = sin(2.*k[0]);  
			s2y = sin(2.*k[1]);  

			const ComplexType ii = ComplexType(0.0,1.0);

		    eigenvects(0,0) = 2*tx(0,0)*cx + 2*ty(0,0)*cy                      
		           + 4*txy(0,0)*cx*cy + 2*txx(0,0)*c2x+2*tyy(0,0)*c2y      
		           + 4*txxy(0,0)*c2x*cy + 4*txyy(0,0)*c2y*cx      
		           + 4*txxyy(0,0)*c2x*c2y + 4*txz(0,0)*(cx+cy)*cz 
		           + 4*txxz(0,0)*(c2x-c2y)*cz 
		           - 0.293586;
		
		    eigenvects(1,1) = 2*ty(0,0)*cx + 2*tx(0,0)*cy                      
		           + 4*txy(0,0)*cx*cy + 2*txx(0,0)*c2y+ 2*tyy(0,0)*c2x      
		           + 4*txyy(0,0)*c2x*cy + 4*txxy(0,0)*c2y*cx      
		           + 4*txxyy(0,0)*c2x*c2y + 4*txz(0,0)*(cx+cy)*cz 
		           - 4*txxz(0,0)*(c2x-c2y)*cz 
		           - 0.293586;
		
		    eigenvects(2,2) = 2*tx(2,2)*(cx+cy) + 4*txy(2,2)*cx*cy + 2*txx(2,2)*(c2x+c2y) 
		              - 0.692525;
		
		    eigenvects(3,3) = 2*tx(3,3)*(cx+cy) + 4*txy(3,3)*cx*cy 
		           + 2*txx(3,3)*(c2x+c2y) + 4*txxy(3,3)*(c2x*cy+c2y*cx) 
		           + 4*txxyy(3,3)*c2x*c2y + 2*tz(3,3)*cz 
		           + 4*txz(3,3)*(cx+cy)*cz + 8*txyz(3,3)*cx*cy*cz 
		           - 0.340752;
		
		    eigenvects(4,4) = 2*tx(4,4)*(cx+cy) + 2*txy(4,4)*(cos(k[0]+k[1]) + cos(k[0]-k[1])) 
		           + 2*txx(4,4)*(c2x+c2y) 
		           + 4*txxy(4,4)*(c2x*cy+c2y*cx) 
		           + 4*txxyy(4,4)*c2x*c2y + 2*tz(4,4)*cz 
		           + 4*txz(4,4)*(cx+cy)*cz 
		           - 0.672753;
		
		    eigenvects(0,1) = 4*txy(0,1)*sx*sy + 4*txxy(0,1)*(s2x*sy+s2y*sx) 
		           + 4*txxyy(0,1)*s2x*s2y + 8*txyz(0,1)*sx*sy*cz;
		
		    eigenvects(0,2) = 2*ii*tx(0,2)*sy + 4*ii*txy(0,2)*sy*cx 
		           - 4*ii*txxy(0,2)*(s2y*cx-c2x*sy);
		
		    eigenvects(1,2) = 2*ii*tx(0,2)*sx + 4*ii*txy(0,2)*sx*cy 
		           - 4*ii*txxy(0,2)*(s2x*cy-c2y*sx);
		
		    eigenvects(0,3) = + 2*ii*tx(0,3)*sx + 4*ii*txy(0,3)*cy*sx 
		           + 4*ii*txxy(0,3)*s2x*cy + 4*ii*txz(0,3)*sx*cz 
		           - 4*txz(1,3)*sx*sz + 8*ii*txyz(0,3)*cy*sx*cz 
		           + 8*ii*txxyz(0,3)*s2x*cy*cz - 8*txxyz(2,4)*s2x*cy*sz;
		
		
		    eigenvects(1,3) = - 2*ii*tx(0,3)*sy - 4*ii*txy(0,3)*cx*sy 
		           - 4*ii*txxy(0,3)*s2y*cx - 4*ii*txz(0,3)*sy*cz 
		           - 4*txz(1,3)*sy*sz - 8*ii*txyz(0,3)*cx*sy*cz 
		           - 8*ii*txxyz(0,3)*s2y*cx*cz - 8*txxyz(1,3)*s2y*cx*sz;
		
		    eigenvects(0,4) = + 2*ii*tx(0,4)*sy - 4*ii*txy(0,4)*sy*cx 
		           - 8*ii*txyz(0,4)*sy*cx*cz;
		
		    eigenvects(1,4) = - 2*ii*tx(0,4)*sx + 4*ii*txy(0,4)*sx*cy 
		           + 8*ii*txyz(0,4)*sx*cy*cz;
		
		    eigenvects(2,3) = 4*txxy(2,3)*(s2y*sx-s2x*sy);
		
		    eigenvects(2,4) = 2*tx(2,4)*(cx-cy) + 4*txxy(2,4)*(c2x*cy-c2y*cx);
		
		    eigenvects(3,4) = 4*txy(3,4)*sx*sy + 4*txxyy(3,4)*s2x*s2y 
		           + 2*ii*tz(3,4)*sz + 4*ii*txz(3,4)*(cx+cy)*sz;

		}

	};
}

#endif
