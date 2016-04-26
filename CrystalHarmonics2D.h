//-*-C++-*-

/** \ingroup DCA */
/*@{*/

/*! \file CrystalHarmonics.h  
 * \author T.A.Maier
 *
 *  Contains a class to represent the crystal harmonics (s-wave, d-wave, etc.)
 */
#ifndef DCA_CrystalHarmonics_H
#define DCA_CrystalHarmonics_H

#include <string>
#include <vector>
#include <fstream>
#include <complex>

namespace rpa {
	
	/*! \brief 			
	 */

	template<typename FieldType>
	FieldType swave(const std::vector<FieldType>& kvector, 
				    const FieldType& a, 
					const FieldType& b, 
					const FieldType& c)
		{
			return a 
				 + b * (cos(kvector[0]) + cos(kvector[1]))
				 + c * (cos(kvector[0]) * cos(kvector[1]));
			
		}

	template<typename FieldType>
	FieldType swaveRot(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b, 
						     const FieldType& c )
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[0]-kvector[1]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			return a 
				 + b * (cos(kRot[0]) + cos(kRot[1]))
				 + c * (cos(kRot[0]) * cos(kRot[1]));
			
		}


	template<typename FieldType>
	FieldType dwave(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b, 
						     const FieldType& c )
		{
			return a * (cos(kvector[0]) - cos(kvector[1])) 
				 + b * (cos(2.*kvector[0]) - cos(2.*kvector[1]))
				 + c * (cos(2.*kvector[0])*cos(kvector[1]) - cos(2.*kvector[1])*cos(kvector[0]));
			
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType dwave(const std::vector<FieldType>& kvector, 
					const MatrixTemplate<FieldType>& w,
				    const std::vector<FieldType>& kz,
					const FieldType& k0,
					const size_t band)
			{
				return 0.5 * (cos(kvector[0]) - cos(kvector[1]));
					 // + w(1,0) * (cos(2.*kvector[0]) - cos(2.*kvector[1]))
					 // + w(2,0) * (cos(2.*kvector[0])*cos(kvector[1]) - cos(2.*kvector[1])*cos(kvector[0]));
			}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType swave(const std::vector<FieldType>& kvector, 
					const MatrixTemplate<FieldType>& w,
				    const std::vector<FieldType>& kz,
					const FieldType& k0,
					const size_t band)
			{
				return w(0,0) 
					 + w(1,0) * (cos(kvector[0]) + cos(kvector[1]))
					 + w(2,0) * (cos(kvector[0]) * cos(kvector[1]));
			}
	
	template<typename FieldType>
	FieldType dwaveRot(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b, 
						     const FieldType& c )
		{
			std::vector<FieldType> k(3,0);
			k[0] = 0.5*(kvector[0]-kvector[1]);
			k[1] = 0.5*(kvector[0]+kvector[1]);
			k[2] = kvector[2];

			return a * (cos(k[0]) - cos(k[1])) 
				 + b * (cos(2.*k[0]) - cos(2.*k[1]))
				 + c * (cos(2.*k[0])*cos(k[1]) - cos(2.*k[1])*cos(k[0]));
			
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType swaveRPALiFeAs(const std::vector<FieldType>& kvector, 
						     const  MatrixTemplate<FieldType>& w,
						     const  std::vector<FieldType>& kz,
						     const  FieldType& k0,
						     const size_t band)
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			kRot[2] = kvector[2];

			// Build amplitudes w_j(kz)
			std::vector<FieldType> wofkz(3,0);
			for (size_t j=0;j<3;j++) {
				wofkz[j] = 0.0;
				for (size_t i=0;i<6;i++) {
					wofkz[j] += w(i,j)*sqrt(pow(fabs(kRot[2]) - kz[i],2) + pow(k0,2));
				}
			}

			return wofkz[0] * (cos(kRot[0]) * cos(kRot[1])) +
			       // wofkz[1] * (cos(2.*kRot[0]) + cos(2.*kRot[1])) +
			       wofkz[1] * (cos(2.*kRot[0]) * cos(2.*kRot[1])) +
			       wofkz[2] * (cos(4.*kRot[0]) * cos(4.*kRot[1]));
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType swaveRPALiFeAs_2(const std::vector<FieldType>& kvector, 
						     	const MatrixTemplate<FieldType>& w,
						     	const std::vector<FieldType>& kz,
						        const FieldType& k0,
						     const size_t band)
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			kRot[2] = kvector[2];

			// Build amplitudes w_j(kz)
			std::vector<FieldType> wofkz(3,0);
			for (size_t j=0;j<4;j++) {
				wofkz[j] = 0.0;
				for (size_t i=0;i<7;i++) {
					wofkz[j] += w(i,j)*sqrt(pow(fabs(kRot[2]) - kz[i],2) + pow(k0,2));
				}
			}

			return wofkz[0] * (cos(kRot[0]) * cos(kRot[1])) +
			       wofkz[1] * (cos(2.*kRot[0]) * cos(2.*kRot[1])) +
			       wofkz[2] * (cos(2.*kRot[0]) + cos(2.*kRot[1])) +
			       wofkz[3] * (cos(4.*kRot[0]) * cos(4.*kRot[1]));
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType dwaveRPAKFe2Se2_1(const std::vector<FieldType>& kvector, 
						     	const MatrixTemplate<FieldType>& w,
						     	const std::vector<FieldType>& kz,
						     	const FieldType& k0,
						     	const size_t band)
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			// kRot[0] = kvector[0];
			// kRot[1] = kvector[1];
			kRot[2] = kvector[2];

			// Build amplitudes w_j(kz)
			std::vector<FieldType> wofkz(4,0);
			for (size_t j=0;j<w.n_col();j++) {
				wofkz[j] = 0.0;
				for (size_t i=0;i<w.n_row();i++) {
					wofkz[j] += w(i,j)*sqrt(pow(fabs(kRot[2]) - kz[i],2) + pow(k0,2));
				}
			}

			return wofkz[0] * (cos(kRot[0]) - cos(kRot[1])) +
			       wofkz[1] * (cos(2.*kRot[0]) - cos(2.*kRot[1])) +
			       wofkz[2] * (cos(3.*kRot[0]) - cos(3.*kRot[1])) +
			       wofkz[3] * (cos(4.*kRot[0]) - cos(4.*kRot[1]));
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType dwavePhenomKFe2Se2(const std::vector<FieldType>& kvector, 
						     	 const MatrixTemplate<FieldType>& w,
						     	 const std::vector<FieldType>& kz,
						     	 const FieldType& k0,
						     	const size_t band)
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			kRot[2] = kvector[2];

			FieldType Pi(3.141592653589793);

			FieldType rkx = floor(kRot[0]/Pi + 0.5);
			FieldType rky = floor(kRot[1]/Pi + 0.5);

			FieldType a1(3.027-2.5596*cos(kRot[2]*0.19242)); // note: hardcoded for FS for \mu=6.1897
			                                                 // Needs to be adjusted for other \mu !!

			FieldType f1s = cos(a1*(kRot[0]-Pi))-cos(kRot[1]);
			FieldType f1(1.0);
			if (f1s < 0.) f1=-1.0;

			FieldType f2s = cos(kRot[0])-cos(a1*(kRot[1]-Pi));
			FieldType f2(1.0);
			if (f2s < 0.) f2=-1.0;

			FieldType f = rkx*f1+rky*f2;

			if (band==6) {
				return f;
			} else if (band==7) {
				return -f;
			} else {
				return 0.0;
			}
		}

		template<typename FieldType, template<typename> class MatrixTemplate>
		FieldType dwavePhenomKFe2Se2_overdoped(const std::vector<FieldType>& kvector, 
									       	   const MatrixTemplate<std::complex<FieldType> >& ak,
											   size_t band
												   )
			{
				std::vector<FieldType> kRot(3,0);
				kRot[0] = 0.5*(kvector[1]-kvector[0]);
				kRot[1] = 0.5*(kvector[1]+kvector[0]);
				kRot[2] = kvector[2];

				FieldType Pi(3.141592653589793);
				FieldType Delta(0.);


				FieldType rkx = floor(kRot[0]/Pi + 0.5);
				FieldType rky = floor(kRot[1]/Pi + 0.5);

				FieldType reXZ1(real(ak(2,band)));
				FieldType imXZ1(imag(ak(2,band)));
				FieldType reXZ2(real(ak(7,band)));
				FieldType imXZ2(imag(ak(7,band)));
				FieldType reYZ1(real(ak(1,band)));
				FieldType imYZ1(imag(ak(1,band)));
				FieldType reYZ2(real(ak(6,band)));
				FieldType imYZ2(imag(ak(6,band)));
				FieldType reXY1(real(ak(3,band)));
				FieldType imXY1(imag(ak(3,band)));
				FieldType reXY2(real(ak(8,band)));
				FieldType imXY2(imag(ak(8,band)));
				// FieldType re551(real(ak(4,band)));
				// FieldType im551(imag(ak(4,band)));
				// FieldType re552(real(ak(9,band)));
				// FieldType im552(imag(ak(9,band)));

				FieldType wxz(0.5*(sqrt(pow(reXZ1,2)+pow(imXZ1,2))+sqrt(pow(reXZ2,2)+pow(imXZ2,2))));
				FieldType wyz(0.5*(sqrt(pow(reYZ1,2)+pow(imYZ1,2))+sqrt(pow(reYZ2,2)+pow(imYZ2,2))));
				FieldType wxy(0.5*(sqrt(pow(reXY1,2)+pow(imXY1,2))+sqrt(pow(reXY2,2)+pow(imXY2,2))));
				// FieldType w55(0.5*(sqrt(pow(re551,2)+pow(im551,2))+sqrt(pow(re552,2)+pow(im552,2))));

				if (rkx==1 && rky==0) { // k near (pi,0)
					if (cos(kRot[0])+cos(kRot[1]) > 0.0) { // inside 2 Fe BZ
						if      (wxy > wxz) Delta =  1.0;
						else                Delta = -1.0; 
					} else { // outside 2 Fe zone
						if      (wyz > wxy) Delta =  1.0;
						else                Delta = -1.0; 
					}
				} else if (rkx==0 && rky==1) { // k near (0,pi)
					if (cos(kRot[0])+cos(kRot[1]) > 0.0) { // inside 2 Fe BZ
						if      (wxy > wyz) Delta = -1.0;
						else                Delta =  1.0; 
					} else { // outside 2 Fe zone
						if      (wxz > wxy) Delta = -1.0;
						else                Delta =  1.0; 
					}
				}
				return Delta;
			}


	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType MazinKFe2Se2(const std::vector<FieldType>& kvector, 
						   const MatrixTemplate<FieldType>& w,
						   const std::vector<FieldType>& kz,
						   const FieldType& k0,
						   const size_t band)
		{

			return w(0,0);
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType MazinKFe2Se2_wZ(const std::vector<FieldType>& kvector, 
						   	  const MatrixTemplate<FieldType>& w,
						      const std::vector<FieldType>& kz,
						      const FieldType& k0,
						      const size_t band)
		{

			// std::vector<FieldType> kRot(3,0);
			// kRot[0] = 0.5*(kvector[1]-kvector[0]);
			// kRot[1] = 0.5*(kvector[0]+kvector[1]);
			// kRot[2] = kvector[2];

			// if (band==6 && cos(kRot[0])*cos(kRot[1]) > 0.0) { // we are dealing with the Z-pocket
			// 	return 0.5;
			// } else return w(0,0);
			return w(0,0);
		}

	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType dwaveRPAKFe2Se2_2(const std::vector<FieldType>& kvector, 
						     	const MatrixTemplate<FieldType>& w,
						     	const std::vector<FieldType>& kz,
						     	const FieldType& k0,
						     	const size_t band)
		{
			MatrixTemplate<FieldType> w_(w);
			std::vector<FieldType>   kz_(kz);
			FieldType                k0_(k0);
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			kRot[2] = kvector[2];

			if (band==6 && cos(kRot[0])*cos(kRot[1]) > 0.0) { // we are dealing with the Z-pocket
				//  Change kRot for kappa^00 Z-pocket
				FieldType Pi(3.141592653589793);
				if (fabs(kRot[0]) < Pi/2) { // shift for Z-pocket at Gamma
					if (kRot[2] < 0.0) kRot[2] += 4*Pi;
					kRot[0] += Pi;
					kRot[1] += Pi;
					kRot[2] -= 2.*Pi;
				}
				// Change w-Matrix (only first column required)
				// std::vector<FieldType> w0(6,0.0);
				w_(0,0) = -0.0324695593378; w_(0,1) = 0.0; w_(0,2) = 0.0; w_(0,3) = 0.0;  
				w_(1,0) = -0.0401810253714; w_(1,1) = 0.0; w_(1,2) = 0.0; w_(1,3) = 0.0;  
				w_(2,0) = -0.040565910487 ; w_(2,1) = 0.0; w_(2,2) = 0.0; w_(2,3) = 0.0;  
				w_(3,0) = -0.492403951557 ; w_(3,1) = 0.0; w_(3,2) = 0.0; w_(3,3) = 0.0;  
				w_(4,0) = 1.65805016676   ; w_(4,1) = 0.0; w_(4,2) = 0.0; w_(4,3) = 0.0;  	
				w_(5,0) = -1.0835968904   ; w_(5,1) = 0.0; w_(5,2) = 0.0; w_(5,3) = 0.0;  	
				// Change kz
                FieldType kzmax(2.33166717);
                FieldType nz(10.);
                // Change kz
                for (size_t i=0;i<6;i++) kz_[i] = float(i)*2.*kzmax/nz; 
                // Change k0
                k0_ = 2*kzmax*0.25/nz;
				
				// FieldType wofkz(0.);
				// for (size_t i=0;i<6;i++) wofkz += w(i,0)*sqrt(pow(fabs(kRot[2]) - kzi[i],2) + pow(k0nu,2));

				// return wofkz * (cos(kRot[0]) - cos(kRot[1]));

			} 
			// Build amplitudes w_j(kz)
			std::vector<FieldType> wofkz(4,0);
			for (size_t j=0;j<4;j++) {
				wofkz[j] = 0.0;
				for (size_t i=0;i<6;i++) {
					wofkz[j] += w_(i,j)*sqrt(pow(fabs(kRot[2]) - kz_[i],2) + pow(k0_,2));
				}
			} 

			return wofkz[0] * (cos(kRot[0]) - cos(kRot[1])) +
			       wofkz[1] * (cos(2.*kRot[0]) - cos(2.*kRot[1])) +
			       wofkz[2] * (cos(3.*kRot[0]) - cos(3.*kRot[1])) +
			       wofkz[3] * (cos(4.*kRot[0]) - cos(4.*kRot[1]));

		}


	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType swaveRPAKFe2Se2_elDoped(const std::vector<FieldType>& kvector, 
						     	      const MatrixTemplate<FieldType>& w,
						     	      const std::vector<FieldType>& kz,
						     		  const FieldType& k0,
						     	      const size_t band)
		{
			MatrixTemplate<FieldType> w_(w);
			std::vector<FieldType>   kz_(kz);
			FieldType                k0_(k0);

			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			kRot[2] = kvector[2];

			if (band==6 && cos(kRot[0])*cos(kRot[1]) > 0.0) { 
			    // we are dealing with the Z-pocket
				//  Change kRot for kappa^00 Z-pocket
				FieldType Pi(3.141592653589793);
				if (fabs(kRot[0]) < Pi/2) { // shift for Z-pocket at Gamma
					if (kRot[2] < 0.0) kRot[2] += 4*Pi;
					kRot[0] += Pi;
					kRot[1] += Pi;
					kRot[2] -= 2.*Pi;
				}
				// Change w-Matrix
				w_(0,0)=-0.760151271561; w_(0,1)=0.885918647045; w_(0,3)= -0.0300845885259;  w_(0,2)=0.0;
				w_(1,0)=0.803808939228;  w_(1,1)=-0.915304101901;w_(1,3)= -0.00599406603862; w_(1,2)=0.0;
				w_(2,0)=-0.753617453353; w_(2,1)=0.960455896346; w_(2,3)= -0.19825354557;    w_(2,2)=0.0;
				w_(3,0)=2.27593468314;   w_(3,1)=-3.41409614667; w_(3,3)=  1.16696394559;    w_(3,2)=0.0;
				w_(4,0)=-1.78152997008;  w_(4,1)=3.35682849276;  w_(4,3)= -1.70420046252;    w_(4,2)=0.0;
				w_(5,0)=0.334852407784;  w_(5,1)=-0.947721922153;w_(5,3)=  0.704058825776;   w_(5,2)=0.0;

				// Change kz
                FieldType kzmax(2.33166717);
                FieldType nz(10.);
                // std::vector<FieldType> kzi(6,0.0);
                for (size_t i=0;i<6;i++) kz_[i] = float(i)*2.*kzmax/nz; 
                // Change k0
                k0_ = 2*kzmax*0.25/nz;
				

			}
			// Build amplitudes w_j(kz)
			std::vector<FieldType> wofkz(4,0);
			for (size_t j=0;j<4;j++) {
				wofkz[j] = 0.0;
				for (size_t i=0;i<6;i++) {
					wofkz[j] += w_(i,j)*sqrt(pow(fabs(kRot[2]) - kz_[i],2) + pow(k0_,2));
				}
			}

			return wofkz[0] * (cos(kRot[0]) * cos(kRot[1])) +
			       wofkz[1] * (cos(2.*kRot[0]) * cos(2.*kRot[1])) +
			       wofkz[2] * (cos(3.*kRot[0]) * cos(3.*kRot[1])) +
			       wofkz[3] * (cos(4.*kRot[0]) * cos(4.*kRot[1]));
			
		}


	template<typename FieldType, template<typename> class MatrixTemplate>
	FieldType dwaveRPAKFe2Se2_3(const std::vector<FieldType>& kvector, 
						     	const MatrixTemplate<FieldType>& w,
						     	const std::vector<FieldType>& kz,
						     	const FieldType& k0,
						     	const size_t band)
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[1]-kvector[0]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			kRot[2] = kvector[2];

			if (band==6 && cos(kRot[0])*cos(kRot[1]) > 0.8) { // we are dealing with the Z-pocket
				//  Change kRot for kappa^00 Z-pocket
				FieldType Pi(3.141592653589793);
				if (fabs(kRot[0]) < Pi/2) { // shift for Z-pocket at Gamma
					if (kRot[2] < 0.0) kRot[2] += 4*Pi;
					kRot[0] += Pi;
					kRot[1] += Pi;
					kRot[2] -= 2.*Pi;
				}
				// Change w-Matrix (only first column required)
				std::vector<FieldType> w0(2,0.0);
				w0[0] = -0.235082868311;
				w0[1] = -0.0768135720907;
				// Change kz
                FieldType nz(10.);
                std::vector<FieldType> kzi(6,0.0);
                for (size_t i=0;i<2;i++) kzi[i] = float(i)*4.*Pi/nz; 
                // Change k0
                FieldType k0nu = Pi/nz;
				
				FieldType wofkz(0.);
				for (size_t i=0;i<w0.size();i++) wofkz += w0[i]*sqrt(pow(fabs(kRot[2]) - kzi[i],2) + pow(k0nu,2));

				// return wofkz * (cos(kRot[0]) - cos(kRot[1]));
				return wofkz * (cos(kvector[0]) - cos(kvector[1]));

			} else { // no Z-pocket
				// Build amplitudes w_j(kz)
				std::vector<FieldType> wofkz(4,0);
				for (size_t j=0;j<4;j++) {
					wofkz[j] = 0.0;
					for (size_t i=0;i<6;i++) {
						wofkz[j] += w(i,j)*sqrt(pow(fabs(kRot[2]) - kz[i],2) + pow(k0,2));
					}
				}

				// return wofkz[0] * (cos(kRot[0]) - cos(kRot[1])) +
				//        wofkz[1] * (cos(2.*kRot[0]) - cos(2.*kRot[1])) +
				//        wofkz[2] * (cos(3.*kRot[0]) - cos(3.*kRot[1])) +
				//        wofkz[3] * (cos(4.*kRot[0]) - cos(4.*kRot[1]));
				return wofkz[0] * (cos(kvector[0]) - cos(kvector[1])) +
				       wofkz[1] * (cos(2.*kvector[0]) - cos(2.*kvector[1])) +
				       wofkz[2] * (cos(3.*kvector[0]) - cos(3.*kvector[1])) +
				       wofkz[3] * (cos(4.*kvector[0]) - cos(4.*kvector[1]));
			}
		}



	
}

#endif
