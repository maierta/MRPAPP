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


namespace rpa {
	
	/*! \brief 			
	 */

	template<typename FieldType>
	class functorSwave {
	public:
		FieldType operator()(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b = 0.0, 
						     const FieldType& c = 0.0) const
		{
			return a 
				 + b * (cos(kvector[0]) + cos(kvector[1]))
				 + c * (cos(kvector[0]) * cos(kvector[1]));
			
		}
	};

	template<typename FieldType>
	class functorSwaveRot {
	public:
		FieldType operator()(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b, 
						     const FieldType& c ) const
		{
			std::vector<FieldType> kRot(3,0);
			kRot[0] = 0.5*(kvector[0]-kvector[1]);
			kRot[1] = 0.5*(kvector[0]+kvector[1]);
			return a 
				 + b * (cos(kRot[0]) + cos(kRot[1]))
				 + c * (cos(kRot[0]) * cos(kRot[1]));
			
		}
	};

	template<typename FieldType>
	class functorDwave {
	public:
		FieldType operator()(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b, 
						     const FieldType& c ) const
		{
			return a * (cos(kvector[0]) - cos(kvector[1])) 
				 + b * (cos(2.*kvector[0]) - cos(2.*kvector[1]))
				 + c * (cos(2.*kvector[0])*cos(kvector[1]) - cos(2.*kvector[1])*cos(kvector[0]));
			
		}
	};
	
	template<typename FieldType>
	class functorDwaveRot {
	public:
		FieldType operator()(const std::vector<FieldType>& kvector, 
						     const FieldType& a, 
						     const FieldType& b, 
						     const FieldType& c ) const
		{
			std::vector<FieldType> k(3,0);
			k[0] = 0.5*(kvector[0]-kvector[1]);
			k[1] = 0.5*(kvector[0]+kvector[1]);
			k[2] = kvector[2];

			return a * (cos(k[0]) - cos(k[1])) 
				 + b * (cos(2.*k[0]) - cos(2.*k[1]))
				 + c * (cos(2.*k[0])*cos(k[1]) - cos(2.*k[1])*cos(k[0]));
			
		}
	};
	
	
}

#endif
