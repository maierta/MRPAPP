//-*-C++-*-

#ifndef PSIMAG_VectorLike_H
#define PSIMAG_VectorLike_H

#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>


namespace psimag {
	
	template<typename T>
	class Vector : public std::vector<T> {
		
		template<typename ScalarType>
		Vector<T>& operator*=(Vector<T>& v,ScalarType x)
		{
			return VectorLike::times(v,x);
		}
	};

}; // psimag

#endif 

/*@}*/
#endif

