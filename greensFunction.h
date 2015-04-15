//-*-C++-*-

#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class greensFunction:
		public MatrixTemplate<std::complex<Field> >
 	{
 	
	private:
		typedef std::complex<Field> ComplexType;
		typedef psimag::Matrix<ComplexType> ComplexMatrixType;
		typedef std::vector<Field> VectorType;
		const rpa::parameters<Field,MatrixTemplate>& param;
		size_t nOrb;
		size_t nktot;

 	public:
 		typedef psimag::Matrix<std::complex<Field> > BaseType;
		momentumDomain<Field,psimag::Matrix> momentumDomain1;
		size_t nwn;
 		BaseType& super;


	greensFunction(const rpa::parameters<Field,MatrixTemplate>& parameters):
		BaseType(size_t(pow(parameters.nkInt,parameters.dimension))*parameters.nwn,parameters.nOrb*parameters.nOrb),
		param(parameters),
		nOrb(param.nOrb),
		nktot(size_t(pow(param.nkInt,param.dimension))),
		momentumDomain1(param,param.nkInt,param.dimension),
		nwn(param.nwn),
		super (*this)
		{
			momentumDomain1.set_momenta();
			std::cout << "Starting initialization of G(k,iw) ... \n"; 
			init();
			std::cout << "G(k,iw) initialized... \n"; 
		}	

	void toZero() {
		size_t dim1(nktot*param.nwn);
		size_t dim2(nOrb*nOrb);
		for (size_t i=0; i<dim1; i++) {
			for (size_t j=0; j<dim2; j++) {
				(*this)(i,j) = ComplexType(0.0);
			}
		}
	}


	void init() {

		size_t nktot(momentumDomain1.nktot);
		VectorType ek(nOrb);
		ComplexMatrixType ak(nOrb,nOrb);
		VectorType k(3);
		Bands<Field,psimag::Matrix,ConcurrencyType> bands(param);


		for (size_t ik = 0; ik < nktot; ++ik)
		{
			momentumDomain1.momenta.getRow(ik,k);
			bands.calcBandsK(k,ek,ak);

			for (size_t in = 0; in < nwn; ++in)
			{
				Field iwn = (2*in+1)*param.pi_f*param.temperature;
				size_t ikw = in + ik * nwn;
				for (size_t l1 = 0; l1 < nOrb; ++l1)
				{
					for (size_t l2 = 0; l2 < nOrb; ++l2)
					{
						size_t ind = l2 + l1 * nOrb;
						ComplexType c2(0.0);	
						for (size_t iband = 0; iband < nOrb; ++iband)
						{
							ComplexType c1 = ak(l1,iband) * conj(ak(l2,iband));
							c2 += ComplexType(1.)/(ComplexType(-ek[iband]+param.mu,iwn)) * c1;						
						}
						(*this)(ikw,ind) = c2;
					}
				}
			}
		}
	}

	};
}

#endif

