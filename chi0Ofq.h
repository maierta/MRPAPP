//-*-C++-*-

#ifndef CHI0OFQ_H
#define CHI0OFQ_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "greensFunction.h"
#include "utilities.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class chi0ofq
	{

	private:
		typedef std::complex<Field> ComplexType;
		typedef psimag::Matrix<ComplexType> ComplexMatrixType;
		typedef std::vector<Field> VectorType;
		const rpa::parameters<Field,MatrixTemplate>& param;
		const rpa::greensFunction<Field,MatrixTemplate,ConcurrencyType>& g;
		ConcurrencyType& conc;
		size_t nOrb;
		size_t msize;
		size_t nktot;

	public:
		typedef std::vector<psimag::Matrix<std::complex<Field> > > ChiqMatrixType;
		const momentumDomain<Field,psimag::Matrix>& momentumDomain1;
		size_t nwn;
		// ChiqMatrixType chi0matrix;
		ComplexMatrixType chi0matrix;
		std::vector<Field> chiq;


		chi0ofq(const rpa::parameters<Field,MatrixTemplate>& parameters, const rpa::greensFunction<Field,MatrixTemplate,ConcurrencyType>& green, ConcurrencyType& concurrency):
		param(parameters),
		g(green),
		conc(concurrency),
		nOrb(param.nOrb),
		msize(nOrb*nOrb),
		nktot(size_t(pow(param.nkInt,param.dimension))),
		momentumDomain1(green.momentumDomain1),
		nwn(param.nwn),
		chi0matrix(msize,msize),
		chiq(nktot,0)
		{
			// for (size_t iq = 0; iq < nktot; ++iq) {
			// 	chi0matrix[iq] = ComplexMatrixType(msize,msize);
			// 	for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) chi0matrix[iq](i,j) = ComplexType(0.0,0.0);
			// }

			typedef PsimagLite::Range<ConcurrencyType> RangeType;
			RangeType range(0,nktot,conc);

	
			for (;!range.end();range.next()) {
				size_t iq = range.index();
				if (conc.rank()==0) std::cout << "iq=" << iq << " of " << nktot << " total\n";
				ComplexMatrixType matrix(msize,msize);
				calcMatrix(iq,matrix);
				chiq[iq] = calcChiPhys(param,matrix);
			}
			conc.reduce(chiq);
		}	

		void calcMatrix(const size_t iq, ComplexMatrixType& matrix) {
			// Calculate chi0(q) using Green's function matrix

			for (size_t ik = 0; ik < nktot; ++ik)
			{
				size_t ikq = momentumDomain1.indexOfAdd(iq,ik);

				for (size_t iw = 0; iw < nwn; ++iw)
				{
					size_t ind1 = iw + ik  * nwn;
					size_t ind2 = iw + ikq * nwn;

					for (size_t l1 = 0; l1 < nOrb; ++l1)
					{
						for (size_t l2 = 0; l2 < nOrb; ++l2)
						{
							for (size_t l3 = 0; l3 < nOrb; ++l3)
							{
								for (size_t l4 = 0; l4 < nOrb; ++l4)
								{
									size_t iorb1 = l2 + l1 * nOrb;
									size_t iorb2 = l4 + l3 * nOrb;
									
									// chi0matrix[iq](iorb1,iorb2) += g(ind1,iorb1) * g(ind2,iorb2);
									matrix(iorb1,iorb2) += g(ind1,iorb1) * g(ind2,iorb2);
								}
							}
						}
					}
				}
			}

			// chi0matrix[iq] *= -2.0*param.temperature/nktot;
			matrix *= -2.0*param.temperature/nktot;
		}


		};
	}

	#endif


