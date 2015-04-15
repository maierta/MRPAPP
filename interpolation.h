//-*-C++-*-

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <string>
#include <vector>
#include <fstream>
// #include "math.h"
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"

namespace rpa{

	template<typename Field,template<typename> class MatrixTemplate,
			 typename SuscType, typename ConcurrencyType>

	class interpolation{

	private:
		typedef Field 						    FieldType;
		typedef std::complex<Field>				ComplexType;
		typedef MatrixTemplate<Field> 			MatrixType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      		VectorType;
		
		const momentumDomain<Field,MatrixTemplate,ConcurrencyType>&	kMesh;
		const std::vector<SuscType>& function;

	public:

		interpolation(const momentumDomain<FieldType,MatrixTemplate,ConcurrencyType>& kMeshIn,
					  const std::vector<SuscType>& functionIn):
			kMesh(kMeshIn),
			function(functionIn)
			{

			}

		void BiLinear(VectorType& q,SuscType& result) const {
			size_t ind00(0),ind10(0),ind01(0),ind11(0);
			kMesh.getSurroundingIndices(q,ind00,ind10,ind01,ind11);
			VectorType dis(2,0);

			FieldType G0Length = sqrt(pow(kMesh.b(0,0),2) + pow(kMesh.b(0,1),2));
			FieldType G1Length = sqrt(pow(kMesh.b(1,0),2) + pow(kMesh.b(1,1),2));
			FieldType dG0 = G0Length/kMesh.nk;
			FieldType dG1 = G1Length/kMesh.nk;

		    // dis[0] = (q[0]-kMesh.momenta(ind00,0))/(kMesh.momenta(ind11,0)-kMesh.momenta(ind00,0));
		    // dis[1] = (q[1]-kMesh.momenta(ind00,1))/(kMesh.momenta(ind11,1)-kMesh.momenta(ind00,1));
		 //    if ((kMesh.momenta(ind11,0)-kMesh.momenta(ind00,0))<0.0) {
		 //    	std::cout << "q: " << q << "\n";
			//     std::cout << "dG0: " << kMesh.momenta(ind11,0)-kMesh.momenta(ind00,0) << "\n";
			//     std::cout << "q11: " << kMesh.momenta(ind11,0) << "," <<kMesh.momenta(ind11,1)<< "\n";
			//     std::cout << "q00: " << kMesh.momenta(ind00,0) << "," <<kMesh.momenta(ind00,1) << "\n";
			// }
		    dis[0] = (q[0]-kMesh.momenta(ind00,0))/dG0;
		    dis[1] = (q[1]-kMesh.momenta(ind00,1))/dG1;
			// std::cout << ind00<<","<<ind10<<","<<ind01<<","<<ind11<<"dis: " << dis << "\n";
			// std::cout << "q,k:" << q << " , (" << kMesh.momenta(ind00,0) 
			// 	<< "," << kMesh.momenta(ind00,1) << ") , (" << kMesh.momenta(ind11,0) << "," << kMesh.momenta(ind11,1) << ")\n";
		    size_t msize(function[0].n_row());
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				ComplexType c00(0.0); ComplexType c10(0.0); ComplexType  c0(0.0);
				c00 = function[ind00](l1,l2) * (1.-dis[0]) + function[ind10](l1,l2) * dis[0];
				c10 = function[ind01](l1,l2) * (1.-dis[0]) + function[ind11](l1,l2) * dis[0];
				result(l1,l2) = c00 * (1.-dis[1]) + c10 * dis[1];
			}
		}

		void BiLinearGeneral(VectorType& q,SuscType& result) const {
			size_t ind00(0),ind10(0),ind01(0),ind11(0);
			kMesh.getSurroundingIndices(q,ind00,ind10,ind01,ind11);
			VectorType dis(2,0);

			FieldType G0Length = sqrt(pow(kMesh.b(0,0),2) + pow(kMesh.b(0,1),2));
			FieldType G1Length = sqrt(pow(kMesh.b(1,0),2) + pow(kMesh.b(1,1),2));
			FieldType dG0 = G0Length/kMesh.nk;
			FieldType dG1 = G1Length/kMesh.nk;

			FieldType dqx(q[0]-kMesh.momenta(ind00,0));
			FieldType dqy(q[1]-kMesh.momenta(ind00,1));

			dis[0] = (dqx*kMesh.b(0,0) + dqy*kMesh.b(0,1)) / G0Length / dG0;
			dis[1] = (dqx*kMesh.b(1,0) + dqy*kMesh.b(1,1)) / G1Length / dG1;

			// std::cout << ind00<<","<<ind10<<","<<ind01<<","<<ind11<<"dis: " << dis << "\n";
			// std::cout << "q,k:" << q << " , (" << kMesh.momenta(ind00,0) 
				// << "," << kMesh.momenta(ind00,1) << ") , (" << kMesh.momenta(ind11,0) << "," << kMesh.momenta(ind11,1) << ")\n";
		    size_t msize(function[0].n_row());
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				ComplexType c00(0.0); ComplexType c10(0.0); ComplexType  c0(0.0);
				c00 = function[ind00](l1,l2) * (1.-dis[0]) + function[ind10](l1,l2) * dis[0];
				c10 = function[ind01](l1,l2) * (1.-dis[0]) + function[ind11](l1,l2) * dis[0];
				result(l1,l2) = c00 * (1.-dis[1]) + c10 * dis[1];
			}
		}

		void TriLinear(VectorType& q,SuscType& result) const {
			size_t ind000(0),ind001(0),ind010(0),ind100(0),ind011(0),ind101(0),ind110(0),ind111(0);
			kMesh.getSurroundingIndices(q,ind000,ind001,ind010,ind100,ind011,ind101,ind110,ind111);
			VectorType dis(3,0);
			for(size_t l=0;l<3;l++) dis[l] = (q[l]-kMesh.momenta(ind000,l)) / 
							                 (kMesh.momenta(ind111,l)-kMesh.momenta(ind000,l));

		    size_t msize(function[0].n_row());
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				ComplexType c00(0.0); ComplexType c10(0.0); ComplexType  c01(0.0);
				ComplexType c11(0.0); ComplexType  c0(0.0); ComplexType   c1(0.0);
				c00 = function[ind000](l1,l2)*(1.-dis[0])+function[ind100](l1,l2)*dis[0];
				c10 = function[ind010](l1,l2)*(1.-dis[0])+function[ind110](l1,l2)*dis[0];
				c01 = function[ind001](l1,l2)*(1.-dis[0])+function[ind101](l1,l2)*dis[0];
				c11 = function[ind011](l1,l2)*(1.-dis[0])+function[ind111](l1,l2)*dis[0];
				c0  = c00 * (1.-dis[1]) + c10 * dis[1];
				c1  = c01 * (1.-dis[1]) + c11 * dis[1];

				result(l1,l2) = c0 * (1.-dis[2])+c1 * dis[2];

			}
		}


		void TriLinearGeneral(VectorType& q,SuscType& result) const {
			size_t ind000(0),ind001(0),ind010(0),ind100(0),ind011(0),ind101(0),ind110(0),ind111(0);
			kMesh.getSurroundingIndices(q,ind000,ind001,ind010,ind100,ind011,ind101,ind110,ind111);
			VectorType dis(3,0);

			FieldType G0Length = sqrt(pow(kMesh.b(0,0),2) + pow(kMesh.b(0,1),2) + pow(kMesh.b(0,2),2));
			FieldType G1Length = sqrt(pow(kMesh.b(1,0),2) + pow(kMesh.b(1,1),2) + pow(kMesh.b(1,2),2));
			FieldType G2Length = sqrt(pow(kMesh.b(2,0),2) + pow(kMesh.b(2,1),2) + pow(kMesh.b(2,2),2));
			FieldType dG0 = G0Length/kMesh.nk;
			FieldType dG1 = G1Length/kMesh.nk;
			FieldType dG2 = G2Length/kMesh.nkz;

			FieldType dqx(q[0]-kMesh.momenta(ind000,0));
			FieldType dqy(q[1]-kMesh.momenta(ind000,1));
			FieldType dqz(q[2]-kMesh.momenta(ind000,2));

			dis[0] = (dqx*kMesh.b(0,0) + dqy*kMesh.b(0,1) + dqz*kMesh.b(0,2)) / G0Length / dG0;
			dis[1] = (dqx*kMesh.b(1,0) + dqy*kMesh.b(1,1) + dqz*kMesh.b(1,2)) / G1Length / dG1;
			dis[2] = (dqx*kMesh.b(2,0) + dqy*kMesh.b(2,1) + dqz*kMesh.b(2,2)) / G2Length / dG2;

		    size_t msize(function[0].n_row());
			for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
				ComplexType c00(0.0); ComplexType c10(0.0); ComplexType  c01(0.0);
				ComplexType c11(0.0); ComplexType  c0(0.0); ComplexType   c1(0.0);
				c00 = function[ind000](l1,l2)*(1.-dis[0])+function[ind100](l1,l2)*dis[0];
				c10 = function[ind010](l1,l2)*(1.-dis[0])+function[ind110](l1,l2)*dis[0];
				c01 = function[ind001](l1,l2)*(1.-dis[0])+function[ind101](l1,l2)*dis[0];
				c11 = function[ind011](l1,l2)*(1.-dis[0])+function[ind111](l1,l2)*dis[0];
				c0  = c00 * (1.-dis[1]) + c10 * dis[1];
				c1  = c01 * (1.-dis[1]) + c11 * dis[1];

				result(l1,l2) = c0 * (1.-dis[2])+c1 * dis[2];

			}
		}



	};


}


#endif
