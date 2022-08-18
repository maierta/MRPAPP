//-*-C++-*-

#ifndef RPA_MOMENTUMDOMAIN_H
#define RPA_MOMENTUMDOMAIN_H

#include "parameters.h"
#include "Matrix.h"
#include "utilities.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class momentumDomain{

	private:
		typedef Field 								FieldType;
		typedef std::complex<Field>					ComplexType;
		typedef MatrixTemplate<Field>			 	MatrixType;
		typedef MatrixTemplate<std::complex<Field> >	ComplexMatrixType;
		typedef MatrixTemplate<size_t> 				MatrixIntType;
		typedef std::vector<Field>      			VectorType;
		typedef std::vector<std::complex<Field> >   ComplexVectorType;
		
		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>&   param;
		ConcurrencyType& conc;
		VectorType a1,a2,a3;
		VectorType shift;
		MatrixType dGInverse;

	public:
		size_t dim;
		size_t nk;
		size_t nkz;
		size_t nktot;
		MatrixType b;
		MatrixType momenta;
		MatrixIntType indexOfAdd;
		MatrixIntType index2Components;
		// VectorType dG;
		// VectorType GNorm;

		size_t index();

		momentumDomain(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
					   ConcurrencyType& concurrency, 
					   const size_t& numberOfkpoints,
					   const size_t& numberOfkzpoints,
					   const size_t& dimension):
			param(parameters),
			conc(concurrency),
			a1(param.a1),
			a2(param.a2),
			a3(param.a3),
			shift(3,-0.5),
			// shift(3,0.0),
			dGInverse(3,3),
			dim(dimension),
			nk(numberOfkpoints),
			nkz((dim==3)?numberOfkzpoints:1),
			nktot((dim==3)?size_t(nk*nk*nkz):size_t(nk*nk)),
			b(3,3),
			momenta(nktot,3),
			// indexOfAdd(nktot,nktot),
			index2Components(nktot,3)
			// dG(3,0),
			// GNorm(3,0)
			{
				set_primitiveVectors();
				if (conc.rank()==0) std::cout << "nktot: " << nktot << "\n";

				// shift[0] = -0.25; shift[1] = -0.25; shift[2] = -0.5;
			}

		momentumDomain(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, 
					   ConcurrencyType& concurrency, 
					   const size_t& numberOfkpoints,
					   const size_t& numberOfkzpoints,
					   const size_t& dimension,
					   const VectorType& p1, 
					   const VectorType& p2, 
					   const VectorType& p3):
			param(parameters),
			conc(concurrency),
			a1(p1),
			a2(p2),
			a3(p3),
			shift(3,-0.5),
			dGInverse(3,3),
			dim(dimension),
			nk(numberOfkpoints),
			nkz((dim==3)?numberOfkzpoints:1),
			nktot((dim==3)?size_t(nk*nk*nkz):size_t(nk*nk)),
			b(3,3),
			momenta(nktot,3),
			indexOfAdd(nktot,nktot),
			index2Components(nktot,3)
			// dG(3,0),
			// GNorm(3,0)
			{
				set_primitiveVectors();
			}


		momentumDomain(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, 
					   ConcurrencyType& concurrency, 
					   const size_t& numberOfkpoints,
					   const size_t& dimension):
			param(parameters),
			conc(concurrency),
			a1(param.a1),
			a2(param.a2),
			a3(param.a3),
			shift(3,-0.5),
			dGInverse(3,3),
			dim(dimension),
			nk(numberOfkpoints),
			nkz(1),
			nktot(nk*nk),
			b(3,3),
			momenta(nktot,3), 
			// indexOfAdd(nktot,nktot),
			index2Components(nktot,3)
			// dG(3,0),
			// GNorm(3,0)
			{
				set_primitiveVectors();
			}

		momentumDomain(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
					   ConcurrencyType& concurrency, 
					   const std::string& path, 
					   size_t numberOfkpoints): 
			param(parameters),
			conc(concurrency),
			a1(param.a1),
			a2(param.a2),
			a3(param.a3),
			shift(3,-0.5),
			dGInverse(3,3),
			nk(numberOfkpoints),
			nktot(nk),
			b(3,3),
			momenta(nk,3)
			{ 
				set_primitiveVectors();

				if (path == "Path1") set_momenta_Path1(); 
				if (path == "Path2") set_momenta_Path2(); 
				if (path == "Path3") set_momenta_Path3(); 
				if (path == "triangular") set_momenta_Path4(); 
				if (path == "triangular2") set_momenta_Path4_b(); 
				if (path == "triangularExtended") set_momenta_Path5(); 
			}



		void set_primitiveVectors() {
			if (conc.rank()==0) {
				std::cout << "Primitive vectors of direct lattice:" << "\n";
				std::cout << "a1: (" << a1[0]<<","<<a1[1]<<","<<a1[2]<<")\n";
				std::cout << "a2: (" << a2[0]<<","<<a2[1]<<","<<a2[2]<<")\n";
				std::cout << "a3: (" << a3[0]<<","<<a3[1]<<","<<a3[2]<<")\n";
			}
			Field volume = a1[0]*(a2[1]*a3[2]-a2[2]*a3[1]) +
						   a1[1]*(a2[2]*a3[0]-a2[0]*a3[2]) +
						   a1[2]*(a2[0]*a3[1]-a2[1]*a3[0]);
			b(0,0) = 2.*param.pi_f/volume * (a2[1]*a3[2]-a2[2]*a3[1]);
			b(0,1) = 2.*param.pi_f/volume * (a2[2]*a3[0]-a2[0]*a3[2]);
			b(0,2) = 2.*param.pi_f/volume * (a2[0]*a3[1]-a2[1]*a3[0]);

			b(1,0) = 2.*param.pi_f/volume * (a3[1]*a1[2]-a3[2]*a1[1]);
			b(1,1) = 2.*param.pi_f/volume * (a3[2]*a1[0]-a3[0]*a1[2]);
			b(1,2) = 2.*param.pi_f/volume * (a3[0]*a1[1]-a3[1]*a1[0]);

			b(2,0) = 2.*param.pi_f/volume * (a1[1]*a2[2]-a1[2]*a2[1]);
			b(2,1) = 2.*param.pi_f/volume * (a1[2]*a2[0]-a1[0]*a2[2]);
			b(2,2) = 2.*param.pi_f/volume * (a1[0]*a2[1]-a1[1]*a2[0]);

			if (conc.rank()==0) {
				std::cout << "Primitive vectors of reciprocal lattice:" << "\n";
				std::cout << "b:" << "\n";
				std::cout << "b1: (" << b(0,0)<<","<<b(0,1)<<","<<b(0,2)<<")\n";
				std::cout << "b2: (" << b(1,0)<<","<<b(1,1)<<","<<b(1,2)<<")\n";
				std::cout << "b3: (" << b(2,0)<<","<<b(2,1)<<","<<b(2,2)<<")\n";
			}

			// Now setting up matrix dGInverse which is needed to map k-point onto its index
			std::vector<FieldType> nki(3);
			nki[0] = nk; nki[1] = nk; nki[2] = nkz;
			for (size_t l1=0; l1<3; l1++) for (size_t l2=0; l2<3; l2++) 
				dGInverse(l1,l2) = b(l2,l1) / std::max(nki[l2],FieldType(1));
			// Now invert
			calcInverse(dGInverse);
		}


		void set_momenta_FS(const VectorType& shift) {
			if (dim==2) {	
				for (size_t ikx = 0; ikx < nk; ++ikx)
				{
					for (size_t iky = 0; iky < nk; ++iky)
					{
						momenta(ikx+iky*nk,0) = float(ikx)*2.*param.pi_f/float(nk) - shift[0];
						momenta(ikx+iky*nk,1) = float(iky)*2.*param.pi_f/float(nk) - shift[1];
						momenta(ikx+iky*nk,2) = 0.0;
					}
				}
			} else if (dim==3) {
				for (size_t ikx = 0; ikx < nk; ++ikx)
				{
					for (size_t iky = 0; iky < nk; ++iky)
					{
						for (size_t ikz = 0; ikz < nkz; ++ikz)
						{
							momenta(ikz+iky*nkz+ikx*nkz*nk,0) = float(ikx)*2.*param.pi_f/float(nk-1) - shift[0];
							momenta(ikz+iky*nkz+ikx*nkz*nk,1) = float(iky)*2.*param.pi_f/float(nk-1) - shift[1];
							momenta(ikz+iky*nkz+ikx*nkz*nk,2) = float(ikz)*2.*param.pi_f/float(nkz-1)- shift[2];
						}
					}
				}

			}
		}

		void set_momenta(const bool& indexation) {
						 	
					if (dim==2) {	
						if (conc.rank()==0) std::cout << "shift=" << shift[0] << ", b[0,0]" << b(0,0) << ", b[1,0]" << b(1,0) << "\n";
						for (size_t ikx = 0; ikx < nk; ++ikx) {
							for (size_t iky = 0; iky < nk; ++iky) {
								size_t ind = index(ikx,iky);
								for (size_t i=0; i<3; i++) {
									momenta(ind,i) = (float(ikx)/float(nk) + shift[0]) * b(0,i) + 
													 (float(iky)/float(nk) + shift[1]) * b(1,i) +
													  						param.kz2D * b(2,i);
								}
								// std::cout << "momenta=" << momenta(ind,0) << "," << momenta(ind,1) << "\n";

								index2Components(ind,0) = ikx;
								index2Components(ind,1) = iky;
								index2Components(ind,2) = 0;
							}
						}
					} else if (dim==3) {
						if (conc.rank()==0) std::cout << "shift=" << shift[0] << "," << shift[1] << "," << shift[2] << "\n";
						for (size_t ikz = 0; ikz < nkz; ++ikz) {
							for (size_t iky = 0; iky < nk; ++iky) {
								for (size_t ikx = 0; ikx < nk; ++ikx) {
									size_t ind = index(ikx,iky,ikz);
									for (size_t i=0; i<3; i++) {
										momenta(ind,i) = (float(ikx)/float(nk)  + shift[0]) * (*this).b(0,i) + 
														 (float(iky)/float(nk)  + shift[1]) * (*this).b(1,i) +
													     (float(ikz)/float(nkz) + shift[2]) * (*this).b(2,i) ;
									}
									
									index2Components(ind,0) = ikx;
									index2Components(ind,1) = iky;
									index2Components(ind,2) = ikz;
								}
							}
						}
					}
			if (indexation) set_indexOfAdd();	

				}

		void set_momenta(const Field& kxmin, const Field& kxmax, 
						const Field& kymin, const Field& kymax,
						const Field& kzmin, const Field& kzmax) {
			// Assumes a tetragonal BZ
			for (size_t ikz = 0; ikz < nkz; ++ikz) {
				for (size_t iky = 0; iky < nk; ++iky) {
					for (size_t ikx = 0; ikx < nk; ++ikx) {
						size_t ind = index(ikx,iky,ikz);
						index2Components(ind,0) = ikx;
						index2Components(ind,1) = iky;
						index2Components(ind,2) = ikz;

						momenta(ind,0) = kxmin + float(ikx)/float(nk)  * (kxmax - kxmin);
						momenta(ind,1) = kymin + float(iky)/float(nk)  * (kymax - kymin);
						momenta(ind,2) = kzmin + float(ikz)/float(nkz) * (kzmax - kzmin);
					}
				}
			}

		}

		void set_momenta(const Field& kxmin, const Field& kxmax, 
						const Field& kymin, const Field& kymax) {
			// Assumes a tetragonal BZ
			for (size_t ikx = 0; ikx < nk; ++ikx) {
				for (size_t iky = 0; iky < nk; ++iky) {
						size_t ind = index(ikx,iky);
					index2Components(ind,0) = ikx;
					index2Components(ind,1) = iky;
					index2Components(ind,2) = 0;

					momenta(ind,0) = kxmin + float(ikx)/float(nk-1) * (kxmax - kxmin);
					momenta(ind,1) = kymin + float(iky)/float(nk-1) * (kymax - kymin);
					momenta(ind,2) = 0;
				}
			}

		}

		void set_momenta_Path1() {

			size_t nks(nk/4);
			size_t ind(0);
			for (size_t ik=0; ik<nks; ik++) { // Gamma -> M
				momenta(ind,0) = float(ik)/float(nks) * param.pi_f;
				momenta(ind,1) = float(ik)/float(nks) * param.pi_f;
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { 
				momenta(ind,0) = param.pi_f + float(ik)/float(nks) * param.pi_f;
				momenta(ind,1) = param.pi_f - float(ik)/float(nks) * param.pi_f;
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) {
				momenta(ind,0) = 2.*param.pi_f - float(ik)/float(nks) * param.pi_f;
				momenta(ind,1) = 0.0;
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) {
				momenta(ind,0) = param.pi_f - float(ik)/float(nks) * param.pi_f;
				momenta(ind,1) = 0.0;
				momenta(ind,2) = 0.0;
				ind += 1;
			}

		}

		void set_momenta_Path2() { // Gamma -> X -> M -> Gamma

			size_t nks(nk/3);
			size_t ind(0);
			for (size_t ik=0; ik<nks; ik++) { // Gamma -> X // 0 --> 0.5*b1
				momenta(ind,0) = float(ik)/float(nks) * 0.5*(*this).b(0,0);
				momenta(ind,1) = float(ik)/float(nks) * 0.5*(*this).b(0,1);
				momenta(ind,2) = float(ik)/float(nks) * 0.5*(*this).b(0,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // X -> M // 0.5*b1 --> 0.5*b1 + 0.5*b2
				momenta(ind,0) = 0.5*(*this).b(0,0) + float(ik)/float(nks) * 0.5*(*this).b(1,0);
				momenta(ind,1) = 0.5*(*this).b(0,1) + float(ik)/float(nks) * 0.5*(*this).b(1,1);
				momenta(ind,2) = 0.5*(*this).b(0,2) + float(ik)/float(nks) * 0.5*(*this).b(1,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // M -> Gamma // 0.5*b1+0.5*b2 --> 0
				momenta(ind,0) = 0.5*(*this).b(0,0) + 0.5*(*this).b(1,0) - float(ik)/float(nks) * (0.5*(*this).b(0,0) + 0.5*(*this).b(1,0));
				momenta(ind,1) = 0.5*(*this).b(0,1) + 0.5*(*this).b(1,1) - float(ik)/float(nks) * (0.5*(*this).b(0,1) + 0.5*(*this).b(1,1));
				momenta(ind,2) = 0.5*(*this).b(0,2) + 0.5*(*this).b(1,2) - float(ik)/float(nks) * (0.5*(*this).b(0,2) + 0.5*(*this).b(1,2));
				ind += 1;
			}
		}

		void set_momenta_Path3() { // Gamma -> X -> M -> Gamma -> Z -> R -> A -> Z

			size_t nks(nk/7);
			size_t ind(0);
			for (size_t ik=0; ik<nks; ik++) { // Gamma -> X // 0 --> 0.5*b1
				momenta(ind,0) = float(ik)/float(nks) * 0.5*(*this).b(0,0);
				momenta(ind,1) = float(ik)/float(nks) * 0.5*(*this).b(0,1);
				momenta(ind,2) = float(ik)/float(nks) * 0.5*(*this).b(0,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // X -> M // 0.5*b1 --> 0.5*b1 + 0.5*b2
				momenta(ind,0) = 0.5*(*this).b(0,0) + float(ik)/float(nks) * 0.5*(*this).b(1,0);
				momenta(ind,1) = 0.5*(*this).b(0,1) + float(ik)/float(nks) * 0.5*(*this).b(1,1);
				momenta(ind,2) = 0.5*(*this).b(0,2) + float(ik)/float(nks) * 0.5*(*this).b(1,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // M -> Gamma // 0.5*b1+0.5*b2 --> 0
				momenta(ind,0) = 0.5*(*this).b(0,0) + 0.5*(*this).b(1,0) - float(ik)/float(nks) * (0.5*(*this).b(0,0) + 0.5*(*this).b(1,0));
				momenta(ind,1) = 0.5*(*this).b(0,1) + 0.5*(*this).b(1,1) - float(ik)/float(nks) * (0.5*(*this).b(0,1) + 0.5*(*this).b(1,1));
				momenta(ind,2) = 0.5*(*this).b(0,2) + 0.5*(*this).b(1,2) - float(ik)/float(nks) * (0.5*(*this).b(0,2) + 0.5*(*this).b(1,2));
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // Gamma -> Z // 0 --> 0.5*b3
				momenta(ind,0) = float(ik)/float(nks) * 0.5*(*this).b(2,0); 
				momenta(ind,1) = float(ik)/float(nks) * 0.5*(*this).b(2,1); 
				momenta(ind,2) = float(ik)/float(nks) * 0.5*(*this).b(2,2); 
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // Z -> R // 0.5*b3 --> 0.5*b1+0.5*b3
				momenta(ind,0) = 0.5*(*this).b(2,0) + float(ik)/float(nks) * 0.5*(*this).b(0,0);
				momenta(ind,1) = 0.5*(*this).b(2,1) + float(ik)/float(nks) * 0.5*(*this).b(0,1);
				momenta(ind,2) = 0.5*(*this).b(2,2) + float(ik)/float(nks) * 0.5*(*this).b(0,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // R -> A // 0.5*b1+0.5*b3 --> 0.5*b1+0.5*b2+0.5*b3
				momenta(ind,0) = 0.5*(*this).b(0,0) + 0.5*(*this).b(2,0) + float(ik)/float(nks) * 0.5*(*this).b(1,0);
				momenta(ind,1) = 0.5*(*this).b(0,1) + 0.5*(*this).b(2,1) + float(ik)/float(nks) * 0.5*(*this).b(1,1);
				momenta(ind,2) = 0.5*(*this).b(0,2) + 0.5*(*this).b(2,2) + float(ik)/float(nks) * 0.5*(*this).b(1,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // A -> Z
				momenta(ind,0) = 0.5*(*this).b(0,0) + 0.5*(*this).b(2,0) + 0.5*(*this).b(1,0) - float(ik)/float(nks) * (0.5*(*this).b(0,0) + 0.5*(*this).b(2,0) + 0.5*(*this).b(1,0)); 
				momenta(ind,1) = 0.5*(*this).b(0,1) + 0.5*(*this).b(2,1) + 0.5*(*this).b(1,1) - float(ik)/float(nks) * (0.5*(*this).b(0,1) + 0.5*(*this).b(2,1) + 0.5*(*this).b(1,1));
				momenta(ind,2) = 0.5*(*this).b(0,2) + 0.5*(*this).b(2,2) + 0.5*(*this).b(1,2) - float(ik)/float(nks) * (0.5*(*this).b(0,2) + 0.5*(*this).b(2,2) + 0.5*(*this).b(1,2));
				ind += 1;
			}
		}

		void set_momenta_Path4() { //Gamma -> M -> K -> Gamma for triangular lattice

			size_t nks(nk/3);
			size_t ind(0);
			for (size_t ik=0; ik<nks; ik++) { // Gamma = 0 -> M = (pi, -pi/sqrt(3))
				momenta(ind,0) = float(ik)/float(nks) * param.pi_f;
				momenta(ind,1) = float(ik)/float(nks) * (-param.pi_f/sqrt(3.));
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // M -> K = (4*pi/3, 0)
				momenta(ind,0) = param.pi_f + float(ik)/float(nks) * param.pi_f/3.;
				momenta(ind,1) = -param.pi_f/sqrt(3.) + float(ik)/float(nks) * param.pi_f/sqrt(3.);
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // K -> Gamma = (0,0)
				momenta(ind,0) = 4.*param.pi_f/3. - float(ik)/float(nks) * 4.*param.pi_f/3.;
				momenta(ind,1) = 0.0;
				momenta(ind,2) = 0.0;
				ind += 1;
			}
		}

		void set_momenta_Path4_b() { //Gamma -> (0.5,0) -> (1/3,1/3) -> Gamma in terms of rec. lattice b vectors

			size_t nks(nk/3);
			size_t ind(0);
			for (size_t ik=0; ik<nks; ik++) { // Gamma = 0 -> (0.5,0) = 0.5*b1
				momenta(ind,0) = float(ik)/float(nks) * 0.5*(*this).b(0,0);
				momenta(ind,1) = float(ik)/float(nks) * 0.5*(*this).b(0,1);
				momenta(ind,2) = float(ik)/float(nks) * 0.5*(*this).b(0,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // (0.5,0) --> (1/3, 1/3): 0.5*b1 --> 1/3*b1+1/3*b2
				momenta(ind,0) = (0.5 - float(ik)/float(nks) * (0.5 - 1./3.))*(*this).b(0,0) + float(ik)/float(nks)*1./3.*(*this).b(1,0);
				momenta(ind,1) = (0.5 - float(ik)/float(nks) * (0.5 - 1./3.))*(*this).b(0,1) + float(ik)/float(nks)*1./3.*(*this).b(1,1);
				momenta(ind,2) = (0.5 - float(ik)/float(nks) * (0.5 - 1./3.))*(*this).b(0,2) + float(ik)/float(nks)*1./3.*(*this).b(1,2);
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // (1/3,1/3) --> Gamma
				momenta(ind,0) = 1./3.*((*this).b(0,0)+(*this).b(1,0)) * (1. - float(ik)/float(nks));
				momenta(ind,1) = 1./3.*((*this).b(0,1)+(*this).b(1,1)) * (1. - float(ik)/float(nks));
				momenta(ind,2) = 1./3.*((*this).b(0,2)+(*this).b(1,2)) * (1. - float(ik)/float(nks));
				ind += 1;
			}
		}

		void set_momenta_Path5() { //Gamma -> M -> K -> Gamma for triangular lattice in extended BZ

			size_t nks(nk/3);
			size_t ind(0);
			for (size_t ik=0; ik<nks; ik++) { // Gamma = 0 -> M = (2*pi, -2*pi/sqrt(3))
				momenta(ind,0) = float(ik)/float(nks) * 2.*param.pi_f;
				momenta(ind,1) = float(ik)/float(nks) * (-2.*param.pi_f/sqrt(3.));
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // M -> K = (8*pi/3, 0)
				momenta(ind,0) = 2.*param.pi_f + float(ik)/float(nks) * 2.*param.pi_f/3.;
				momenta(ind,1) = -2.*param.pi_f/sqrt(3.) + float(ik)/float(nks) * 2.*param.pi_f/sqrt(3.);
				momenta(ind,2) = 0.0;
				ind += 1;
			}
			for (size_t ik=0; ik<nks; ik++) { // K -> Gamma = (0,0)
				momenta(ind,0) = 8.*param.pi_f/3. - float(ik)/float(nks) * 8.*param.pi_f/3.;
				momenta(ind,1) = 0.0;
				momenta(ind,2) = 0.0;
				ind += 1;
			}
		}


		size_t index(size_t ikx,size_t iky,size_t ikz=0) const {
			// return ikz+iky*nkz+ikx*nkz*nk;
			return ikx+iky*nk+ikz*nk*nk;
		}

		// size_t index(size_t ikx,size_t iky) const {
		// 	return iky+ikx*nk;
		// }

		// void kToik(VectorType& k, size_t& ik, Field& residue)  const {
		// 	std::vector<size_t> iki(3,0);
		// 	std::vector<Field > fki(3,0);
		// 	std::vector<size_t> nki(3,0); nki[0] = (*this).nk; nki[1] = (*this).nk; nki[2] = (*this).nkz;
		// 	residue = 0.0;
		// 	mapTo1BZ(k);
		// 	VectorType bi(3,0);
		// 	for(size_t l=0;l<(*this).dim;l++) {
		// 		(*this).b.getRow(l,bi);
		// 		fki[l] = (scalarProd(k,bi)/scalarProd(bi,bi) - shift[l]) * float(nki[l]);
		// 		iki[l] = size_t(fki[l] + 1.0e-5); 
		// 		residue += fabs(fki[l]-float(iki[l]));
		// 		iki[l] = iki[l] % nki[l];
		// 	}

		// 	// std::cout << "k,ik,residue="<<k[0]<<","<<k[1]<<","<<k[2]<<","<<iki[0]<<","<<iki[1]<<","<<iki[2]<<","<<fki[0]<<","<<fki[1]<<","<<fki[2]<<","<<residue<<"\n";
		// 	ik = index(iki[0],iki[1],iki[2]);
		// }

		void kToik(VectorType& k, size_t& ik, Field& residue, bool map1BZ = false) const {
			// std::vector<FieldType> k(kIn);
			std::vector<size_t> iki(3,0);
			std::vector<Field > fki(3,0);
			std::vector<size_t> nki(3,0); nki[0] = (*this).nk; nki[1] = (*this).nk; nki[2] = (*this).nkz;
			residue = 0.0;
			if (map1BZ) mapTo1BZ(k);
			for (size_t l=0; l<(*this).dim; l++) {
				fki[l] = dGInverse(l,0) * k[0] + dGInverse(l,1) * k[1] + dGInverse(l,2) * k[2] - shift[l] * nki[l];	
				iki[l] = size_t(fki[l] + 1.0e-5);
				residue += fabs(fki[l]-float(iki[l]));
				iki[l] = iki[l] % nki[l];
			}
			ik = index(iki[0],iki[1],iki[2]);
		}

		FieldType scalarProd(const VectorType& a, const VectorType& b) const {
			size_t n(a.size());
			FieldType r1(0.0);
			for (size_t i=0;i<n;i++) r1 += a[i]*b[i];
			return r1;
			// return r1/(4.0*pow(param.pi_f,2));
		}

		void mapTo1BZ(VectorType& k) const {
			FieldType small(1.0e-5);
			VectorType bi(3,0);
			for (size_t l=0;l<3;l++) {
				(*this).b.getRow(l,bi);
				FieldType r1 = scalarProd(k,bi)/scalarProd(bi,bi);
				// std::cout << "bi="<<bi[0]<<","<<bi[1]<<"\n";
				while (r1 <= shift[l]-small) {
					// std::cout << "r1="<<r1<<"\n";
					for (size_t i=0;i<3;i++) k[i] += bi[i];
					r1 = scalarProd(k,bi)/scalarProd(bi,bi);
				}
				while (r1 >= 1.0+shift[l]-small) {
					for (size_t i=0;i<3;i++) k[i] -= bi[i];
					r1 = scalarProd(k,bi)/scalarProd(bi,bi);
				}
			}
		}


		// void interpolateBiLinear(const VectorType& q, 
		// 						 const ComplexVectorType& complexField, 
		// 						 ComplexType& result) const {

		// 	// First find indices of vertices of cube surrounding q
		// 	size_t index(0); FieldType residue(0.0);
		// 	kToik(q,index,residue);
			
		// 	// std::cout<<"q,index,q[index]"<<q[0]<<","<<q[1]<<","<<index<<","<<(*this).momenta(index,0)<<","<<(*this).momenta(index,1)<<"\n";

		// 	if(fabs(residue)<=1.0e-5) { // no interpolation needed!
		// 		result = complexField[index];
		// 		return;
		// 	}

		// 	std::vector<size_t> n(2,0);
		// 	n[0] = (*this).index2Components(index,0);
		// 	n[1] = (*this).index2Components(index,1);

		// 	// std::cout<<"q,index,n[0],n[1]"<<q[0]<<","<<q[1]<<","<<index<<"  ,   n: "<<n[0]<<","<<n[1]<<"\n";

		// 	// // for(size_t l=0;l<2;l++) n[l]=size_t((qMesh.nk-1)*q[l]/param.pi_f); // n[l] now contains the indices of the q-point below q
		// 	// for(size_t l=0;l<2;l++) n[l]=size_t((qMesh.nk-1)*(q[l]+param.pi_f)/(2.0*param.pi_f)); // n[l] now contains the indices of the q-point below q

		// 	size_t n0p1((n[0]+1) % ((*this).nk-1));
		// 	size_t n1p1((n[1]+1) % ((*this).nk-1));

		// 	size_t ind00 = (*this).index(n[0],n[1]);
		// 	size_t ind01 = (*this).index(n[0],n1p1);
		// 	size_t ind10 = (*this).index(n0p1,n[1]);
		// 	size_t ind11 = (*this).index(n0p1,n1p1);

		// 	// Now interpolate, see <http://en.wikipedia.org/wiki/Trilinear_interpolation>
		// 	VectorType dis(2,0);
		//     dis[0] = (q[0]-(*this).momenta(ind00,0))/((*this).momenta(ind11,0)-(*this).momenta(ind00,0));
		//     dis[1] = (q[1]-(*this).momenta(ind00,1))/((*this).momenta(ind11,1)-(*this).momenta(ind00,1));
			
		// 	ComplexType c00(0.0); ComplexType c10(0.0); ComplexType  c0(0.0);
		// 	c00 = complexField[ind00] * (1.-dis[0]) + complexField[ind10] * dis[0];
		// 	c10 = complexField[ind01] * (1.-dis[0]) + complexField[ind11] * dis[0];
		// 	result = c00 * (1.-dis[1]) + c10 * dis[1];

		// 	std::cout << "q,q00,q10,q01,q11:  "<<q[0]<<","<<q[1]<<" , "<<(*this).momenta(ind00,0)<<","<<(*this).momenta(ind00,1)<<
		// 										 " , "<<(*this).momenta(ind10,0)<<","<<(*this).momenta(ind10,1)<<
		// 										 " , "<<(*this).momenta(ind01,0)<<","<<(*this).momenta(ind01,1)<<
		// 										 " , "<<(*this).momenta(ind11,0)<<","<<(*this).momenta(ind11,1)<< "\n";
		// }

		void getSurroundingIndices(VectorType& q, 
					    		   size_t& ind00, size_t& ind10, size_t& ind01, size_t& ind11) const {

			// Find indices of vertices of cube surrounding q
			size_t index(0); FieldType residue(0.0);
			kToik(q,index,residue,true);

			std::vector<size_t> n(2,0);
			n[0] = (*this).index2Components(index,0);
			n[1] = (*this).index2Components(index,1);

			size_t n0p1((n[0]+1) % (*this).nk);
			size_t n1p1((n[1]+1) % (*this).nk);

			ind00 = (*this).index(n[0],n[1]);
			ind01 = (*this).index(n[0],n1p1);
			ind10 = (*this).index(n0p1,n[1]);
			ind11 = (*this).index(n0p1,n1p1);
		}

		void getSurroundingIndices(VectorType& q, 
					    		   size_t& ind000, size_t& ind001, size_t& ind010, size_t& ind100, 
					    		   size_t& ind011, size_t& ind101, size_t& ind110, size_t& ind111) const {

			// Find indices of vertices of cube surrounding q
			size_t index(0); FieldType residue(0.0);
			kToik(q,index,residue,true);

			std::vector<size_t> n(3,0);
			n[0] = (*this).index2Components(index,0);
			n[1] = (*this).index2Components(index,1);
			n[2] = (*this).index2Components(index,2);

			size_t n0p1((n[0]+1) % (*this).nk);
			size_t n1p1((n[1]+1) % (*this).nk);
			size_t n2p1((n[2]+1) % (*this).nkz);

			ind000 = (*this).index(n[0],n[1],n[2]);
			ind010 = (*this).index(n[0],n1p1,n[2]);
			ind100 = (*this).index(n0p1,n[1],n[2]);
			ind110 = (*this).index(n0p1,n1p1,n[2]);
			ind001 = (*this).index(n[0],n[1],n2p1);
			ind101 = (*this).index(n0p1,n[1],n2p1);
			ind011 = (*this).index(n[0],n1p1,n2p1);
			ind111 = (*this).index(n0p1,n1p1,n2p1);
		}



		void set_indexOfAdd() {
			for (size_t i1 = 0; i1 < nktot; ++i1)
			{
				for (size_t i2 = 0; i2 < nktot; ++i2)
				{
					indexOfAdd(i1,i2) = calcIndexOfAdd(i1,i2);
					// std::cout << "i1,i2,index: " << i1 << "," << i2 << "," << indexOfAdd(i1,i2) << "\n";
				}
			}
		}	
						 
		size_t calcIndexOfAdd(const size_t i1, const size_t i2) {
			std::vector<size_t> m(3);
			for (size_t i=0; i<3; i++) {
				m[i] = index2Components(i1,i) + index2Components(i2,i);
				if (m[i] > nk-1) m[i] -= nk;
			}
			size_t index(0);
			if (dim==2) index = m[0]+m[1]*nk;
			else if (dim==3) index = m[0]+m[1]*nk*nk+m[2]*nk;
			return index;
		}

		// size_t calcIndexOfAdd_Old(const size_t i1, const size_t i2) {
		// 	// find index of q = k[i1]+k[i2]
		// 	// Note: We assume that reciprocal lattice vectors are given by  b1, b2 and b3
		// 	Field eps(1.0e-5);
		// 	VectorType q(3);
		// 	// std::cout << "q=" << q ;
			
		// 	for (size_t i=0; i<dim; i++) {
		// 		q[i]=momenta(i1,i)+momenta(i2,i);
		// 	}

		// 	mapBack(q);

		// 	for (size_t ik = 0; ik < nktot; ++ik)
		// 	{
		// 		Field diff(0.0);
		// 		for (size_t i=0; i<dim; i++) diff += abs(q[i]-momenta(ik,i));
		// 		// std::cout << "ik=" <<ik<< "," << diff << "\n";
		// 		if (diff < eps) return ik;
		// 	}
		// 	std::cout << "Index not found, bailing out! q="<<q<<"\n"; exit(0);
		// }

		// void mapBack_Old(std::vector<Field>& v) {
		// 	// Projects v back into the 1. BZ
		// 	// Assumes that v can be mapped back into 1. BZ by single translation
		// 	Field eps(1.0e-5);
		// 	Field vProj;
		// 	// std::cout << "v before map " << v << "\n";
		// 	for (size_t i = 0; i < dim; ++i) {
		// 		Field r1=pow(b(i,0),2)+pow(b(i,1),2)+pow(b(i,2),2);
		// 		if (r1==0) r1=eps;
		// 		vProj = (v[0]*b(i,0)+v[1]*b(i,1)+v[2]*b(i,2))/r1;
		// 		// std::cout << "i,vProj " << i << "," << vProj << "\n";
		// 		if (vProj < -eps)        for (size_t j=0; j<dim; j++) v[j] += b(i,j);
		// 		else if (vProj > 1-eps) for (size_t j=0; j<dim; j++) v[j] -= b(i,j);
		// 	}
		// 	// std::cout << "v after map " << v << "\n";

		// }

	};

}

#endif
