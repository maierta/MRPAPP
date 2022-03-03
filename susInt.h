//-*-C++-*-

#ifndef SUSINT_H
#define SUSINT_H

#include "math.h"

namespace rpa {

	// General Function needed in suscept. calculation
	template<typename FieldType>
	inline FieldType susInt(const FieldType& e1, const FieldType& e2, 
						const FieldType& invT) {
		FieldType sus(0);
		if (e1!=e2)
		{
			FieldType xx1(e1*invT);
			FieldType xx2(e2*invT);
			FieldType r1 = PsimagLite::fermi(xx1)-PsimagLite::fermi(xx2);
			sus = r1/(e2-e1);
		} else 
		{
			FieldType r1(2.*cosh(0.5*e1*invT));
			sus = invT/(r1*r1);
		}
		return sus;
	}

	template<typename FieldType>
	inline std::complex<FieldType> susInt(const FieldType& e1, const FieldType& e2, 
					                	  const FieldType& invT, const FieldType& omega,
					                	  const FieldType& damp=FieldType(1.0e-3)) {
		if (omega==0) return susInt(e1, e2, invT);

		std::complex<FieldType> sus(0);
		FieldType xx1(e1*invT);
		FieldType xx2((e2)*invT);
		FieldType r1 = PsimagLite::fermi(xx1)-PsimagLite::fermi(xx2);
		sus = r1/(e2-e1+omega+std::complex<FieldType>(0.0,damp));
		return sus;
	}

	// template<typename FieldType>
	// inline std::complex<FieldType> susIntBCS(const FieldType& e1, const FieldType& e2,
	// 									     const std::complex<FieldType>& gap1, 
	// 									     const std::complex<FieldType>& gap2, 
	// 					                	 const FieldType& invT, const FieldType& omega,
	// 				                	     const FieldType& damp=FieldType(1.0e-3),
	// 				                	     const FieldType& signF = -1) {
	// 	std::complex<FieldType> sus(0);

	// 	FieldType r1(std::norm(gap1));
	// 	FieldType r2(std::norm(gap2));

	// 	FieldType EnergyBCS1(sqrt(pow(e1,2)+r1));
	// 	FieldType EnergyBCS2(sqrt(pow(e2,2)+r2));

	// 	r1 = e1/EnergyBCS1;
	// 	r2 = e2/EnergyBCS2;

	// 	FieldType uk1(0),vk1(0),uk2(0),vk2(0);
		
	// 	if (EnergyBCS1==0.0) {uk1=0.0; vk1=1.0;}
	// 	else {uk1 = 0.5*(1.0+r1); vk1 = 0.5*(1.0-r1);}

	// 	if (EnergyBCS2==0.0) {uk2=0.0; vk2=1.0;}
	// 	else {uk2 = 0.5*(1.0+r2); vk2 = 0.5*(1.0-r2);}

	// 	r1 = -signF * 0.25 * real(gap1*conj(gap2))/(EnergyBCS1*EnergyBCS2);
		
	// 	FieldType skq1(uk1*uk2+r1);
	// 	FieldType skq2(vk1*vk2+r1);
	// 	FieldType skq3(uk1*vk2-r1);
	// 	FieldType skq4(vk1*uk2-r1);

	// 	sus = susInt(+EnergyBCS1,+EnergyBCS2,invT,omega,damp) * skq1 
	// 		+ susInt(-EnergyBCS1,-EnergyBCS2,invT,omega,damp) * skq2
	// 		+ susInt(+EnergyBCS1,-EnergyBCS2,invT,omega,damp) * skq3
	// 		+ susInt(-EnergyBCS1,+EnergyBCS2,invT,omega,damp) * skq4;

	// 	return sus;
	// }

	template<typename FieldType>
	inline std::complex<FieldType> susIntBCS(const FieldType& e1, const FieldType& e2,
						 const std::complex<FieldType>& gap1, 
					         const std::complex<FieldType>& gap2, 
					         const FieldType& invT, const FieldType& omega,
					         const FieldType& damp=FieldType(1.0e-3),
					         const FieldType& signF = 1) {
		std::complex<FieldType> sus(0);
		FieldType uk1(0),vk1(0),uk2(0),vk2(0),r1(0),r2(0);
		FieldType Delta1s(std::norm(gap1));
		FieldType Delta2s(std::norm(gap2));


		if (Delta1s == 0 && Delta2s == 0) {

			return susInt(e1,e2,invT,omega,damp);

		} else if (Delta1s == 0 && Delta2s != 0) {
	
			FieldType EnergyBCS2(sqrt(pow(e2,2)+Delta2s));
			r2 = e2 / EnergyBCS2;
			uk2 = 0.5*(1.0+r2); vk2 = 0.5*(1.0-r2);
			return susInt(e1,EnergyBCS2,invT,omega,damp) * uk2 + susInt(e1,-EnergyBCS2,invT,omega,damp) * vk2;

		} else if (Delta1s != 0 && Delta2s == 0) {
	
			FieldType EnergyBCS1(sqrt(pow(e1,2)+Delta1s));
			r1 = e1 / EnergyBCS1;
			uk1 = 0.5*(1.0+r1); vk1 = 0.5*(1.0-r1);
			return susInt(EnergyBCS1,e2,invT,omega,damp) * uk1 + susInt(-EnergyBCS1,e2,invT,omega,damp) * vk1;

		} else {

			FieldType EnergyBCS1(sqrt(pow(e1,2)+Delta1s));
			FieldType EnergyBCS2(sqrt(pow(e2,2)+Delta2s));
			r1 = e1 / EnergyBCS1; r2 = e2 / EnergyBCS2;
			uk1 = 0.5*(1.0+r1); vk1 = 0.5*(1.0-r1);
			uk2 = 0.5*(1.0+r2); vk2 = 0.5*(1.0-r2);

			r1 = signF * 0.25 * real(gap1*conj(gap2))/(EnergyBCS1*EnergyBCS2);

			FieldType skq1(uk1*uk2+r1);
			FieldType skq2(vk1*vk2+r1);
			FieldType skq3(uk1*vk2-r1);
			FieldType skq4(vk1*uk2-r1);

			sus = susInt(+EnergyBCS1,+EnergyBCS2,invT,omega,damp) * skq1 
				+ susInt(-EnergyBCS1,-EnergyBCS2,invT,omega,damp) * skq2
				+ susInt(+EnergyBCS1,-EnergyBCS2,invT,omega,damp) * skq3
				+ susInt(-EnergyBCS1,+EnergyBCS2,invT,omega,damp) * skq4;

			return sus;

		}

	}

	template<typename FieldType>
	inline std::complex<FieldType> susIntBCSGG(const FieldType& e1, const FieldType& e2,
											   const std::complex<FieldType>& gap1, 
											   const std::complex<FieldType>& gap2, 
											   const FieldType& invT, const FieldType& omega,
											   const FieldType& damp=FieldType(1.0e-3)) {
		
		std::complex<FieldType> sus(0);
		FieldType uk1(0),vk1(0),uk2(0),vk2(0),r1(0),r2(0);
		FieldType Delta1s(std::norm(gap1));
		FieldType Delta2s(std::norm(gap2));


		if (Delta1s == 0 && Delta2s == 0) {

			if (omega==0) return susInt(e1,e2,invT);
			else return susInt(e1,e2,invT,omega,damp);

		} else if (Delta1s == 0 && Delta2s != 0) {
	
			FieldType EnergyBCS2(sqrt(pow(e2,2)+Delta2s));
			r2 = e2 / EnergyBCS2;
			uk2 = 0.5*(1.0+r2); vk2 = 0.5*(1.0-r2);
			if (omega==0) return susInt(e1,EnergyBCS2,invT) * uk2 + susInt(e1,-EnergyBCS2,invT) * vk2;
			else return susInt(e1,EnergyBCS2,invT,omega,damp) * uk2 + susInt(e1,-EnergyBCS2,invT,omega,damp) * vk2;

		} else if (Delta1s != 0 && Delta2s == 0) {
	
			FieldType EnergyBCS1(sqrt(pow(e1,2)+Delta1s));
			r1 = e1 / EnergyBCS1;
			uk1 = 0.5*(1.0+r1); vk1 = 0.5*(1.0-r1);
			if (omega==0) return susInt(EnergyBCS1,e2,invT) * uk1 + susInt(-EnergyBCS1,e2,invT) * vk1;
			else return susInt(EnergyBCS1,e2,invT,omega,damp) * uk1 + susInt(-EnergyBCS1,e2,invT,omega,damp) * vk1;

		} else {

			FieldType EnergyBCS1(sqrt(pow(e1,2)+Delta1s));
			FieldType EnergyBCS2(sqrt(pow(e2,2)+Delta2s));
			r1 = e1 / EnergyBCS1; r2 = e2 / EnergyBCS2;
			uk1 = 0.5*(1.0+r1); vk1 = 0.5*(1.0-r1);
			uk2 = 0.5*(1.0+r2); vk2 = 0.5*(1.0-r2);

			FieldType skq1(uk1*uk2);
			FieldType skq2(vk1*vk2);
			FieldType skq3(uk1*vk2);
			FieldType skq4(vk1*uk2);

			if (omega==0) { 
				sus = susInt(+EnergyBCS1,+EnergyBCS2,invT,omega,damp) * skq1 
					+ susInt(-EnergyBCS1,-EnergyBCS2,invT,omega,damp) * skq2
					+ susInt(+EnergyBCS1,-EnergyBCS2,invT,omega,damp) * skq3
					+ susInt(-EnergyBCS1,+EnergyBCS2,invT,omega,damp) * skq4;
			} else {
				sus = susInt(+EnergyBCS1,+EnergyBCS2,invT) * skq1 
					+ susInt(-EnergyBCS1,-EnergyBCS2,invT) * skq2
					+ susInt(+EnergyBCS1,-EnergyBCS2,invT) * skq3
					+ susInt(-EnergyBCS1,+EnergyBCS2,invT) * skq4;
			}


			return sus;

		}

	}


	template<typename FieldType>
	inline std::complex<FieldType> susIntBCSFF(const FieldType& e1, const FieldType& e2,
						     				   const std::complex<FieldType>& gap1, 
					         				   const std::complex<FieldType>& gap2, 
					         				   const FieldType& invT, const FieldType& omega,
					         				   const FieldType& damp=FieldType(1.0e-3),
					         				   const FieldType& signF = 1) {
		
		std::complex<FieldType> sus(0);
		std::complex<FieldType> r1(0);
		FieldType Delta1s(std::norm(gap1));
		FieldType Delta2s(std::norm(gap2));


		if (Delta1s == 0 || Delta2s == 0) {

			return std::complex<FieldType>(0.0,0.0);


		} else {

			FieldType EnergyBCS1(sqrt(pow(e1,2)+Delta1s));
			FieldType EnergyBCS2(sqrt(pow(e2,2)+Delta2s));
			// std::cout << "signF " << signF << "\n";
			r1 = signF * 0.25 * conj(gap1)*gap2/(EnergyBCS1*EnergyBCS2);
			// r1 = signF * 0.25 * real(conj(gap1)*gap2)/(EnergyBCS1*EnergyBCS2);
			// std::cout << "gap1="<<gap1<<"\n";
			// r1 = 1;

			sus = susInt(+EnergyBCS1,+EnergyBCS2,invT,omega,damp) 
				+ susInt(-EnergyBCS1,-EnergyBCS2,invT,omega,damp)
				- susInt(+EnergyBCS1,-EnergyBCS2,invT,omega,damp)
				- susInt(-EnergyBCS1,+EnergyBCS2,invT,omega,damp);

			return sus * r1;

		}

	}


}

#endif
