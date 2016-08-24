#include "Matrix.h"
// #include "bands.h"
#include "bandstructure.h"
#include "susceptibility.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "utilities.h"
#include "Range.h"
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<double> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<double> ConcurrencyType;
#endif
typedef double FieldType;



#include <vector>
#include "rpa.h"
// #include "greensFunction.h"
// #include "chi0Ofq.h"
#include "chi0.h"
#include "pairing.h"

typedef PsimagLite::Range<ConcurrencyType> RangeType;

template <typename Field,  template<typename> class MatrixTemplate, typename ConcurrencyType> 
void calcBands(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param, ConcurrencyType& conc) {
	rpa::momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
	kmesh.set_momenta(false);	
	rpa::bandstructure<Field,psimag::Matrix,ConcurrencyType> bands(param,conc,kmesh,false);
	// std::vector<FieldType> w(10);
	// ComplexMatrixType v(10,10);
	// std::vector<FieldType> k(3,0.0); k[0] = 3.90704907; k[1] = param.pi_f;
	// bands.getBands(k,w,v);
	// std::cout << "k = " << k << "\n";
	// std::cout << "w = " << w << "\n";
	// std::cout << "v(1,6),v(1,7) = " << v(1,6) << " , " << v(1,7) << "\n";

	bands.calcBandStructure("ek.dat",true);

	// Now calculate bands along high-symmetry direction
	rpa::momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh2(param,conc,"Path2",66);
	// kmesh2.set_momenta_Path2();
	rpa::bandstructure<Field,psimag::Matrix,ConcurrencyType> bands2(param,conc,kmesh2,false);
	bands2.calcBandStructure("ek_high_sym.dat",false);

}

int main(int argc,char *argv[])
{
    if (argc<2) {
        std::cerr<<"At least one argument needed\n";
        return 1;
    }

    std::string filename="";
    filename=argv[1];

	using namespace rpa;


	// std::string fileName(argv[1]);


	// typedef std::complex<FieldType>	 			ComplexType;
	// typedef psimag::Matrix<ComplexType> 	    ComplexMatrixType;
	// typedef psimag::Matrix<FieldType> 	        MatrixType;


	ConcurrencyType concurrency(argc,argv);

	parameters<FieldType,psimag::Matrix,ConcurrencyType> param(concurrency);
	param.readFromInputFile(filename);
	if (concurrency.rank()==0) param.writeParameters(std::cout);
	param.setupOrbitalIndices();


	calcBands(param,concurrency);

	// if(param.options.find("calcBands")!=std::string::npos) {
	// 	if (concurrency.rank()==0) std::cout << "Now calculating Bands \n";
	// 	Bands<FieldType,psimag::Matrix,ConcurrencyType> bands(param,concurrency);
	// 	FieldType kmin(-0.5*param.pi_f); FieldType kmax(1.5*param.pi_f);
	// 	FieldType kzmin(-0.5*param.pi_f); FieldType kzmax(1.5*param.pi_f);
	// 	bands.calcBandsKMesh(65,1,2,kmin,kmax,kmin,kmax,kzmin,kzmax);
	// }

	typedef bandstructure<FieldType,psimag::Matrix,ConcurrencyType> BandsType;
	if (param.options.find("calcFS")!=std::string::npos) {
		ferminator<FieldType,BandsType,psimag::Matrix,ConcurrencyType> FSpoints(param,concurrency,1);
	}

	if (param.options.find("calcSus")!=std::string::npos) {
		if (concurrency.rank()==0) std::cout << "Now calculating susceptibility \n";
		if (concurrency.rank()==0) std::cout << "qxmin,qxmax,qymin,qymax="
											<<param.qxmin <<","<<param.qxmax<<","<<param.qymin<<","<<param.qymax<<"\n";


		// size_t nq=8;
		// momentumDomain<FieldType,psimag::Matrix> qMesh(param,nq,1,2,param.chia1,param.chia2,param.chia3);
		// qMesh.set_momenta(false);
		typedef susc<FieldType,psimag::Matrix,ConcurrencyType> SuscType;
		// chi0q<FieldType,SuscType,BandsType,psimag::Matrix,ConcurrencyType> susq(param,qMesh,concurrency);
		susceptibility<FieldType,SuscType,BandsType,psimag::Matrix,ConcurrencyType> 
					   susq(param,concurrency,
					   	    param.qxmin,param.qxmax,param.nqx,
					   	    param.qymin,param.qymax,param.nqy,
					   	    param.qzmin,param.qzmax,param.nqz,
					   	    param.wmin,param.wmax,param.nw
					   	    );
	}



	if(param.options.find("calcPairing")!=std::string::npos) {
		if (concurrency.rank()==0) std::cout << "Now calculating Pairing \n";
		typedef bandstructure<FieldType,psimag::Matrix,ConcurrencyType> BandsType;
		typedef susc<FieldType,psimag::Matrix,ConcurrencyType> SuscType;
		size_t nq(0),nqz(0);
		if (param.interpolateChi==1) {
			if (param.readChi==1) {
				std::vector<FieldType> data;
				loadVector(data,param.chifile);
				nq=size_t(data[0]);
				nqz=size_t(data[1]);
				std::cout << "nq,nqz="<<nq<<","<<nqz<<"\n";
			} else {
				nq=param.nkInt; nqz=param.nkIntz;
			}
		}
		momentumDomain<FieldType,psimag::Matrix,ConcurrencyType> qMesh(param,concurrency,nq,nqz,param.dimension,param.chia1,param.chia2,param.chia3);
		qMesh.set_momenta(false);
		typedef rpa::gap2D<FieldType,psimag::Matrix,ConcurrencyType> GapType; // this is not really needed
		pairing<FieldType,BandsType,SuscType,GapType,psimag::Matrix,ConcurrencyType> pairing(param,concurrency,param.interpolateChi,qMesh);
	}



	// if(param.options.find("calcTest")!=std::string::npos) {
	// 	susc<FieldType,psimag::Matrix,ConcurrencyType> chi0(param,concurrency);
	// 	typedef bandstructure<FieldType,psimag::Matrix,ConcurrencyType> BandsType;
	// 	momentumDomain<FieldType,psimag::Matrix> kMesh(param,param.nkInt,param.nkIntz,param.dimension);
	// 	kMesh.set_momenta(false);
	// 	BandsType bands(param,concurrency,kMesh);
	// 	std::vector<FieldType> q(3); q[0] = 0.0; q[1] = 0.*param.pi_f; q[2] = 0.*param.pi_f;
	// 	calcChi0Matrix<FieldType,BandsType,psimag::Matrix,ConcurrencyType> calcChi0(param,kMesh,bands,q,concurrency);
	// 	size_t msize(param.nOrb*param.nOrb);
	// 	for (size_t i=0;i<msize;i++) for (size_t j=i;j<msize;j++) chi0(i,j) = calcChi0(i,j);
	// 	chi0.setLowerTriangle();
	// 	std::cout << "q,Chi0:" << q << ", " << chi0.calcSus() << "\n";
	// }

	// greensFunction<FieldType> g(param);
	// chi0ofq<FieldType,ConcurrencyType> chi0(param,g,concurrency);
	// if (concurrency.rank()==0) {
	// 	std::ofstream os("chi0ofq.dat");
	// 	os << chi0.chiq;
	// 	os.close();
	// }



	// size_t nq=32,dim=2;
	// momentumDomain<FieldType,psimag::Matrix> momenta(param,nq,dim);
	// size_t nq=64;
	// std::string path("Path");
	// momentumDomain<FieldType,psimag::Matrix> momenta(param,path,nq);
	// // momenta.set_momenta_standard();
	// momenta.set_momenta_Path1();
	// // momenta.set_momenta(false);
	// susceptibility<FieldType,psimag::Matrix> chi0(param);



	// RangeType range(0,momenta.nktot,concurrency);
	// std::vector<FieldType> chi0q(momenta.nktot);

	// interaction<FieldType,psimag::Matrix> rpa(param);
	// std::vector<FieldType> chiUq(momenta.nktot);

	// size_t msize(param.nOrb*param.nOrb);
	// psimag::Matrix<std::complex<FieldType> > chiRPA(msize,msize);


	// for (;!range.end();range.next()) {
		
	// 		size_t iq = range.index();
	// 		std::vector<FieldType> q(param.dimension);
	// 		for (size_t i=0; i<param.dimension; i++) q[i] = momenta.momenta(iq,i);

	// 		if (concurrency.rank()==0) std::cout << "iq = " << iq << " of " << momenta.nktot << " total." << "\n";
	// 		chi0.chi0OfQ(q);
	// 		chi0q[iq] = chi0.calcChiPhys(chi0.chi0,q);

	// 		rpa.calcRPAResult(chi0.chi0,rpa.spinMatrix,chiRPA);
	// 		chiUq[iq]=chi0.calcChiPhys(chiRPA,q);
			
	// 		// std::cout << "chi0: " << chi0.chi0(0,0) <<  ", chiU: " << chiUq[iq] << "\n"; 

	// 	}
	// 	concurrency.reduce(chi0q);
	// 	concurrency.reduce(chiUq);
	// 	if (concurrency.rank()==0) {
	// 		std::ofstream os("chi0Ofq.dat");
	// 		std::ofstream os2("chiUOfq.dat");
	// 		os << chi0q;
	// 		os2 << chiUq;
	// 		os.close();
	// 		os2.close();
	// 	}

		return 0;

}


