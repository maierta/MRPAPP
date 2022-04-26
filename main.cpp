#include "Matrix.h"
// #include "bands.h"
#include "bandstructure.h"
#include "susceptibility.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "utilities.h"
#include "model.h"
#include "Range.h"
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<double> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<double> ConcurrencyType;
#endif




#include <vector>
#include "chi0.h"
#include "pairing.h"

typedef double FieldType;
typedef PsimagLite::Range<ConcurrencyType> RangeType;
typedef rpa::model<FieldType, psimag::Matrix, ConcurrencyType> ModelType;


template <typename Field,  template<typename> class MatrixTemplate, typename ModelType, typename ConcurrencyType> 
void calcBands(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param, ModelType& model, ConcurrencyType& conc) {

	rpa::momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkBands,param.nkIntz,param.dimension);
	kmesh.set_momenta(false);	

	rpa::bandstructure<Field,psimag::Matrix, ModelType, ConcurrencyType> bands(param,model,conc,kmesh,false);


	// Simple test
	// std::vector<FieldType> ek(param.nOrb,0), k(3);
	// ComplexMatrixType ak(param.nOrb,param.nOrb);
        //
	// k[0] = 0.6283185307;
	// k[1] = 0.9424777961;
	// k[2] = 1.2566370614;
	//
	// bands.getEkAndAk(k, ek, ak);
	// std::cout << ek << "\n";


	FieldType filling(bands.calcFilling());

	if (conc.rank()==0) std::cout << "\n\nFilling = " << filling << " \n";
	if (conc.rank()==0) std::cout << "Target filling = " << param.nTarget << " \n";

	if (param.adjustChemicalPotential) {
		// std::cout << "adjusting chem. pot." << abs(filling) - param.nTarget << "\n";

		FieldType delta = 1.0e-3;
		while (abs(filling - param.nTarget) > delta) {
			// adjust chemical potential to yield target density
			param.mu -= 0.1*(filling - param.nTarget);
			// if (conc.rank()==0) std::cout << "mu adjusted to " << param.mu << "\n";
			filling = bands.calcFilling();
		}

		if (conc.rank()==0) {
			std::cout << "\n\nChemical potential = " << param.mu << " \n";
			std::cout << "Filling = " << filling << " \n\n\n";
		}
	}

	// Now that chemical potential is fixed to give target filling, calculate bandstructure
	std::string filename = "ek_" + param.fileID + ".txt";
	bands.calcBandStructure(filename);

	// Now calculate bands along high-symmetry direction
	std::string path("Path2");
	size_t nkPath(1080);
	if (param.dimension == 3) {
		path = "Path3";
		nkPath = 2520;
	}
	rpa::momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh2(param,conc,path,nkPath);

	// kmesh2.set_momenta_Path2();
	rpa::bandstructure<Field,psimag::Matrix,ModelType,ConcurrencyType> bands2(param,model,conc,kmesh2,false);
	std::string filename2 = "ek_high_sym_" + param.fileID + ".txt";
	bands2.calcBandStructure(filename2);

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
	param.setupOrbitalIndices();

	if (concurrency.rank()==0) std::cout << "Now setting up model\n";
	ModelType model(param, concurrency);
	if (concurrency.rank()==0) std::cout << "Model has been setup\n";
	calcBands(param,model,concurrency);

	// if(param.options.find("calcBands")!=std::string::npos) {
	// 	if (concurrency.rank()==0) std::cout << "Now calculating Bands \n";
	// 	calcBands(param,concurrency);
	// // 	Bands<FieldType,psimag::Matrix,ConcurrencyType> bands(param,concurrency);
	// // 	FieldType kmin(-0.5*param.pi_f); FieldType kmax(1.5*param.pi_f);
	// // 	FieldType kzmin(-0.5*param.pi_f); FieldType kzmax(1.5*param.pi_f);
	// // 	bands.calcBandsKMesh(65,1,2,kmin,kmax,kmin,kmax,kzmin,kzmax);
	// }

	typedef bandstructure<FieldType,psimag::Matrix,ModelType,ConcurrencyType> BandsType;
	if (param.options.find("calcFS")!=std::string::npos) {
		ferminator<FieldType,BandsType,psimag::Matrix,ModelType,ConcurrencyType> FSpoints(param,model,concurrency,1);
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
		susceptibility<FieldType,SuscType,BandsType,psimag::Matrix,ModelType,ConcurrencyType> 
						   susq(param,model,concurrency,
					   	    param.qxmin,param.qxmax,param.nqx,
					   	    param.qymin,param.qymax,param.nqy,
					   	    param.qzmin,param.qzmax,param.nqz,
					   	    param.wmin,param.wmax,param.nw
					   	    );
	}



	if(param.options.find("calcPairing")!=std::string::npos) {
		if (concurrency.rank()==0) std::cout << "Now calculating Pairing \n";
		typedef bandstructure<FieldType,psimag::Matrix,ModelType,ConcurrencyType> BandsType;
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
		pairing<FieldType,BandsType,SuscType,GapType,psimag::Matrix,ModelType,ConcurrencyType> pairing(param,model,concurrency,param.interpolateChi,qMesh);
	}

	if (concurrency.rank()==0) {
		std::string cstr = "parameters_" + param.fileID + ".txt";
		const char *filename1 = cstr.c_str();
		std::ofstream os1(filename1);
		param.writeParameters(os1);
		os1.close();
	}
	// concurrency.barrier();


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


