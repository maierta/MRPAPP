//-*-C++-*-

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include "ConcurrencyMpi.h"



namespace rpa {
	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class parameters {
	
	private:
		ConcurrencyType& conc;
	public:
		typedef MatrixTemplate<Field> 			MatrixType;
		typedef MatrixTemplate<size_t> 			IntMatrixType;

		Field temperature;
		const Field pi_f;
		Field U,Up,J,Jp,U_d_s,U_d_c,U_p_s,U_p_c,U_pd_s,U_pd_c,U_pp_s,U_pp_c,U_d_coupl,U_p_coupl,U_pd_coupl,U_pp_coupl;
		Field staticUFactor,chargeFactor,spinFactor;
		std::vector<Field> deltaU;
		std::vector<Field> a1,a2,a3;
		std::vector<Field> chia1,chia2,chia3;
		size_t nqx,nqy,nqz;
		Field qxmin,qxmax,qymin,qymax,qzmin,qzmax;
		size_t nw;
		Field wmin,wmax;
		size_t scState;
		size_t printGap;
		std::string gAmpl;
		Field Delta0;
		Field Omega0;
		Field signF; // sign of FF term in BCS chi0 calculation
		size_t nwn;
		std::string Case;
		size_t dimension;
		size_t nkInt; 
		size_t nkIntz; 
		Field kz2D;
		size_t FSnkz;
		size_t nOrb;
		Field  mu;
		std::string tbfile;
		std::string fsfile;
		bool readFSFromFile;
		size_t nkPerSheet;
		bool complexHopping;
		size_t pairingSpinParity;
		size_t pairingFromSpin;
		size_t pairingFromCharge;
		size_t storeChi;
		size_t storeGammaOrb;
		size_t readChi;
		size_t readChiForSus;
		std::string chifile;
		size_t interpolateChi;
		size_t interpolateNqx;
		size_t interpolateNqz;
		size_t sublattice;
		Field deltax;
		Field deltay;
		Field deltaz;
		size_t kTrafo;
		std::string options;
		std::string subOptions;
		MatrixType WTrafo;
		IntMatrixType indexToOrb;
		std::vector<int> orbToSite;
		bool hyb;
		bool LS;
		Field hybStrength;
		size_t hybBand1;
		size_t hybBand2;
		size_t nSitesPerUnitCell;
		std::vector<int> nOrbAtom;
		std::string nOrbAtomStr;
		Field damp;
		bool calcOnlyDiagonal;
		bool writeFullChi0;
		bool fixEvecs;

		// size_t nktot;


		parameters(ConcurrencyType& concurrency):
			conc(concurrency), 
			temperature(1.0),
			pi_f(3.141592653589793238462643383),
			U(1.1),
			Up(0.55),
			J(0.275),
			Jp(0.275),
			U_d_s(0.0),
			U_d_c(0.0),
			U_p_s(0.0),
			U_p_c(0.0),
			U_pd_s(0.0),
			U_pd_c(0.0),
			U_pp_s(0.0),
			U_pp_c(0.0),
			U_d_coupl(0.0),
			U_p_coupl(0.0),
			U_pd_coupl(0.0),
			U_pp_coupl(0.0),
			staticUFactor(1.0),
			chargeFactor(1.0),
			spinFactor(1.0),
			deltaU(5,0.0),
			a1(3,0),
			a2(3,0),
			a3(3,0),
			chia1(3,0),
			chia2(3,0),
			chia3(3,0),
			nqx(1),
			nqy(1),
			nqz(1),
			qxmin(0.0),
			qxmax(0.0),
			qymin(0.0),
			qymax(0.0),
			qzmin(0.0),
			qzmax(0.0),
			nw(1),
			wmin(0.0),
			wmax(0.0),
			scState(0),
			printGap(0),
			gAmpl("LaOFeAs_s_1"),
			Delta0(0.04),
			Omega0(0.1),
			signF(-1.0),
			nwn(100),
			Case(""),
			dimension(2),
			nkInt(64),
			nkIntz(16),
			kz2D(0.0),
			FSnkz(10),
			nOrb(1),
			mu(-0.1),
			tbfile(""),
			fsfile("FSforPairing.dat"),
			readFSFromFile(0),
			nkPerSheet(40),
			complexHopping(0),
			pairingSpinParity(0),
			pairingFromSpin(1),
			pairingFromCharge(1),
			storeChi(0),
			storeGammaOrb(0),
			readChi(0),
			readChiForSus(0),
			chifile("none"),
			interpolateChi(0),
			interpolateNqx(17),
			interpolateNqz(5),
			sublattice(0),
			deltax(0.0), 
			deltay(0.0), 
			deltaz(0.0),
			kTrafo(0),
			options(""),
			subOptions(""),
			WTrafo(3,3),
			indexToOrb(0,0),
			orbToSite(0),
			hyb(0),
			LS(0),
			hybStrength(0.0),
			hybBand1(0),
			hybBand2(0),
			nSitesPerUnitCell(1),
			nOrbAtom(0,0),
			nOrbAtomStr(""),
			damp(1.0e-3),
			calcOnlyDiagonal(0),
			writeFullChi0(0),
			fixEvecs(0)

			// single-band model in 2-sub-lattice formulation
			// dimension(2),
			// nkInt(64),
			// nOrb(2),
			// mu(-0.1),
			// tbfile("1band_AB_t.csv"),
			// sublattice(1),
			// deltax(0.0), 
			// deltay(0.0), 
			// deltaz(0.0)

			// 2-orbital model (by Raghu)
			// dimension(2),
			// nkInt(64),
			// nOrb(2),
			// mu(1.45),
			// tbfile("2band.csv"),
			// sublattice(0),
			// deltax(0.0), 
			// deltay(0.0), 
			// deltaz(0.0)

			// 10-orbital model for LaOFeAs
			 // dimension(3),
			 // nkInt(4),
			 // nOrb(10),
			 // mu(11.60),
			 // tbfile("tij_LaOFeAs.csv"),
			 // sublattice(1),
			 // deltax(0.5), 
			 // deltay(0.5), 
			 // deltaz(0.0)

			// 10-orbital model for FeSe
			// dimension(2),
			// nkInt(8),
			// nOrb(10),
			// mu(0.0),
			// tbfile("fese_ab_strain_dx_00_000pc.dat"),
			// sublattice(1),
			// deltax(0.0), 
			// deltay(0.0), 
			// deltaz(0.0)

			{
				// 1 site / unit cell BZ
				a1[0]=1.0; a2[1]=1.0; a3[2]=1.0;
				chia1=a1; chia2=a2; chia3=a3;
				// 2 sites / unit cell BZ
				// a1[0]= 1.0  ; a1[1]=1.0;
				// a2[0]= -1.0 ; a2[1]=1.0;
				// a3[2]=1.0;

			};

			void readFromInputFile(const std::string& file) {
				if (conc.rank()==0) {
					std::ifstream data(file.c_str());
					
					std::string line;
				    while(std::getline(data,line)) {
				        std::stringstream    str(line);
				        std::string          text;
				        std::getline(str,text,'=');
				        setParamBasedOnText(text,str);
					}
				}

				
				conc.barrier();
				broadcastParam();
		        loadVector(nOrbAtom,nOrbAtomStr);
				// bcTest();

				if(kTrafo==0) {
					for (int i = 0; i < 3; ++i)
					{
						WTrafo(i,i) = 1.0;
					}
				} else if(kTrafo==1) {
					WTrafo(0,0) =  0.5; WTrafo(0,1) = 0.5; WTrafo(0,2) = -0.5;
					WTrafo(1,0) = -0.5; WTrafo(1,1) = 0.5; WTrafo(1,2) = -0.5;
					WTrafo(2,0) =  0.5; WTrafo(2,1) = 0.5; WTrafo(2,2) =  0.5;
					// WTrafo(0,0) =  -0.5   ; WTrafo(0,1) =  0.5   ; WTrafo(0,2) = 0.5;
					// WTrafo(1,0) =   0.5   ; WTrafo(1,1) = -0.5   ; WTrafo(1,2) = 0.5;
					// WTrafo(2,0) =   1.8039; WTrafo(2,1) =  1.8039; WTrafo(2,2) =-1.8039;
				}

		        if (dimension==2) {
		        	nkIntz=1;
		        	kz2D *= pi_f;
		        }
			}

			void setParamBasedOnText(std::string& text, std::stringstream& str) {
		        if      (text.find("dimension")!=std::string::npos) str >> (*this).dimension; 
		        else if (text.find("temperature")!=std::string::npos) str >> (*this).temperature; 
		        else if (text.find("numberOfOrbitals")!=std::string::npos) str >> (*this).nOrb; 
		        else if (text.find("chemicalPotential")!=std::string::npos) str >> (*this).mu; 
		        else if (text.find("tbParametersFile")!=std::string::npos) str >> (*this).tbfile; 
		        else if (text.find("complexHopping")!=std::string::npos) str >> (*this).complexHopping; 
		        else if (text.find("FSforPairingFile")!=std::string::npos) str >> (*this).fsfile; 
		        else if (text.find("ChiForPairingFile")!=std::string::npos) str >> (*this).chifile; 
		        else if (text.find("interpolateChi")!=std::string::npos) str >> (*this).interpolateChi; 
		        else if (text.find("interpolateNqx")!=std::string::npos) str >> (*this).interpolateNqx; 
		        else if (text.find("interpolateNqz")!=std::string::npos) str >> (*this).interpolateNqz; 
		        else if (text.find("Coulomb1U")!=std::string::npos) str >> (*this).U; 
		        else if (text.find("Coulomb2Up")!=std::string::npos) str >> (*this).Up; 
		        else if (text.find("Coulomb3J")!=std::string::npos) str >> (*this).J; 
		        else if (text.find("Coulomb4Jp")!=std::string::npos) str >> (*this).Jp; 
		        else if (text.find("Coulomb5U_d")!=std::string::npos) str >> (*this).U_d_c; 
		        else if (text.find("Coulomb6U_p")!=std::string::npos) str >> (*this).U_p_c; 
		        else if (text.find("Coulomb7U_pd")!=std::string::npos) str >> (*this).U_pd_c; 
		        else if (text.find("Coulomb8U_pp")!=std::string::npos) str >> (*this).U_pp_c; 
		        else if (text.find("Coulomb9U_d")!=std::string::npos) str >> (*this).U_d_s; 
		        else if (text.find("Coulomb10U_p")!=std::string::npos) str >> (*this).U_p_s; 
		        else if (text.find("Coulomb11U_pd")!=std::string::npos) str >> (*this).U_pd_s; 
		        else if (text.find("Coulomb12U_pp")!=std::string::npos) str >> (*this).U_pp_s; 
		        else if (text.find("Coulomb13U_d_coupl")!=std::string::npos) str >> (*this).U_d_coupl; 
		        else if (text.find("Coulomb14U_p_coupl")!=std::string::npos) str >> (*this).U_p_coupl; 
		        else if (text.find("Coulomb15U_pd_coupl")!=std::string::npos) str >> (*this).U_pd_coupl; 
		        else if (text.find("Coulomb16U_pp_coupl")!=std::string::npos) str >> (*this).U_pp_coupl; 
		        else if (text.find("sublattice")!=std::string::npos) str >> (*this).sublattice; 
		        else if (text.find("deltaU0")!=std::string::npos) str >> (*this).deltaU[0]; 
		        else if (text.find("deltaU1")!=std::string::npos) str >> (*this).deltaU[1]; 
		        else if (text.find("deltaU2")!=std::string::npos) str >> (*this).deltaU[2]; 
		        else if (text.find("deltaU3")!=std::string::npos) str >> (*this).deltaU[3]; 
		        else if (text.find("deltaU4")!=std::string::npos) str >> (*this).deltaU[4]; 
		        else if (text.find("staticUFactor")!=std::string::npos) str >> (*this).staticUFactor; 				
		        else if (text.find("chargeFactor")!=std::string::npos) str >> (*this).chargeFactor; 				
		        else if (text.find("spinFactor")!=std::string::npos) str >> (*this).spinFactor; 				
		        else if (text.find("nkIntegration")!=std::string::npos) str >> (*this).nkInt; 				
		        else if (text.find("nkzIntegration")!=std::string::npos) str >> (*this).nkIntz; 				
		        else if (text.find("kz2D")!=std::string::npos) str >> (*this).kz2D; 				
		        else if (text.find("a1x")!=std::string::npos) str >> (*this).a1[0]; 				
		        else if (text.find("a1y")!=std::string::npos) str >> (*this).a1[1]; 				
		        else if (text.find("a1z")!=std::string::npos) str >> (*this).a1[2]; 				
		        else if (text.find("a2x")!=std::string::npos) str >> (*this).a2[0]; 				
		        else if (text.find("a2y")!=std::string::npos) str >> (*this).a2[1]; 				
		        else if (text.find("a2z")!=std::string::npos) str >> (*this).a2[2]; 				
		        else if (text.find("a3x")!=std::string::npos) str >> (*this).a3[0]; 				
		        else if (text.find("a3y")!=std::string::npos) str >> (*this).a3[1]; 				
		        else if (text.find("a3z")!=std::string::npos) str >> (*this).a3[2]; 				
		        else if (text.find("c1x")!=std::string::npos) str >> (*this).chia1[0]; 				
		        else if (text.find("c1y")!=std::string::npos) str >> (*this).chia1[1]; 				
		        else if (text.find("c1z")!=std::string::npos) str >> (*this).chia1[2]; 				
		        else if (text.find("c2x")!=std::string::npos) str >> (*this).chia2[0]; 				
		        else if (text.find("c2y")!=std::string::npos) str >> (*this).chia2[1]; 				
		        else if (text.find("c2z")!=std::string::npos) str >> (*this).chia2[2]; 				
		        else if (text.find("c3x")!=std::string::npos) str >> (*this).chia3[0]; 				
		        else if (text.find("c3y")!=std::string::npos) str >> (*this).chia3[1]; 				
		        else if (text.find("c3z")!=std::string::npos) str >> (*this).chia3[2]; 				
		        else if (text.find("nqx")!=std::string::npos) str >> (*this).nqx; 				
		        else if (text.find("nqy")!=std::string::npos) str >> (*this).nqy; 				
		        else if (text.find("nqz")!=std::string::npos) str >> (*this).nqz; 				
		        else if (text.find("qxmin")!=std::string::npos) {str >> (*this).qxmin; qxmin *= pi_f;}
		        else if (text.find("qxmax")!=std::string::npos) {str >> (*this).qxmax; qxmax *= pi_f;}
		        else if (text.find("qymin")!=std::string::npos) {str >> (*this).qymin; qymin *= pi_f;}				
		        else if (text.find("qymax")!=std::string::npos) {str >> (*this).qymax; qymax *= pi_f;}				
		        else if (text.find("qzmin")!=std::string::npos) {str >> (*this).qzmin; qzmin *= pi_f;}				
		        else if (text.find("qzmax")!=std::string::npos) {str >> (*this).qzmax; qzmax *= pi_f;}				
		        else if (text.find("nw")!=std::string::npos) str >> (*this).nw; 				
		        else if (text.find("wmin")!=std::string::npos) str >> (*this).wmin; 				
		        else if (text.find("wmax")!=std::string::npos) str >> (*this).wmax; 				
		        else if (text.find("scState")!=std::string::npos) str >> (*this).scState; 				
		        else if (text.find("printGap")!=std::string::npos) str >> (*this).printGap; 				
		        else if (text.find("gAmpl")!=std::string::npos) str >> (*this).gAmpl; 				
		        else if (text.find("Delta0")!=std::string::npos) str >> (*this).Delta0; 				
		        else if (text.find("deltax")!=std::string::npos) str >> (*this).deltax; 				
		        else if (text.find("deltay")!=std::string::npos) str >> (*this).deltay; 				
		        else if (text.find("deltaz")!=std::string::npos) str >> (*this).deltaz; 				
		        else if (text.find("kTrafo")!=std::string::npos) str >> (*this).kTrafo; 				
		        else if (text.find("options")!=std::string::npos) str >> (*this).options; 				
		        else if (text.find("subOptions")!=std::string::npos) str >> (*this).subOptions; 				
		        else if (text.find("Case")!=std::string::npos) str >> (*this).Case;
		        else if (text.find("pairingSpinParity")!=std::string::npos) str >> (*this).pairingSpinParity;
		        else if (text.find("pairingFromSpin")!=std::string::npos) str >> (*this).pairingFromSpin;
		        else if (text.find("pairingFromCharge")!=std::string::npos) str >> (*this).pairingFromCharge;
		        else if (text.find("storeChi")!=std::string::npos) str >> (*this).storeChi;
		        else if (text.find("storeGammaOrb")!=std::string::npos) str >> (*this).storeGammaOrb;
		        else if (text.find("readChiForPairing")!=std::string::npos) str >> (*this).readChi;
		        else if (text.find("readChiForSus")!=std::string::npos) str >> (*this).readChiForSus;
		        else if (text.find("hybridization")!=std::string::npos) str >> (*this).hyb;
		        else if (text.find("spinOrbit")!=std::string::npos) str >> (*this).LS;
		        else if (text.find("hybStrength")!=std::string::npos) str >> (*this).hybStrength;
		        else if (text.find("hybBand1")!=std::string::npos) str >> (*this).hybBand1;
		        else if (text.find("hybBand2")!=std::string::npos) str >> (*this).hybBand2;
		        else if (text.find("nSitesPerUnitCell")!=std::string::npos) str >> (*this).nSitesPerUnitCell;
		        else if (text.find("nOrbAtom")!=std::string::npos) str >> (*this).nOrbAtomStr;
		        else if (text.find("nkPerSheet")!=std::string::npos) str >> (*this).nkPerSheet;
		        else if (text.find("FSnkz")!=std::string::npos) str >> (*this).FSnkz;
		        else if (text.find("Omega0")!=std::string::npos) str >> (*this).Omega0;
		        else if (text.find("signF")!=std::string::npos) str >> (*this).signF;
		        else if (text.find("damp")!=std::string::npos) str >> (*this).damp;
		        else if (text.find("calcOnlyDiagonal")!=std::string::npos) str >> (*this).calcOnlyDiagonal;
		        else if (text.find("writeFullChi0")!=std::string::npos) str >> (*this).writeFullChi0;
		        else if (text.find("fixEvecs")!=std::string::npos) str >> (*this).fixEvecs;
			}

			void writeParameters(std::ostream& os) {
				os << "Case = " << (*this).Case << "\n";
				os << "dimension = " << (*this).dimension << "\n";
				os << "temperature = " << (*this).temperature << "\n";
				os << "numberOfOrbitals = " << (*this).nOrb << "\n";
				os << "chemicalPotential = " << (*this).mu << "\n";
				os << "tbParametersFile = " << (*this).tbfile << "\n";
				os << "complexHopping = " << (*this).complexHopping << "\n";
				os << "FSforPairingFile = " << (*this).fsfile << "\n";
				os << "pairingSpinParity = " << (*this).pairingSpinParity << "\n";
				os << "pairingFromSpin = " << (*this).pairingFromSpin << "\n";
				os << "pairingFromCharge = " << (*this).pairingFromCharge << "\n";
				os << "storeChi = " << (*this).storeChi << "\n";
				os << "storeGammaOrb = " << (*this).storeGammaOrb << "\n";
				os << "readChiForPairing = " << (*this).readChi << "\n";
				os << "readChiForSus = " << (*this).readChiForSus << "\n";
				os << "ChiforPairingFile = " << (*this).chifile << "\n";
				os << "interpolateChi = " << (*this).interpolateChi << "\n";
				os << "interpolateNqx = " << (*this).interpolateNqx << "\n";
				os << "interpolateNqz = " << (*this).interpolateNqz << "\n";
				os << "Coulomb1U = " << (*this).U << "\n";
				os << "Coulomb2Up = " << (*this).Up << "\n";
				os << "Coulomb3J = " << (*this).J << "\n";
				os << "Coulomb4Jp = " << (*this).Jp << "\n";
				os << "Coulomb5U_d = " << (*this).U_d_c << "\n";
				os << "Coulomb6U_p = " << (*this).U_p_c << "\n";
				os << "Coulomb7U_pd = " << (*this).U_pd_c << "\n";
				os << "Coulomb8U_pp = " << (*this).U_pp_c << "\n";
				os << "Coulomb9U_d = " << (*this).U_d_s << "\n";
				os << "Coulomb10U_p = " << (*this).U_p_s << "\n";
				os << "Coulomb11U_pd = " << (*this).U_pd_s << "\n";
				os << "Coulomb12U_pp = " << (*this).U_pp_s << "\n";
				os << "Coulomb13U_d_coupl = " << (*this).U_d_coupl << "\n";
				os << "Coulomb14U_p_coupl = " << (*this).U_p_coupl << "\n";
				os << "Coulomb15U_pd_coupl = " << (*this).U_pd_coupl << "\n";
				os << "Coulomb16U_pp_coupl = " << (*this).U_pp_coupl << "\n";
				os << "sublattice = " << (*this).sublattice << "\n";
				os << "deltaU0 = " << (*this).deltaU[0] << "\n";
				os << "deltaU1 = " << (*this).deltaU[1] << "\n";
				os << "deltaU2 = " << (*this).deltaU[2] << "\n";
				os << "deltaU3 = " << (*this).deltaU[3] << "\n";
				os << "deltaU4 = " << (*this).deltaU[4] << "\n";
				os << "staticUFactor = " << (*this).staticUFactor << "\n";
				os << "chargeFactor = " << (*this).chargeFactor << "\n";
				os << "spinFactor = " << (*this).spinFactor << "\n";
				os << "nkIntegration = " << (*this).nkInt << "\n";
				os << "nkzIntegration = " << (*this).nkIntz << "\n";
				os << "kz2D = " << (*this).kz2D << "\n";
				os << "a1x = " << (*this).a1[0] << "\n";
				os << "a1y = " << (*this).a1[1] << "\n";
				os << "a1z = " << (*this).a1[2] << "\n";
				os << "a2x = " << (*this).a2[0] << "\n";
				os << "a2y = " << (*this).a2[1] << "\n";
				os << "a2z = " << (*this).a2[2] << "\n";
				os << "a3x = " << (*this).a3[0] << "\n";
				os << "a3y = " << (*this).a3[1] << "\n";
				os << "a3z = " << (*this).a3[2] << "\n";
				os << "chia1x = " << (*this).chia1[0] << "\n";
				os << "chia1y = " << (*this).chia1[1] << "\n";
				os << "chia1z = " << (*this).chia1[2] << "\n";
				os << "chia2x = " << (*this).chia2[0] << "\n";
				os << "chia2y = " << (*this).chia2[1] << "\n";
				os << "chia2z = " << (*this).chia2[2] << "\n";
				os << "chia3x = " << (*this).chia3[0] << "\n";
				os << "chia3y = " << (*this).chia3[1] << "\n";
				os << "chia3z = " << (*this).chia3[2] << "\n";
				os << "nqx = " << (*this).nqx << "\n";
				os << "nqy = " << (*this).nqy << "\n";
				os << "nqz = " << (*this).nqz << "\n";
				os << "qxmin = " << (*this).qxmin << "\n";
				os << "qxmax = " << (*this).qxmax << "\n";
				os << "qymin = " << (*this).qymin << "\n";
				os << "qymax = " << (*this).qymax << "\n";
				os << "qzmin = " << (*this).qzmin << "\n";
				os << "qzmax = " << (*this).qzmax << "\n";
				os << "nw = " << (*this).nw << "\n";
				os << "wmin = " << (*this).wmin << "\n";
				os << "wmax = " << (*this).wmax << "\n";
				os << "scState = " << (*this).scState << "\n";
				os << "printGap = " << (*this).printGap << "\n";
				os << "gAmpl = " << (*this).gAmpl << "\n";
				os << "Delta0 = " << (*this).Delta0 << "\n";
				os << "deltax = " << (*this).deltax << "\n";
				os << "deltay = " << (*this).deltay << "\n";
				os << "deltaz = " << (*this).deltaz << "\n";
				os << "kTrafo = " << (*this).kTrafo << "\n";
				os << "options = " << (*this).options << "\n";
				os << "subOptions = " << (*this).subOptions << "\n";
				os << "hybridization = " << (*this).hyb << "\n";
				os << "spinOrbit = " << (*this).LS << "\n";
				os << "hybStrength = " << (*this).hybStrength << "\n";
				os << "hybBand1 = " << (*this).hybBand1 << "\n";
				os << "hybBand2 = " << (*this).hybBand2 << "\n";
				os << "nSitesPerUnitCell = " << (*this).nSitesPerUnitCell << "\n";
				os << "nOrbAtom = ";
				for (int i=0; i<nOrbAtom.size();i++) os << (*this).nOrbAtom[i] << " "; 
				os << "\n";
				os << "nkPerSheet = " << (*this).nkPerSheet << "\n";
				os << "FSnkz = " << (*this).FSnkz << "\n";
				os << "Omega0 = " << (*this).Omega0 << "\n";
				os << "signF = " << (*this).signF << "\n";
				os << "damp = " << (*this).damp << "\n";
				os << "calcOnlyDiagonal = " << (*this).calcOnlyDiagonal << "\n";
				os << "writeFullChi0 = " << (*this).writeFullChi0 << "\n";
				os << "fixEvecs = " << (*this).fixEvecs << "\n";
			}


			void bcTest() {
				int h(0);
				if (conc.rank()==0) h=5;
		        conc.broadcast(h); 
			}
			void broadcastParam() {

		        if (conc.rank()==0) std::cout << "Now broadcasting parameters! \n";
		        conc.broadcast((*this).Case);
		        conc.broadcast((*this).dimension); 
		        conc.broadcast((*this).temperature); 
		        conc.broadcast((*this).nOrb); 
		        conc.broadcast((*this).mu); 
		        conc.broadcast((*this).tbfile); 
		        conc.broadcast((*this).complexHopping);
		        conc.broadcast((*this).fsfile); 
		        conc.broadcast((*this).chifile); 
		        conc.broadcast((*this).interpolateChi); 
		        conc.broadcast((*this).interpolateNqx); 
		        conc.broadcast((*this).interpolateNqz); 
		        conc.broadcast((*this).U); 
		        conc.broadcast((*this).Up); 
		        conc.broadcast((*this).J); 
		        conc.broadcast((*this).Jp); 
		        conc.broadcast((*this).U_d_c); 
		        conc.broadcast((*this).U_p_c); 
		        conc.broadcast((*this).U_pd_c); 
		        conc.broadcast((*this).U_pp_c); 
		        conc.broadcast((*this).U_d_s); 
		        conc.broadcast((*this).U_p_s); 
		        conc.broadcast((*this).U_pd_s); 
		        conc.broadcast((*this).U_pp_s); 
		        conc.broadcast((*this).U_d_coupl); 
		        conc.broadcast((*this).U_p_coupl); 
		        conc.broadcast((*this).U_pd_coupl); 
		        conc.broadcast((*this).U_pp_coupl); 
		        conc.broadcast((*this).sublattice); 
		        conc.broadcast((*this).deltaU[0]); 
		        conc.broadcast((*this).deltaU[1]); 
		        conc.broadcast((*this).deltaU[2]); 
		        conc.broadcast((*this).deltaU[3]); 
		        conc.broadcast((*this).deltaU[4]); 
		        conc.broadcast((*this).staticUFactor); 				
		        conc.broadcast((*this).chargeFactor); 				
		        conc.broadcast((*this).spinFactor); 				
		        conc.broadcast((*this).nkInt); 				
		        conc.broadcast((*this).nkIntz); 				
		        conc.broadcast((*this).kz2D); 				
		        conc.broadcast((*this).a1[0]); 				
		        conc.broadcast((*this).a1[1]); 				
		        conc.broadcast((*this).a1[2]); 				
		        conc.broadcast((*this).a2[0]); 				
		        conc.broadcast((*this).a2[1]); 				
		        conc.broadcast((*this).a2[2]); 				
		        conc.broadcast((*this).a3[0]); 				
		        conc.broadcast((*this).a3[1]); 				
		        conc.broadcast((*this).a3[2]); 				
		        conc.broadcast((*this).chia1[0]); 				
		        conc.broadcast((*this).chia1[1]); 				
		        conc.broadcast((*this).chia1[2]); 				
		        conc.broadcast((*this).chia2[0]); 				
		        conc.broadcast((*this).chia2[1]); 				
		        conc.broadcast((*this).chia2[2]); 				
		        conc.broadcast((*this).chia3[0]); 				
		        conc.broadcast((*this).chia3[1]); 				
		        conc.broadcast((*this).chia3[2]); 				
		        conc.broadcast((*this).nqx); 				
		        conc.broadcast((*this).nqy); 				
		        conc.broadcast((*this).nqz); 				
		        conc.broadcast((*this).qxmin);
		        conc.broadcast((*this).qxmax);				
		        conc.broadcast((*this).qymin);				
		        conc.broadcast((*this).qymax);				
		        conc.broadcast((*this).qzmin);				
		        conc.broadcast((*this).qzmax);				
		        conc.broadcast((*this).nw); 				
		        conc.broadcast((*this).wmin); 				
		        conc.broadcast((*this).wmax); 				
		        conc.broadcast((*this).scState); 				
		        conc.broadcast((*this).printGap); 				
		        conc.broadcast((*this).gAmpl); 				
		        conc.broadcast((*this).Delta0); 				
		        conc.broadcast((*this).deltax); 			
		        conc.broadcast((*this).deltay); 				
		        conc.broadcast((*this).deltaz); 				
		        conc.broadcast((*this).kTrafo); 				
		        conc.broadcast((*this).options); 				
		        conc.broadcast((*this).subOptions); 				
		        conc.broadcast((*this).pairingSpinParity);
		        conc.broadcast((*this).pairingFromSpin);
		        conc.broadcast((*this).pairingFromCharge);
		        conc.broadcast((*this).storeChi);
		        conc.broadcast((*this).storeGammaOrb);
		        conc.broadcast((*this).readChi);
		        conc.broadcast((*this).readChiForSus);
		        conc.broadcast((*this).hyb);
		        conc.broadcast((*this).LS);
		        conc.broadcast((*this).hybStrength);
		        conc.broadcast((*this).hybBand1);
		        conc.broadcast((*this).hybBand2);
		        conc.broadcast((*this).nSitesPerUnitCell);
		        conc.broadcast((*this).nOrbAtomStr);
		        conc.broadcast((*this).nkPerSheet);
		        conc.broadcast((*this).FSnkz);
		        conc.broadcast((*this).Omega0);
		        conc.broadcast((*this).signF);
		        conc.broadcast((*this).damp);
		        conc.broadcast((*this).writeFullChi0);
		        conc.broadcast((*this).fixEvecs);
			}

			void setupOrbitalIndices(){
				indexToOrb.resize(nOrb*nOrb,2);
				for (size_t l1 = 0; l1 < nOrb; ++l1){
					for (size_t l2 = 0; l2 < nOrb; ++l2){
						size_t ind=l2+l1*nOrb;
						indexToOrb(ind,0) = l1;
						indexToOrb(ind,1) = l2;
					}
				}
				for (size_t l=0;l<nOrb; l++) {
					int ll(l);
					for (size_t site=0;site<nSitesPerUnitCell;site++) {
						ll -= nOrbAtom[site];
						if (ll < 0) {
							orbToSite.push_back(site);
							break;
						}
					}
				}
				if (conc.rank()==0) {
					std::cout << "orbToSite: ";
					for (size_t l=0;l<nOrb;l++) std::cout << l << "->" << orbToSite[l] << " ";
					std::cout << "\n";
				}
			}


			template<typename FieldType>
			void loadVector(std::vector<FieldType>& v,const std::string& vstring)
			{
				std::stringstream ss(vstring);
				int value;
				while (ss >> value)
				{
					v.push_back(value);
					if (ss.peek() == ',') ss.ignore();
				}
			}


	};

}

#endif
