//-*-C++-*-

#ifndef FERMINATOR_H
#define FERMINATOR_H

#include <string>
#include <vector>
#include <fstream>
#include <complex>
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "utilities.h"
#include "bandstructure.h"

namespace rpa {

	template<typename Field, typename BandsType,
	         template<typename> class MatrixTemplate,
	         typename ConcurrencyType>
	class ferminator {


	private:
		typedef Field 							FieldType;
		typedef std::complex<Field>				ComplexType;
		typedef MatrixTemplate<Field> 			MatrixType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      		VectorType;
		typedef std::vector<int>      		    VectorIntType;
		typedef std::vector<ComplexType>      	ComplexVectorType;


		rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		ConcurrencyType& conc;
		FieldType Pi;
		std::string Case_;

	public:
		size_t nkPerSheet,nSheets,nTotal;
		VectorType kFx,kFy,kFz;
		VectorIntType kFtoBand;
		VectorType deltakF,vkF,gammaB1GkF;
		std::vector<std::vector<FieldType> > FSCenters;
		std::vector<size_t> FSBand;
		BandsType bands;
		bool calcOW_;
		std::vector<std::vector<ComplexType> > weights;
		std::vector<std::vector<ComplexType> > weights2;


	ferminator(rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
			   ConcurrencyType& concurrency,bool calcOW=0):
		param(parameters),
		conc(concurrency),
		Pi(param.pi_f),
		Case_(param.Case),
		nkPerSheet(param.nkPerSheet),
		nSheets(0),
		nTotal(0),
		kFx(0,0),
		kFy(0,0),
		kFz(0,0),
		kFtoBand(0,0),
		deltakF(0,0),
		vkF(0,0),
		gammaB1GkF(0,0),
		FSCenters(0),
		FSBand(0),
		bands(param,conc),
		calcOW_(calcOW),
		weights(0),
		weights2(0)


		{
			// if (param.readFSFromFile) readFromFile(); // if data has additional columns for deltakf and vkf
			if (param.readFSFromFile) {
				readFromFile2();// if data only has kFx, kFy, kFz, band
				if (conc.rank()==0) writeKF();
			} else {
				if (conc.rank()==0) std::cout << "Now setting up Fermi surface for case " << Case_ << "\n";
				setFSValues();
				if (conc.rank()==0) writeKF();
				if (conc.rank()==0) std::cout << "Done setting up Fermi surface \n";
			}
		}

	void setFSValues() {

		if (Case_ == "KFeSe_10Orbit_underdoped_2D") {
			// 4 FS sheets total, two around (pi,pi), two around (-pi,pi)
			nSheets = 4;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  6;
			FSCenters[2][0] =-Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =-Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;

			size_t nkSearch(128);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);

		} else if (Case_ == "2D_5orbit_only_electron_pockets") {
			// 2 FS sheets total, 1 around (pi,0), one around (0,pi)
			nSheets = 2;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = 0;
			FSBand   [0]    =  3;
			FSCenters[1][0] = 0; FSCenters[1][1] = Pi;
			FSBand   [1]    =  3;

			size_t nkSearch(128);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);


		} else if (Case_ == "2D_4orbit_2_hole_2_electron") {
			// 4-orbital model: 4 FS sheets total, 1 around (pi,0), one around (0,pi), 2 at Gamma
			nSheets = 4;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0; FSCenters[0][1] = 0;
			FSBand   [0]    =  0;
			FSCenters[1][0] = 0; FSCenters[1][1] = 0;
			FSBand   [1]    =  1;
			FSCenters[2][0] = Pi; FSCenters[2][1] = 0;
			FSBand   [2]    =  2;
			FSCenters[3][0] = 0; FSCenters[3][1] = Pi;
			FSBand   [3]    =  2;

			size_t nkSearch(128);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);


		} else if (Case_ == "2D_5orbit_5Sheets") {
			// 5 FS sheets total, 2 at (0,0), 1 at (pi,0), one at (0,pi), 1 at (pi,pi)
			nSheets = 5;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] =  0; FSCenters[0][1] = 0;
			FSBand   [0]    =  1;
			FSCenters[1][0] =  0; FSCenters[1][1] = 0;
			FSBand   [1]    =  2;
			FSCenters[2][0] = Pi; FSCenters[2][1] = 0;
			FSBand   [2]    =  3;
			FSCenters[3][0] = 0; FSCenters[3][1] = Pi;
			FSBand   [3]    =  3;
			FSCenters[4][0] = Pi; FSCenters[4][1] = Pi;
			FSBand   [4]    =  2;

			size_t nkSearch(128);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);

		} else if (Case_ == "2D_5orbit_4Sheets") {
			// 4 FS sheets total, 2 at (0,0), 1 at (pi,0), one at (0,pi)
			nSheets = 4;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] =  0; FSCenters[0][1] = 0;
			FSBand   [0]    =  1;
			FSCenters[1][0] =  0; FSCenters[1][1] = 0;
			FSBand   [1]    =  2;
			FSCenters[2][0] = Pi; FSCenters[2][1] = 0;
			FSBand   [2]    =  3;
			FSCenters[3][0] = 0; FSCenters[3][1] = Pi;
			FSBand   [3]    =  3;

			size_t nkSearch(128);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);

		} else if (Case_ == "KFeSe_10Orbit_overdoped_2D") {
			// 4 FS sheets total, two around (pi,pi), two around (-pi,pi)
			nSheets = 5;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  6;
			FSCenters[2][0] =-Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =-Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;
			FSCenters[4][0] =2*Pi; FSCenters[4][1] = 0.0;
			FSBand   [4]    =  6;

			size_t nkSearch(128);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);

		} else if (Case_ == "KFeSe_10Orbit_underdoped_3D") {
			// 4 FS sheets total, two around (pi,pi), two around (-pi,pi)
			nSheets = 4;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  6;
			FSCenters[2][0] =-Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =-Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(float(ikz)*4.*param.pi_f/float(param.FSnkz)-2.*param.pi_f);
					calcKF(nkSearch,iSheet,kz,3);
				}
			}
		} else if (Case_ == "10Orbit_3D_10_01") {
			// 4 FS sheets total, two around (pi,pi), two around (-pi,pi)
			if (conc.rank()==0) std::cout << "Case: " << "10Orbit_3D_10_01 \n";
			nSheets = 4;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = 0;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = 0;
			FSBand   [1]    =  6;
			FSCenters[2][0] =  0; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =  0; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					if (conc.rank()==0) std::cout << "iSheet,ikz="<<iSheet<<","<<ikz<<"\n";
					FieldType kz(float(ikz)*2.*param.pi_f/float(param.FSnkz)-1.*param.pi_f);
					calcKF(nkSearch,iSheet,kz,3);
				}
			}
		} else if (Case_ == "KFeSe_10Orbit_overdoped_3D") {
			// 5 FS sheets total, two around (pi,pi), two around (-pi,pi), one (Z) around (0,0,2pi)
			// Note that the Z pocket is closed along kz

			nSheets = 5;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  6;
			FSCenters[2][0] =-Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =-Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;
			FSCenters[4][0] =2*Pi;FSCenters[4][1] = 0.0;
			FSBand   [4]    =  6;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(float(ikz)*4.*param.pi_f/float(param.FSnkz)-2.*param.pi_f);
					calcKF(nkSearch,iSheet,kz,3);
				}
			}

		} else if (Case_ == "BaFeAs_10Orbit_3D_6Sheets") {
			// 6 FS sheets total, two around (pi,pi), two around (-pi,pi), two around (0,0,0)
			// Note: One of the hole-pockets around Gamma may be closed along kz

			nSheets = 6;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  6;
			FSCenters[2][0] =-Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =-Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;
			FSCenters[4][0] =0.0;FSCenters[4][1] = 0.0;
			FSBand   [4]    =  4;
			FSCenters[5][0] =0.0;FSCenters[5][1] = 0.0;
			FSBand   [5]    =  5;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(float(ikz)*4.*param.pi_f/float(param.FSnkz)-2.*param.pi_f);
					calcKF(nkSearch,iSheet,kz,3);
				}
			}

		} else if (Case_ == "BaFeAs_10Orbit_3D_7Sheets" || Case_ == "1111_10orbit") {
			// 7 FS sheets total, two around (pi,pi), two around (-pi,pi), 3 around (0,0,0)
			// Note: One of the hole-pockets around Gamma may be closed along kz

			nSheets = 7;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  6;
			FSCenters[2][0] =-Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =-Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;
			FSCenters[4][0] =0.0;FSCenters[4][1] = 0.0;
			FSBand   [4]    =  4;
			FSCenters[5][0] =0.0;FSCenters[5][1] = 0.0;
			FSBand   [5]    =  5;
			FSCenters[6][0] =0.0;FSCenters[6][1] = 0.0;
			FSBand   [6]    =  3;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				FieldType kz(param.kz2D);
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					if (param.FSnkz > 1)  kz = float(ikz)*4.*param.pi_f/float(param.FSnkz)-2.*param.pi_f;
					calcKF(nkSearch,iSheet,kz,3);
				}
			}

		} else if (Case_ == "10Orbit_3D_7Sheets_00_10_01" || Case_ == "1111_10orbit") {
			// 7 FS sheets total, two around (pi,pi), two around (-pi,pi), 3 around (0,0,0)
			// Note: One of the hole-pockets around Gamma may be closed along kz

			nSheets = 7;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = 0;
			FSBand   [0]    =  7;
			FSCenters[1][0] = Pi; FSCenters[1][1] = 0;
			FSBand   [1]    =  6;
			FSCenters[2][0] =  0; FSCenters[2][1] = Pi;
			FSBand   [2]    =  7;
			FSCenters[3][0] =  0; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;
			FSCenters[4][0] =0.0;FSCenters[4][1] = 0.0;
			FSBand   [4]    =  4;
			FSCenters[5][0] =0.0;FSCenters[5][1] = 0.0;
			FSBand   [5]    =  5;
			FSCenters[6][0] =0.0;FSCenters[6][1] = 0.0;
			FSBand   [6]    =  3;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				FieldType kz(param.kz2D);
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					if (param.FSnkz > 1)  kz = float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f;
					calcKF(nkSearch,iSheet,kz,3);
				}
			}

		}  else if (Case_ == "5Orbit_3D_6Sheets_00_11_10_01" || Case_ == "1111_10orbit") {
			// 7 FS sheets total, two around (pi,pi), two around (-pi,pi), 3 around (0,0,0)
			// Note: One of the hole-pockets around Gamma may be closed along kz

			nSheets = 6;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] =  0; FSCenters[0][1] = 0;
			FSBand   [0]    =  1;
			FSCenters[1][0] =  0; FSCenters[1][1] = 0;
			FSBand   [1]    =  2;
			FSCenters[2][0] =  Pi; FSCenters[2][1] = Pi;
			FSBand   [2]    =  1;
			FSCenters[3][0] =  Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  2;
			FSCenters[4][0] = Pi;FSCenters[4][1] = 0.0;
			FSBand   [4]    =  3;
			FSCenters[5][0] =  0;FSCenters[5][1] = Pi;
			FSBand   [5]    =  3;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				FieldType kz(param.kz2D);
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					if (param.FSnkz > 1)  kz = float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f;
					calcKF(nkSearch,iSheet,kz,3);
				}
			}

		} else if (Case_ == "BaFeAs_5Orbit_5Sheets") {
			// 5 FS sheets total, two around (0,0), one each around (pi,0) and (0,pi), 1 around (pi,pi)

			nSheets = 5;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			std::cout << nTotal << " k-points expected \n"; 
			resizeContainers();

			FSCenters[0][0] = 0.0; FSCenters[0][1] = 0.0;
			FSBand   [0]    =   1;
			FSCenters[1][0] = 0.0; FSCenters[1][1] = 0.0;
			FSBand   [1]    =   2;
			FSCenters[2][0] =  Pi; FSCenters[2][1] = 0.0;
			FSBand   [2]    =   3;
			FSCenters[3][0] = 0.0; FSCenters[3][1] = Pi;
			FSBand   [3]    =   3;
			FSCenters[4][0] =  Pi;FSCenters[4][1] =  Pi;
			FSBand   [4]    =   2;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				// for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					// FieldType kz(float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f);
					FieldType kz(0.0);
					calcKF(nkSearch,iSheet,kz,3);
				// }
			}
		} else if (Case_ == "BaFeAs_5Orbit_4Sheets") {
			// 4 FS sheets total, two around (0,0), one each around (pi,0) and (0,pi)

			nSheets = 4;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = 0.0; FSCenters[0][1] = 0.0;
			FSBand   [0]    =   1;
			FSCenters[1][0] = 0.0; FSCenters[1][1] = 0.0;
			FSBand   [1]    =   2;
			FSCenters[2][0] =  Pi; FSCenters[2][1] = 0.0;
			FSBand   [2]    =   3;
			FSCenters[3][0] = 0.0; FSCenters[3][1] = Pi;
			FSBand   [3]    =   3;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(param.kz2D);
					if (param.FSnkz > 1) kz = float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f;
					calcKF(nkSearch,iSheet,kz,3);
				}
			}

		} else if (Case_ == "3D_5Orbit_electron_sheets") {
			// 2 FS sheets total, one each around (pi,0) and (0,pi)

			nSheets = 2;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] =  Pi; FSCenters[0][1] = 0.0;
			FSBand   [0]    =   3;
			FSCenters[1][0] = 0.0; FSCenters[1][1] = Pi;
			FSBand   [1]    =   3;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(param.kz2D);
					if (param.FSnkz > 1) kz = float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f;
					calcKF(nkSearch,iSheet,kz,3);
				}
			}
		} else if (Case_ == "BSCCObilayer_OD_1band") {
			// 1 FS sheets total with two different kz, one around (pi,pi,0), one around (pi,pi,pi)

			nSheets = 1;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] =  Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =   0;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(param.kz2D);
					if (param.FSnkz > 1) kz = float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f;
					calcKF(nkSearch,iSheet,kz,3);
				}
			}
		} else if (Case_ == "BSCCObilayer_OD_1band_onlyA") {
			// 1 FS sheets total with only one kz given by param.kz2D

			nSheets = 1;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] =  Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =   0;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				FieldType kz(param.kz2D);
				if (conc.rank()==0) std::cout << "FS for kz = " << kz << "\n";
				calcKF(nkSearch,iSheet,kz,3);
			}
		} else if (Case_ == "KFeAs_10Orbit_2D") {
			// 3 FS sheets total around (0,0), neglecting tiny electron pockets near (pi,pi)
			nSheets = 3;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0.; FSCenters[0][1] = 0.;
			FSBand   [0]    =  3;
			FSCenters[1][0] = 0.; FSCenters[1][1] = 0.;
			FSBand   [1]    =  4;
			FSCenters[2][0] = 0.; FSCenters[2][1] = 0.;
			FSBand   [2]    =  5;

			size_t nkSearch(512);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);
		}   else if (Case_ == "FeSe_monolayer") {
			// 5 FS sheets total, 3 hole around (0,0) + 2 electron pockets
			nSheets = 7;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0.; FSCenters[0][1] = 0.;
			FSBand   [0]    =  3;
			FSCenters[1][0] = 0.; FSCenters[1][1] = 0.;
			FSBand   [1]    =  4;
			FSCenters[2][0] = 0.; FSCenters[2][1] = 0.;
			FSBand   [2]    =  5;
			FSCenters[3][0] = Pi; FSCenters[3][1] = Pi;
			FSBand   [3]    =  6;
			FSCenters[4][0] = Pi; FSCenters[4][1] = Pi;
			FSBand   [4]    =  7;
			FSCenters[5][0] = -Pi; FSCenters[5][1] = Pi;
			FSBand   [5]    =  6;
			FSCenters[6][0] = -Pi; FSCenters[6][1] = Pi;
			FSBand   [6]    =  7;

			size_t nkSearch(512);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,param.kz2D,2);
		} else if (Case_ == "YBCO_orthoII_perpStripes") {
			// 1-band Hubbard bilayer with perp. running Stripess; 8 bands total. 24 FS sheets total.
			nSheets = 24;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi/2.  ; FSCenters[0][1] = Pi/2.  ; FSBand[0]  =  0;
			FSCenters[1][0] = Pi/2.  ; FSCenters[1][1] =-Pi/2.  ; FSBand[1]  =  0;
			FSCenters[2][0] =-Pi/2.  ; FSCenters[2][1] = Pi/2.  ; FSBand[2]  =  0;
			FSCenters[3][0] =-Pi/2.  ; FSCenters[3][1] =-Pi/2.  ; FSBand[3]  =  0;

			FSCenters[4][0] = Pi/2.  ; FSCenters[4][1] = Pi/2.  ; FSBand[4]  =  1;
			FSCenters[5][0] = Pi/2.  ; FSCenters[5][1] =-Pi/2.  ; FSBand[5]  =  1;
			FSCenters[6][0] =-Pi/2.  ; FSCenters[6][1] = Pi/2.  ; FSBand[6]  =  1;
			FSCenters[7][0] =-Pi/2.  ; FSCenters[7][1] =-Pi/2.  ; FSBand[7]  =  1;

			FSCenters[8][0] = Pi/2.  ; FSCenters[8][1] = Pi/2.  ; FSBand[8]  =  2;
			FSCenters[9][0] = Pi/2.  ; FSCenters[9][1] =-Pi/2.  ; FSBand[9]  =  2;
			FSCenters[10][0] =-Pi/2. ; FSCenters[10][1] = Pi/2. ; FSBand[10] =  2;
			FSCenters[11][0] =-Pi/2. ; FSCenters[11][1] =-Pi/2. ; FSBand[11] =  2;

			FSCenters[12][0] = Pi/2. ; FSCenters[12][1] = Pi/2. ; FSBand[12] =  3;
			FSCenters[13][0] = Pi/2. ; FSCenters[13][1] =-Pi/2. ; FSBand[13] =  3;
			FSCenters[14][0] =-Pi/2. ; FSCenters[14][1] = Pi/2. ; FSBand[14] =  3;
			FSCenters[15][0] =-Pi/2. ; FSCenters[15][1] =-Pi/2. ; FSBand[15] =  3;

			FSCenters[16][0] = 0. ; FSCenters[16][1] = 0. ; FSBand[16] =  4;
			FSCenters[17][0] = Pi ; FSCenters[17][1] = 0. ; FSBand[17] =  4;
			FSCenters[18][0] = 0. ; FSCenters[18][1] = Pi ; FSBand[18] =  4;
			FSCenters[19][0] = Pi ; FSCenters[19][1] = Pi ; FSBand[19] =  4;

			FSCenters[20][0] = 0. ; FSCenters[20][1] = 0. ; FSBand[20] =  5;
			FSCenters[21][0] = Pi ; FSCenters[21][1] = 0. ; FSBand[21] =  5;
			FSCenters[22][0] = 0. ; FSCenters[22][1] = Pi ; FSBand[22] =  5;
			FSCenters[23][0] = Pi ; FSCenters[23][1] = Pi ; FSBand[23] =  5;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "bilayer") {
			// 1-band Hubbard bilayer; 2 bands total. 2 FS sheets. Here we assume that they are both closed around (pi,pi)
			nSheets = 2;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  0;
			FSCenters[1][0] = Pi  ; FSCenters[1][1] = Pi  ; FSBand[1]  =  1;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "BSCCObilayer_OD") {
			// 1-band Hubbard bilayer; 2 bands total. 2 FS sheets, both closed around (pi,pi)
			nSheets = 2;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  0;
			FSCenters[1][0] = Pi  ; FSCenters[1][1] = Pi  ; FSBand[1]  =  1;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "BSCCObilayer_beyOD") {
			// beyond overdoped --> 
			// 1-band Hubbard bilayer; 2 bands total. 2 FS sheets, 1 closed around (pi,pi)
			// 1 closed around Gamma
			nSheets = 2;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  0;
			FSCenters[1][0] = 0   ; FSCenters[1][1] = 0   ; FSBand[1]  =  1;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "Raghu") {
			// 2 bands, 4 FS sheets, 1 around Gamma, 1 around (pi,pi), 1 around X, 1 around Y 
			nSheets = 4;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0   ; FSCenters[0][1] = 0   ; FSBand[0]  =  0;
			FSCenters[1][0] = Pi  ; FSCenters[1][1] = Pi  ; FSBand[1]  =  0;
			FSCenters[2][0] = Pi  ; FSCenters[2][1] = 0   ; FSBand[2]  =  1;
			FSCenters[3][0] = 0   ; FSCenters[3][1] = Pi  ; FSBand[3]  =  1;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "Ba2CuO3") {
			// 2 bands, 2 FS sheets, 1 around Gamma, 1 around (pi,pi)
			nSheets = 2;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0   ; FSCenters[0][1] = 0   ; FSBand[0]  =  1;
			FSCenters[1][0] = Pi  ; FSCenters[1][1] = Pi  ; FSBand[1]  =  0;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "Ba2CuO3_1sheet") {
			// 2 bands, 2 FS sheets, 1 around Gamma, 1 around (pi,pi)
			nSheets = 1;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi   ; FSCenters[0][1] = Pi   ; FSBand[0]  =  1;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

		} else if (Case_ == "1band_M") {
			// 2D 1-band Hubbard; 1 band, 1 FS sheet. Here we assume that it is closed around (pi,pi)
			nSheets = 1;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  0;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);
		} else if (Case_ == "1band_Gamma") {
			// 2D 1-band Hubbard; 1 band, 1 FS sheet. Here we assume that it is closed around (0,0)
			nSheets = 1;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0  ; FSCenters[0][1] = 0  ; FSBand[0]  =  0;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);
		}  else if (Case_ == "Emery" || Case_ == "EmeryOnlyUd") {
			// 2D 3-band Hubbard; 3 bands, 1 FS sheet from band 2.
			nSheets = 1;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			if (param.mu > 3.609) {
				FSCenters[0][0] = Pi  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  2;
			} else {
				FSCenters[0][0] = 0  ; FSCenters[0][1] = 0  ; FSBand[0]  =  2;
			}

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);

			// Re-calculate deltakF
			calcDeltaKF(kFx,kFy,deltakF);

		} else if (Case_ == "1band_Gamma") {
			// 2D 1-band Hubbard; 1 band, 1 FS sheet. Here we assume that it is closed around (0,0)
			nSheets = 1;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0  ; FSCenters[0][1] = 0  ; FSBand[0]  =  0;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);
		} else if (Case_ == "Checkerboard") {
			// 2D Checkerboard Hubbard; 4 bands, 2 FS sheet. Here we assume that it is closed around (pi,pi)
			nSheets = 8;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = 0.5*Pi	;	FSCenters[0][1] = 0			; FSBand[0]  =  1;
			FSCenters[1][0] = 0   		;	FSCenters[1][1] = 0.5*Pi	; FSBand[1]  =  1;
			FSCenters[2][0] =-0.5*Pi	;	FSCenters[2][1] = 0			; FSBand[2]  =  1;
			FSCenters[3][0] = 0   		;	FSCenters[3][1] =-0.5*Pi	; FSBand[3]  =  1;
			FSCenters[4][0] =-0.5*Pi	;	FSCenters[4][1] = Pi		; FSBand[4]  =  1;
			FSCenters[5][0] = 0.5*Pi  	;	FSCenters[5][1] = Pi		; FSBand[5]  =  1;
			FSCenters[6][0] =  Pi		;	FSCenters[6][1] = 0.5*Pi	; FSBand[6]  =  1;
			FSCenters[7][0] =  Pi   	;	FSCenters[7][1] =-0.5*Pi	; FSBand[7]  =  1;

			size_t nkSearch(256);
			for (size_t iSheet=0;iSheet<nSheets;iSheet++) calcKF(nkSearch,iSheet,0.0,2);
		} else if (Case_ == "Cuprate_Bilayer_3D") {
			// 2 FS sheets total both cylidrical around (pi,pi)

			nSheets = 2;
			nTotal = nSheets * nkPerSheet * param.FSnkz;
			resizeContainers();

			FSCenters[0][0] = Pi; FSCenters[0][1] = Pi;
			FSBand   [0]    =  0;
			FSCenters[1][0] = Pi; FSCenters[1][1] = Pi;
			FSBand   [1]    =  1;

			size_t nkSearch(256);

			for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
				if (conc.rank()==0) std::cout << "Sheet nr. " << iSheet << "\n";
				for (size_t ikz=0;ikz<param.FSnkz;ikz++) {
					FieldType kz(float(ikz)*2.*param.pi_f/float(param.FSnkz)-param.pi_f);
					calcKF(nkSearch,iSheet,kz,3);
				}
			}
		} else if (Case_ == "trilayer") {
			// 2D 1-band trilayer; 3 bands, 4 FS sheets.
			nSheets = 4;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = Pi  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  0;
			FSCenters[1][0] = Pi  ; FSCenters[1][1] = Pi  ; FSBand[1]  =  1;
			FSCenters[2][0] =-Pi  ; FSCenters[2][1] = 0.  ; FSBand[2]  =  2;

			size_t nkSearch(10000);
			for (size_t iSheet=0;iSheet<2;iSheet++) calcKF(nkSearch,iSheet,0.0,2);
			for (size_t iSheet=2;iSheet<3;iSheet++) calcKFOpen2D(nkSearch,iSheet,1);
			// Now copy isheet = 2 to isheet=3 but with -kFy to preserve symmetry
			size_t ind0(kFx.size()-nkPerSheet);
			for (size_t ikF=0; ikF<nkPerSheet; ikF++) {
				kFx.push_back(kFx[ind0+ikF]);
				kFy.push_back(-kFy[ind0+ikF]);
				kFz.push_back(kFz[ind0+ikF]);
				kFtoBand.push_back(kFtoBand[ind0+ikF]);
				deltakF.push_back(deltakF[ind0+ikF]);
				vkF.push_back(vkF[ind0+ikF]);
			}

		} else if (Case_ == "coupled_ladders") {
			// 1-band model, coupled ladders along x
			nSheets = 4;
			nTotal = nSheets * nkPerSheet;
			resizeContainers();

			FSCenters[0][0] = -Pi/2.  ; FSCenters[0][1] = Pi  ; FSBand[0]  =  0;
			FSCenters[1][0] =  Pi/2.  ; FSCenters[1][1] = Pi  ; FSBand[1]  =  0;
			FSCenters[2][0] =  -Pi    ; FSCenters[2][1] = 0.  ; FSBand[2]  =  1;

			size_t nkSearch(10000);
			for (size_t iSheet=0;iSheet<2;iSheet++) calcKF(nkSearch,iSheet,0.0,2);
			for (size_t iSheet=2;iSheet<3;iSheet++) calcKFOpen2D(nkSearch,iSheet,1);
			// // Now copy isheet = 2 to isheet=3 but with -kFy to preserve symmetry
			size_t ind0(kFx.size()-nkPerSheet);
			for (size_t ikF=0; ikF<nkPerSheet; ikF++) {
				kFx.push_back(kFx[ind0+ikF]);
				kFy.push_back(-kFy[ind0+ikF]);
				kFz.push_back(kFz[ind0+ikF]);
				kFtoBand.push_back(kFtoBand[ind0+ikF]);
				deltakF.push_back(deltakF[ind0+ikF]);
				vkF.push_back(vkF[ind0+ikF]);
				if (calcOW_) {
					size_t ic(kFx.size()-1);
					weights[ic]  = weights[ind0+ikF];
					weights2[ic] = weights2[ind0+ikF];
				}
			}
		}


		if (conc.rank()==0) std::cout << "Number of kF points: " << nTotal << "\n";

		if (kFx.size() != nTotal) {
			if (conc.rank()==0) std::cout << kFx.size() << " kF points instead of " << nTotal << ". Perhaps closed FS?\n";
			nTotal = kFx.size();
		}
	}

	void resizeContainers() {
			FSCenters.resize(nSheets,VectorType(2,0));
			FSBand.resize(nSheets);
			if (calcOW_) {
				weights. resize(nTotal,std::vector<ComplexType>(param.nOrb,0));
				weights2.resize(nTotal,std::vector<ComplexType>(param.nOrb,0));
			}
	}

	void calcKF(const size_t nkSearch,const size_t iSheet, const FieldType& kz,int dim) {
		VectorType w(param.nOrb);
		ComplexMatrixType v(param.nOrb,param.nOrb);
		// for (size_t iSheet=0;iSheet<nSheets;iSheet++) {
			for (size_t ik=0;ik<nkPerSheet;ik++) {
				FieldType Theta = 2.*Pi/float(nkPerSheet)*float(ik);
				size_t ikF = 0;
				VectorType k(3,0);
				for (size_t il=0;il<2;il++) k[il] = FSCenters[iSheet][il];
				k[2] = kz;
				bands.getBands(k,w,v);
				FieldType ek(w[FSBand[iSheet]]);
				FieldType sgnW = (ek > 0) - (ek < 0);
				FieldType sgnOld(sgnW);
				VectorType kOld(k);
				while (ikF <= nkSearch && sgnW == sgnOld) {
					sgnOld = sgnW;
					kOld = k;
					ikF += 1;
					FieldType radius(1.2*Pi/float(nkSearch) * float(ikF));
					k[0] = FSCenters[iSheet][0] + radius * cos(Theta);
					k[1] = FSCenters[iSheet][1] + radius * sin(Theta);
					k[2] = kz;
					bands.getBands(k,w,v);
					FieldType ek(w[FSBand[iSheet]]);
					sgnW = (ek > 0) - (ek < 0);
				}
				if (ikF <= nkSearch) {
					FieldType kx(0.5*(k[0]+kOld[0]));
					FieldType ky(0.5*(k[1]+kOld[1]));
					// FieldType kz(0.5*(k[2]+kOld[2]));
					kFx.push_back(kx);
					kFy.push_back(ky);
					kFz.push_back(kz);
					kFtoBand.push_back(FSBand[iSheet]);
					size_t ic(kFx.size()-1);
					deltakF.push_back(calcDeltaKF(kFx[ic],kFy[ic],kFz[ic],iSheet,dim));
					vkF.push_back(calcVkF(kFx[ic],kFy[ic],kFz[ic],FSBand[iSheet],dim));
					gammaB1GkF.push_back(calcgammaB1GkF(kFx[ic],kFy[ic],kFz[ic],FSBand[iSheet],dim));
					if (calcOW_) calcWeights(ic);
				}

			}
		// }
	}

	void calcWeights(const size_t ic) {
			VectorType k(3,0);
			VectorType mk(3,0);
			VectorType w(param.nOrb);
			ComplexMatrixType v(param.nOrb,param.nOrb);
			ComplexMatrixType mv(param.nOrb,param.nOrb);
			k[0] = kFx[ic] ; k[1] = kFy[ic] ; k[2] = kFz[ic];
			mk[0] =-kFx[ic]; mk[1] =-kFy[ic]; mk[2] =-kFz[ic];
			bands.getBands(k,w,v);
			bands.getBands(mk,w,mv);
			for (size_t iorb=0;iorb<param.nOrb;iorb++) {
				FieldType re(real(v(iorb,kFtoBand[ic])));
				FieldType im(imag(v(iorb,kFtoBand[ic])));
				FieldType re2(real(mv(iorb,kFtoBand[ic])));
				FieldType im2(imag(mv(iorb,kFtoBand[ic])));
				weights[ic][iorb] =   ComplexType(re,im);
				weights2[ic][iorb] =  ComplexType(re2,im2);
			}
	}

	void calcKFOpen2D(const size_t nkSearch, const size_t iSheet, const size_t scanAlongDir) {
		VectorType w(param.nOrb);
		ComplexMatrixType v(param.nOrb,param.nOrb);

		for (size_t ik = 0; ik < nkPerSheet; ++ik) {
			size_t ikF(0);
			VectorType k(3,0);
			k[0] = FSCenters[iSheet][0];
			k[1] = FSCenters[iSheet][1];
			k[scanAlongDir==1?0:1] += float(ik)*2.*Pi/float(nkPerSheet);
			bands.getBands(k,w,v);
			FieldType ek(w[FSBand[iSheet]]);
			FieldType sgnW = (ek > 0) - (ek < 0);
			FieldType sgnOld(sgnW);
			VectorType kOld(k);
			while (ikF <= nkSearch && sgnW == sgnOld) {
				sgnOld = sgnW;
				kOld = k;
				ikF += 1;
				FieldType dist(Pi/float(nkSearch)*float(ikF));
				k[scanAlongDir] += dist;
				bands.getBands(k,w,v);
				FieldType ek(w[FSBand[iSheet]]);
				sgnW = (ek > 0) - (ek < 0);
			}
			if (ikF <= nkSearch) {
				FieldType kx(0.5*(k[0]+kOld[0]));
				FieldType ky(0.5*(k[1]+kOld[1]));
				kFx.push_back(kx);
				kFy.push_back(ky);
				kFz.push_back(0.0);
				kFtoBand.push_back(FSBand[iSheet]);
				size_t ic(kFx.size()-1);
				vkF.push_back(calcVkF(kFx[ic],kFy[ic],kFz[ic],FSBand[iSheet],2));
				if (calcOW_) calcWeights(ic);
			}
		}
		// Now calc. deltakF
		size_t ind0(kFx.size()-nkPerSheet);
		for (size_t ikF = 0; ikF < nkPerSheet-1; ++ikF)
		{
			FieldType dkF(sqrt(pow(kFx[ind0+ikF+1]-kFx[ind0+ikF],2)+pow(kFy[ind0+ikF+1]-kFy[ind0+ikF],2)));
			deltakF.push_back(dkF);
		}
		// Special treatment for last point
		FieldType dkF(0.0);
		if (scanAlongDir == 1) dkF = sqrt(pow(Pi-kFx[ind0+nkPerSheet-1],2)+pow(kFy[ind0]-kFy[ind0+nkPerSheet-1],2));
		if (scanAlongDir == 0) dkF = sqrt(pow(kFx[ind0]-kFx[ind0+nkPerSheet-1],2)+pow(Pi-kFy[ind0+nkPerSheet-1],2));
		deltakF.push_back(dkF);
	}

	FieldType calcDeltaKF(const FieldType& kFx, const FieldType& kFy, const FieldType& kFz,
						  const size_t iSheet,size_t dim) {
		FieldType rx(kFx-FSCenters[iSheet][0]);
		FieldType ry(kFy-FSCenters[iSheet][1]);
		// deltakF.push_back(sqrt(pow(rx,2)+pow(ry,2))*2.*Pi/float(nkPerSheet));
		if (dim==2) {
			return sqrt(pow(rx,2)+pow(ry,2))*2.*Pi/float(nkPerSheet);
		} else {
			FieldType vz (calcVkFz(kFx,kFy,kFz,FSBand[iSheet]));
			FieldType v2D(calcVkF(kFx,kFy,kFz,FSBand[iSheet],2));
			FieldType Theta(asin(vz/v2D));
			FieldType dkz(4.*param.pi_f/float(param.FSnkz));
			return sqrt(pow(rx,2)+pow(ry,2))
			       *2.*Pi/float(nkPerSheet)
			       *dkz/cos(Theta);
		}

	}

	void calcDeltaKF(const VectorType& kFx,const VectorType& kFy, VectorType& deltakF) {
		size_t nk(kFx.size());
		for (size_t ik = 0; ik < nk-1; ++ik)
		{
			deltakF[ik] = sqrt(pow(kFx[ik+1]-kFx[ik],2)+pow(kFy[ik+1]-kFy[ik],2));
		}
		deltakF[nk-1] = sqrt(pow(kFx[0]-kFx[nk-1],2)+pow(kFy[0]-kFy[nk-1],2));
	}

	void calcDeltaKF2D() {
		size_t nk(kFx.size());
		FieldType length, vecx,vecy;
		for (size_t ik = 0; ik < nk-1; ++ik)
		{
			vecx = kFx[ik+1]-kFx[ik]; 
			vecy = kFy[ik+1]-kFy[ik];
			length = sqrt(pow(vecx,2)+pow(vecy,2));
			if (length > 0.25*Pi) {
				if (conc.rank()==0) std::cout << "FS segment larger than pi/4. FS points not ordered? \n";
				exit(0);
			}
			deltakF[ik] = length;
		}
		vecx = kFx[0]-kFx[nk-1]; 
		vecy = kFy[0]-kFy[nk-1]; 
		// deltakF[nk-1] = sqrt(pow(vecx,2)+pow(vecy,2));
		deltakF[nk-1] = 0.0; // for open Fermi surface
	
	}
	// void calcDeltaKF2D() {
	// 	size_t nk(kFx.size());
	// 	FieldType length, vecx,vecy;
	// 	for (size_t ik = 1; ik < nk-1; ++ik)
	// 	{
	// 		vecx = 0.5*(kFx[ik+1]-kFx[ik-1]); 
	// 		vecy = 0.5*(kFy[ik+1]-kFy[ik-1]);
	// 		length = sqrt(pow(vecx,2)+pow(vecy,2));
	// 		deltakF[ik] = length;
	// 	}
	// 	vecx = kFx[1]-kFx[0]; 
	// 	deltakF[0] = sqrt(pow(vecx,2)+pow(vecy,2));
	// 	vecx = kFx[nk-1]-kFx[nk-2]; 
	// 	deltakF[nk-1] = sqrt(pow(vecx,2)+pow(vecy,2));
	// }

	FieldType calcVkF(const FieldType& kFx, const FieldType& kFy, const FieldType& kFz,
		              const size_t iband, size_t dim) {
		FieldType dk(0.0001);
		VectorType k(3,0);
		FieldType ekpx,ekpy,ekmx,ekmy,ekpz,ekmz;
		VectorType w(param.nOrb);
		ComplexMatrixType v(param.nOrb,param.nOrb);
		k[0] = kFx+dk; k[1] = kFy;
		bands.getBands(k,w,v) ; ekpx = w[iband];
		k[0] = kFx; k[1] = kFy + dk;
		bands.getBands(k,w,v) ; ekpy = w[iband];
		k[0] = kFx-dk; k[1] = kFy;
		bands.getBands(k,w,v) ; ekmx = w[iband];
		k[0] = kFx; k[1] = kFy - dk;
		bands.getBands(k,w,v) ; ekmy = w[iband];
		k[0] = kFx; k[1] = kFy; k[2] = kFz + dk;
		bands.getBands(k,w,v) ; ekpz = w[iband];
		k[0] = kFx; k[1] = kFy; k[2] = kFz - dk;
		bands.getBands(k,w,v) ; ekmz = w[iband];
		FieldType rx((ekpx-ekmx)/(2.*dk));
		FieldType ry((ekpy-ekmy)/(2.*dk));
		FieldType rz((ekpz-ekmz)/(2.*dk));
		if (dim==2) {
			return sqrt(pow(rx,2)+pow(ry,2));
		} else {
			return sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
		}
	}

	FieldType calcgammaB1GkF(const FieldType& kFx, const FieldType& kFy, const FieldType& kFz,
		           		   const size_t iband, size_t dim) {
		FieldType dk(0.0001);
		VectorType k(3,0);
		FieldType ekpx,ekpy,ekmx,ekmy,ek;
		VectorType w(param.nOrb);
		ComplexMatrixType v(param.nOrb,param.nOrb);
		k[0] = kFx; k[1] = kFy;
		bands.getBands(k,w,v) ;   ek = w[iband];
		k[0] = kFx+dk; k[1] = kFy;
		bands.getBands(k,w,v) ; ekpx = w[iband];
		k[0] = kFx; k[1] = kFy + dk;
		bands.getBands(k,w,v) ; ekpy = w[iband];
		k[0] = kFx-dk; k[1] = kFy;
		bands.getBands(k,w,v) ; ekmx = w[iband];
		k[0] = kFx; k[1] = kFy - dk;
		bands.getBands(k,w,v) ; ekmy = w[iband];

		FieldType fxx((ekpx+ekmx-2.*ek)/(pow(dk,2)));
		FieldType fyy((ekpy+ekmy-2.*ek)/(pow(dk,2)));

		return fxx - fyy;
	}

	FieldType calcVkFz(const FieldType& kFx, const FieldType& kFy, const FieldType& kFz,
		              const size_t iband) {
		FieldType dk(0.0001);
		VectorType k(3,0);
		FieldType ekmz,ekpz;
		VectorType w(param.nOrb);
		ComplexMatrixType v(param.nOrb,param.nOrb);
		k[0] = kFx; k[1] = kFy; k[2] = kFz + dk;
		bands.getBands(k,w,v) ; ekpz = w[iband];
		k[0] = kFx; k[1] = kFy; k[2] = kFz - dk;
		bands.getBands(k,w,v) ; ekmz = w[iband];
		FieldType rz((ekpz-ekmz)/(2.*dk));
		return rz;
	}

	void readFromFile() {
		if (conc.rank()==0) {
			std::string file(param.fsfile);
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format kFx,kFy,kFz,band,dk,vkF
			size_t step=6;
			// size_t step=5;
			std::cout << "File="<<file<<"\n";
			nTotal=data.size()/step;
			if (conc.rank()==0) std::cout << "nTotal=" << nTotal << "\n";
			for (size_t i = 0; i < nTotal; i++)
			{
				kFx.push_back     (data[i*step]);
				kFy.push_back     (data[i*step+1]);
				kFz.push_back     (data[i*step+2]);
				kFtoBand.push_back(size_t(data[i*step+3]));
				deltakF.push_back (data[i*step+4]);
				vkF.push_back     (data[i*step+5]);
			}
		}
		conc.broadcast(nTotal);
		if (conc.rank()!=0) {
			kFx.resize(nTotal);
			kFy.resize(nTotal);
			kFz.resize(nTotal);
			kFtoBand.resize(nTotal);
			deltakF.resize(nTotal);
			vkF.resize(nTotal);
		}
		conc.broadcast(kFx);
		conc.broadcast(kFy);
		conc.broadcast(kFz);
		conc.broadcast(kFtoBand);
		conc.broadcast(deltakF);
		conc.broadcast(vkF);
	}

	void readFromFile2() {
		if (conc.rank()==0) {
			std::string file(param.fsfile);
			VectorType data;
			loadVector(data,file);
			// We assume that each line has the format kFx,kFy,kFz,band, dk, vF
			size_t step=6;
			std::cout << "File="<<file<<"\n";
			nTotal=data.size()/step;
			std::cout << "nTotal=" << nTotal << "\n";
			for (size_t i = 0; i < nTotal; i++)
			{
				FieldType kx(data[i*step]);
				FieldType ky(data[i*step+1]);
				FieldType kz(data[i*step+2]);
				size_t band(data[i*step+3]);
				if (band+1 > param.nOrb) {
					if (conc.rank() == 0) std::cout << "Too many bands! Exiting ... \n";
					exit(0);
				}
				kFx.push_back(kx);
				kFy.push_back(ky);
				kFz.push_back(kz);
				kFtoBand.push_back(band);
				vkF.push_back(calcVkF(kx,ky,kz,band,param.dimension));
			}
			deltakF.resize(nTotal);
			gammaB1GkF.resize(nTotal);
			calcDeltaKF2D();
		}
		conc.broadcast(nTotal);
		if (conc.rank()!=0) {
			kFx.resize(nTotal);
			kFy.resize(nTotal);
			kFz.resize(nTotal);
			kFtoBand.resize(nTotal);
			deltakF.resize(nTotal);
			vkF.resize(nTotal);
		}
		conc.broadcast(kFx);
		conc.broadcast(kFy);
		conc.broadcast(kFz);
		conc.broadcast(kFtoBand);
		conc.broadcast(deltakF);
		conc.broadcast(vkF);

		std::cout << "Done reading in the FSfile\n";
	}

	void writeKF() {
		std::string filename = "FSpoints_" + param.fileID + ".txt";
		std::ofstream os(filename);
		int width(8);
		os.precision(width);
		os << std::fixed;
		for (size_t i=0;i<nTotal;i++) {
			os << kFx[i] << " , " << kFy[i] << " , " << kFz[i] << " , ";
			os << kFtoBand[i] << " , " << deltakF[i] << " , " << vkF [i] ; // << "\n";
			if (calcOW_) {
				// for (size_t iorb=0;iorb<param.nOrb;iorb++) os << " , " << abs(weights[i][iorb]) ;
				for (size_t iorb=0;iorb<param.nOrb;iorb++) os << " , " << real(weights[i][iorb])<< " , " << imag(weights[i][iorb]);
				// for (size_t iorb=0;iorb<param.nOrb;iorb++) os << " , " << real(weights2[i][iorb])<< " , " << imag(weights2[i][iorb]);
			}

			os << "\n";
		}
	}

};

}
#endif
