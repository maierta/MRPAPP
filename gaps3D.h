//-*-C++-*-

#ifndef GAPS3D_H
#define GAPS3D_H

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib> // for atof and atoi
#include "Matrix.h"
#include "parameters.h"
#include "CrystalHarmonics2D.h"

namespace rpa {

    template<typename Field, template<typename> class MatrixTemplate,
			 // typename BandsType,
             typename ConcurrencyType>
    class gap3D {    // simple s+- for 5-orbital 1111 model 2D

    private:
        typedef Field                   FieldType;
        typedef MatrixTemplate<Field>   MatrixType;
        typedef MatrixTemplate<std::complex<Field> >   ComplexMatrixType;
        typedef std::complex<Field>     ComplexType;
        typedef std::vector<Field>      VectorType;
        typedef PsimagLite::Range<ConcurrencyType> RangeType;

        const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
        ConcurrencyType& conc;
        size_t nbands;
        std::vector<MatrixType> w;
        std::vector<std::vector<FieldType> > kz;
        std::vector<FieldType> k0;
		// momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh;
		// BandsType bands; // to get orbital weights
		// VectorType ek;
		// ComplexMatrixType ak;

        FieldType (*crystHarm)(const std::vector<FieldType>&, 
                               const MatrixType&,
                               const std::vector<FieldType>&,
                               const FieldType&,
                               const size_t);
        FieldType (*crystHarm2)(const std::vector<FieldType>&, 
                               const ComplexMatrixType&,
                               const size_t);

    public:
        gap3D() {}
        
        gap3D(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
              ConcurrencyType& concurrency):
            param(parameters),
            conc(concurrency),
            nbands(param.nOrb),
            w(0),
            kz(0),
            k0(0)  //,
			// kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension),
			// bands(param,conc,kmesh,false),
			// ek(param.nOrb,0),
			// ak(param.nOrb,param.nOrb)
        {           
            setParams3D();
        }

        ComplexType operator()(VectorType& k, size_t band,ComplexMatrixType& ak) {
            FieldType Delta(0.0);
			if (param.gAmpl=="KFeSe_overdoped_PhenD" ||
                param.gAmpl=="KFe2Se2_underdoped_PhenD") { // take the sign(Delta)
				// bands.getEkAndAk(k,ek,ak);	
	            // Delta = crystHarm2(k,ak,band) * param.Delta0;
				Delta = crystHarm(k,w[band],kz[band],k0[band],band); 
				FieldType sgnD = (Delta > 0) - (Delta < 0);
				Delta = sgnD*param.Delta0;
			} else if (param.gAmpl=="dwave_ladders") { // coupled 2-leg ladders in 2-orbital description
                if (real(ak(0,band))*real(ak(1,band))>0.0) { // bonding band
                    return param.Delta0 * (cos(k[0]) - cos(k[1]));
                } else { // antibonding band (shift kx by pi)
                    return param.Delta0 * (-cos(k[0]) - cos(k[1]));
                }

            } else {
				Delta = crystHarm(k,w[band],kz[band],k0[band],band) * param.Delta0;
			}
            return ComplexType(Delta,0.0);
        }
    
// This is the place where different gaps are hardcoded 
// in terms of amplitudes a,b,c for the s- or d-wave crystal harmonics

        void setParams3D() {
			if (param.gAmpl == "dwave" || param.gAmpl == "dwave_ladders") { // simple coskx-cosky for all bands
                crystHarm = &dwave;

                w.resize(param.nOrb,MatrixType(3,1));
				for (size_t i=0;i<param.nOrb;i++) w[i](0,0) =  1.0;
                kz.resize(param.nOrb,VectorType(0,0.0));
                k0.resize(param.nOrb,FieldType(0.0));
            } else if (param.gAmpl == "dwave_2D_KFe2Se2_RPA") { // coskx-cosky + 1.62*(cos2kx-cos2ky); see Maier et al., PRB 83, 100515(R)
                crystHarm = &dwave;

                w.resize(param.nOrb,MatrixType(3,1));
                w[3](0,0) =  1.0; w[3](1,0) = 1.62; // band 3 has the electron pocket
                kz.resize(param.nOrb,VectorType(0,0.0));
                k0.resize(param.nOrb,FieldType(0.0));
            } else if (param.gAmpl == "dwave_2D_KFe2Se2_iso") { // coskx-cosky (isotropic d-wave)
                crystHarm = &dwave;

                w.resize(param.nOrb,MatrixType(3,1));
                w[3](0,0) =  1.0; // band 3 has the electron pocket
                kz.resize(param.nOrb,VectorType(0,0.0));
                k0.resize(param.nOrb,FieldType(0.0));
            } else if (param.gAmpl == "spm-5band") { // simple s+ for all hole pockets and s- for all electron pockets
                crystHarm = &swave;

                w.resize(param.nOrb,MatrixType(3,1));
                w[1](0,0) =  1.0;
                w[2](0,0) =  1.0;
                w[3](0,0) =  -1.0;
                kz.resize(param.nOrb,VectorType(0,0.0));
                k0.resize(param.nOrb,FieldType(0.0));
			} else if (param.gAmpl == "LiFeAs_s_1") { // Yan's parametrization of 3D RPA gap for LiFeAs DFT model
                crystHarm = &swaveRPALiFeAs;

                w.resize(10,MatrixType(6,3));

                // beta_1^out (#1 in Yan's notation)
                w[6](0,0)=0.0164188716876;  w[6](0,1)=0.0205463228788;      w[6](0,2)=-0.000695636205009;
                w[6](1,0)=-0.0246717789076; w[6](1,1)=-0.0670007120592;     w[6](1,2)=0.00133341588286;
                w[6](2,0)=0.0291220104904;  w[6](2,1)=0.0565291453398;      w[6](2,2)=0.00427785462515;
                w[6](3,0)=-0.00708733403803;w[6](3,1)=0.00530496565299;     w[6](3,2)=-0.008279722707;
                w[6](4,0)=-0.00207173475993;w[6](4,1)=-0.00611523851534;    w[6](4,2)=0.00293178285281;
                w[6](5,0)=0.0147530260819;  w[6](5,1)=0.00955704309504;     w[6](5,2)=-0.00112472531463;
                
                // beta_1^in (#2 in Yan's notation)
                w[7](0,0)=0.00978761990877; w[7](0,1)=0.00363839310632;     w[7](0,2)=-0.00469193933527;
                w[7](1,0)=-0.00575725168811;w[7](1,1)=-0.0230025895224;     w[7](1,2)=0.0137676318191;
                w[7](2,0)=0.00370733772768; w[7](2,1)=-0.0176469302143;     w[7](2,2)=-0.00798590782462;
                w[7](3,0)=0.0108625157363;  w[7](3,1)=0.0630149051892;      w[7](3,2)=-0.0036242924301;
                w[7](4,0)=-0.00680517796442;w[7](4,1)=-0.0191650356472;     w[7](4,2)=0.00202674681228;
                w[7](5,0)=0.0150163953251;  w[7](5,1)=0.00764310716167;     w[7](5,2)=-0.000497278711838;
                
                // gamma^00 (#3 in Yan's notation)
                w[5](0,0)=0.00657796389308; w[5](0,1)=-0.00126703842849;    w[5](0,2)=0.000651263699903;
                w[5](1,0)=-0.00104849040636;w[5](1,1)=-0.000237020049388;   w[5](1,2)=-7.36207111786e-005;
                w[5](2,0)=0.00102196177054; w[5](2,1)=6.965951275e-005;     w[5](2,2)=-8.23913495473e-005;
                w[5](3,0)=0.000343614546984;w[5](3,1)=-0.00171290028559;    w[5](3,2)=0.000530338327369;
                w[5](4,0)=-0.00187202342238;w[5](4,1)=-0.00181313122651;    w[5](4,2)=0.000475275095703;
                w[5](5,0)=0.0080066663732;  w[5](5,1)=0.00290330559656;     w[5](5,2)=-0.000529553457368;
                
                // alpha_2^00 (#4 in Yan's notation)
                w[4](0,0)=-0.0888760887281; w[4](0,1)=0.189019944018;       w[4](0,2)=-0.0294760823422;
                w[4](1,0)=0.134154727592;   w[4](1,1)=-0.269317076461;      w[4](1,2)=0.0413057000783;
                w[4](2,0)=-0.0154969018701; w[4](2,1)=0.0347961329169;      w[4](2,2)=-0.00480436220377;
                w[4](3,0)=0.00676182079918; w[4](3,1)=-0.0153777152771;     w[4](3,2)=0.00209370957207;
                w[4](4,0)=0.00420109877936; w[4](4,1)=-0.0384892536849;     w[4](4,2)=0.00258000350961;
                w[4](5,0)=0.0264538289761;  w[4](5,1)=-0.00762314476293;    w[4](5,2)=0.00556621962384;
                
                // alpha_1^00 (#5 in Yan's notation)
                w[3](0,0)=-0.373835551892;  w[3](0,1)=0.387872011932;       w[3](0,2)=0.0226766191301;
                w[3](1,0)=0.463328086015;   w[3](1,1)=-0.465283614164;      w[3](1,2)=-0.0458555723719;
                w[3](2,0)=-0.0488232659442; w[3](2,1)=0.043464245297;       w[3](2,2)=0.00932837400899;
                w[3](3,0)=-0.0171986845078; w[3](3,1)=0.0178982711254;      w[3](3,2)=-0.00271541902697;
                w[3](4,0)=0.206661925228;   w[3](4,1)=-0.310493530687;      w[3](4,2)=0.0693309000131;
                w[3](5,0)=-0.172968867498;  w[3](5,1)=0.3044511094;         w[3](5,2)=-0.0952729099662;

                kz.resize(10,VectorType(6,0));
                for (size_t i=0;i<6;i++) {
                    FieldType r1(2.0*param.pi_f/10.*float(i));
                    kz[6][i] = r1;
                    kz[7][i] = r1;
                    kz[5][i] = r1;
                    kz[4][i] = r1;
                }
                FieldType kzmax(1.80448863);
                for (size_t i=0;i<6;i++) kz[3][i] = 2.0*kzmax/10.*float(i); 
                k0.resize(10,FieldType(0.0));
                FieldType r1(0.25*2.*param.pi_f/10.);
                FieldType r2(0.25*2.*kzmax/10.);
                k0[6] = r1; k0[7] = r1; k0[5] = r1; k0[4] = r1;
                k0[3] = r2;
            } else if (param.gAmpl == "LiFeAs_ARPES") { // Yan's parametrization of 3D RPA gap for LiFeAs ARPES model
                crystHarm = &swaveRPALiFeAs_2;

                w.resize(10,MatrixType(7,4));

                // beta_1^out (#1 in Yan's notation)
                w[6](0,0)=-6.56695815482e-005; w[6](0,1)=0.00171223921426  ; w[6](0,2)=-0.0114817219308  ; w[6](0,3)=-0.00248099230716;
                w[6](1,0)=-0.0224037715483   ; w[6](1,1)=0.0158773832694   ; w[6](1,2)=-0.0185523064403  ; w[6](1,3)=-0.00406282205477;
                w[6](2,0)=0.0100196626674    ; w[6](2,1)=-0.0110287097468  ; w[6](2,2)=0.009506903271    ; w[6](2,3)=0.00338179815868;
                w[6](3,0)=0.017222141969     ; w[6](3,1)=-0.0159647116753  ; w[6](3,2)=0.0160534551557   ; w[6](3,3)=0.00527565613274;
                w[6](4,0)=0.0129128476276    ; w[6](4,1)=-0.0130846741481  ; w[6](4,2)=0.0112541867994   ; w[6](4,3)=0.00518843236881;
                w[6](5,0)=0.00531922489282   ; w[6](5,1)=0.000812009842013 ; w[6](5,2)= 0.00607131322874 ; w[6](5,3)=-0.00143746021959;
                w[6](6,0)=0.0422787672406    ; w[6](6,1)=-0.0353666083809  ; w[6](6,2)=0.0273792086753   ; w[6](6,3)=0.00911667869103;
                
                // beta_1^in (#2 in Yan's notation)
                w[7](0,0)=0.0152000605901       ; w[7](0,1)=-0.0180498678213    ; w[7](0,2)=0.0050227222566     ; w[7](0,3)=0;
                w[7](1,0)=-0.0112638257915      ; w[7](1,1)=-0.00643401925198   ; w[7](1,2)=-0.00529518896178   ; w[7](1,3)=0;
                w[7](2,0)=-0.00979028187138     ; w[7](2,1)=0.00565022920957    ; w[7](2,2)=-0.00954145593514   ; w[7](2,3)=0;
                w[7](3,0)=0.000296808006564     ; w[7](3,1)=0.0326562720763     ; w[7](3,2)=-0.00800943456726   ; w[7](3,3)=0;
                w[7](4,0)=0.00630768008102      ; w[7](4,1)=-0.0223791703763    ; w[7](4,2)=0.0102415536241     ; w[7](4,3)=0;
                w[7](5,0)=0.0081781175208       ; w[7](5,1)=-0.00383800446676   ; w[7](5,2)=0.00891569046887    ; w[7](5,3)=0;
                w[7](6,0)=0.02487722431         ; w[7](6,1)=-0.0190985159494    ; w[7](6,2)=0.012654546818      ; w[7](6,3)=0;
                
                // gamma^00 (#3 in Yan's notation)
                w[5](0,0)=0.012867761808        ; w[5](0,1)=0.00519028557836    ; w[5](0,2)=-0.00164246689608   ; w[5](0,3)=0.000470167473084;
                w[5](1,0)=0.00096934048074      ; w[5](1,1)=-0.0021357928789    ; w[5](1,2)=-0.00244604760663   ; w[5](1,3)=0.000566119353277;
                w[5](2,0)=-0.00163485728495     ; w[5](2,1)=0.000803659912961   ; w[5](2,2)=0.00188532091318    ; w[5](2,3)=-0.000243787756411;
                w[5](3,0)=0.00231131377891      ; w[5](3,1)=-0.00115769628909   ; w[5](3,2)=-0.00364970719156   ; w[5](3,3)=0.000132424421841;
                w[5](4,0)=-0.00218565194971     ; w[5](4,1)=0.00100412607556    ; w[5](4,2)=0.00301487284463    ; w[5](4,3)=-0.000158767048429;
                w[5](5,0)=0.00234146150019      ; w[5](5,1)=-0.00453968125309   ; w[5](5,2)=-0.0085443641168    ; w[5](5,3)=0.000618057487816;
                w[5](6,0)=0.0105937655445       ; w[5](6,1)=0.00877043729364    ; w[5](6,2)=0.00749362478675    ; w[5](6,3)=0.000474607100088;
                
                // alpha_2^00 (#4 in Yan's notation)
                w[4](0,0)=0.0287596792928       ; w[4](0,1)=0.00148940264102    ; w[4](0,2)=-0.0131966265751    ; w[4](0,3)=0;
                w[4](1,0)=-0.0441773750004      ; w[4](1,1)=-0.0119456370984    ; w[4](1,2)=0.0279630984433     ; w[4](1,3)=0;
                w[4](2,0)=0.0132858890607       ; w[4](2,1)=-0.000132920774999  ; w[4](2,2)= -0.00648869199961  ; w[4](2,3)=0;
                w[4](3,0)=0.103768076003        ; w[4](3,1)=0.0604251724981     ; w[4](3,2)=-0.0820231587683    ; w[4](3,3)=0;
                w[4](4,0)=-0.658764695497       ; w[4](4,1)=-0.353724213692     ; w[4](4,2)=0.506119771688      ; w[4](4,3)=0;
                w[4](5,0)=1.10649429827         ; w[4](5,1)=0.649370630588      ; w[4](5,2)=-0.878185101079     ; w[4](5,3)=0;
                w[4](6,0)=-0.51976415934        ; w[4](6,1)=-0.347657947609     ; w[4](6,2)=0.435378594994      ; w[4](6,3)=0;
                
                // alpha_1^00 (#5 in Yan's notation)
                w[3](0,0)=0.00523174507501      ; w[3](0,1)=0                   ; w[3](0,2)=0                   ; w[3](0,3)=0;
                w[3](1,0)=-0.00194848010731     ; w[3](1,1)=0                   ; w[3](1,2)=0                   ; w[3](1,3)=0;
                w[3](2,0)=0.000400554956348     ; w[3](2,1)=0                   ; w[3](2,2)=0                   ; w[3](2,3)=0;
                w[3](3,0)=-0.000250287997696    ; w[3](3,1)=0                   ; w[3](3,2)=0                   ; w[3](3,3)=0;
                w[3](4,0)=6.74819527168e-005    ; w[3](4,1)=0                   ; w[3](4,2)=0                   ; w[3](4,3)=0;
                w[3](5,0)=-0.00042654993804     ; w[3](5,1)=0                   ; w[3](5,2)=0                   ; w[3](5,3)=0;
                w[3](6,0)=0.00582454646185      ; w[3](6,1)=0                   ; w[3](6,2)=0                   ; w[3](6,3)=0;
                
                kz.resize(10,VectorType(7,0));
                FieldType r1(0.0);
                for (size_t i=0;i<7;i++) {
                    r1 = param.pi_f/6.*float(i);
                    kz[6][i] = r1;
                    kz[7][i] = r1;
                    kz[5][i] = r1;
                }
                FieldType kzmax(1.626826);
                for (size_t i=0;i<7;i++) {
                    r1 = param.pi_f - (float(i)*(param.pi_f-kzmax))/6.; 
                    kz[4][i] = r1; 
                    kz[3][i] = r1;
                }
                k0.resize(10,FieldType(0.0));
                r1 = param.pi_f/24.;
                FieldType r2((param.pi_f-kzmax)/24.);
                k0[6] = r1; k0[7] = r1; k0[5] = r1; 
                k0[4] = r2; k0[3] = r2;

            } else if (param.gAmpl == "KFe2Se2_underdoped_d" ||
                       param.gAmpl == "KFe2Se2_underdoped_PhenD") { // Yan's parametrization of 3D d-wave RPA gap for KFe2Se2 underdoped n=6.12 case
                crystHarm = &dwaveRPAKFe2Se2_1;

                w.resize(10,MatrixType(9,4));

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0)=-0.000146835029238 ; w[6](0,1)=-0.0464858830987 ; w[6](0,2)=0.000623916663848  ; w[6](0,3)=0.0233728507709 ; 
                w[6](1,0)=0.00124153919088   ; w[6](1,1)=0.0374087132367  ; w[6](1,2)=0.00679307869731   ; w[6](1,3)=-0.01681046124  ;
                w[6](2,0)=0.000560116252357  ; w[6](2,1)=0.0613145235329  ; w[6](2,2)=-0.00315404851892  ; w[6](2,3)=-0.0262617656602;
                w[6](3,0)=0.00264572508287   ; w[6](3,1)=0.0498196080521  ; w[6](3,2)=-0.0244187510132   ; w[6](3,3)=-0.0200433454023;
                w[6](4,0)=-0.000177994442999 ; w[6](4,1)=-0.266242037149  ; w[6](4,2)=0.00108480820769   ; w[6](4,3)=0.111780249425  ;
                w[6](5,0)=-0.00252649855587  ; w[6](5,1)=0.0553204990297  ; w[6](5,2)=0.0234111501851    ; w[6](5,3)=-0.022261248878 ;
                w[6](6,0)=-0.000323964143232 ; w[6](6,1)=0.0526000431378  ; w[6](6,2)=0.00258058503177   ; w[6](6,3)=-0.0222280805654;
                w[6](7,0)=-0.00159254008946  ; w[6](7,1)=0.0436710570808  ; w[6](7,2)=-0.00737678739033  ; w[6](7,3)=-0.0196829744804;
                w[6](8,0)=0.000352177594024  ; w[6](8,1)=-0.0445591815957 ; w[6](8,2)=-9.69269741202e-005; w[6](8,3)=0.0222585151083 ;
                
                // delta_1^in (#2 in Yan's notation)
                w[7](0,0)=-0.00171186885659  ; w[7](0,1)=0.00886028356862 ; w[7](0,2)=0.00476041297724   ; w[7](0,3)=-0.004750891894  ;
                w[7](1,0)=-0.00339590923358  ; w[7](1,1)=-0.139253267293  ; w[7](1,2)=0.0103647077018    ; w[7](1,3)=0.0510998288678  ;
                w[7](2,0)=-5.14117117339e-005; w[7](2,1)=0.114327595797   ; w[7](2,2)=-0.00419305933221  ; w[7](2,3)=-0.0418855530609 ;
                w[7](3,0)=6.20741656683e-005 ; w[7](3,1)=0.0402686206156  ; w[7](3,2)=-0.00163863513336  ; w[7](3,3)=-0.0152032965508 ;
                w[7](4,0)=0.000553136529395  ; w[7](4,1)=0.0171686542699  ; w[7](4,2)=-0.00136016844513  ; w[7](4,3)=-0.00580222839948;
                w[7](5,0)=-0.000740425816951 ; w[7](5,1)=0.0424038734303  ; w[7](5,2)=0.00336932153652   ; w[7](5,3)=-0.0160047208895 ;
                w[7](6,0)=-0.000194254673806 ; w[7](6,1)=0.115548222159   ; w[7](6,2)=0.0050645517791    ; w[7](6,3)=-0.042600574478  ;
                w[7](7,0)=0.00358010361814   ; w[7](7,1)=-0.146645416407  ; w[7](7,2)=-0.0110783950203   ; w[7](7,3)=0.0538594120761  ;
                w[7](8,0)=0.00175134755737   ; w[7](8,1)=0.0140586738551  ; w[7](8,2)=-0.00477555701052  ; w[7](8,3)=-0.00639546151349;

                
                kz.resize(10,VectorType(9,0.0));
                FieldType r1(0.0);
                for (size_t i=0;i<9;i++) {
                    r1 = 4.*param.pi_f/16.*float(i);
                    kz[6][i] = r1;
                    kz[7][i] = r1;
                }
                k0.resize(10,FieldType(0.0));
                r1 = param.pi_f/16.;
                k0[6] = r1; k0[7] = r1; 
            } else if (param.gAmpl == "KFe2Se2_underdoped_Mazin") { // Pure 3D Mazin state for KFe2Se2 underdoped n=6.12 case
                crystHarm = &MazinKFe2Se2;

                w.resize(10,MatrixType(1,1));

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0) =  1.0;
                // delta_1^in (#2 in Yan's notation)
                w[7](0,0) = -1.0;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            } else if (param.gAmpl == "KFe2Se2_gap_all") { // in phase gaps on all bands
                crystHarm = &MazinKFe2Se2;

                w.resize(10,MatrixType(1,1));

                // delta_1^out (#1 in Yan's notation)
                for (size_t ib=0;ib<nbands;ib++) w[ib](0,0) =  1.0;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            } else if (param.gAmpl == "KFe2Se2_edoped_Mazin") {
                crystHarm = &MazinKFe2Se2_wZ;

                w.resize(10,MatrixType(1,1));

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0) =  1.0;
                // delta_1^in (#2 in Yan's notation)
                w[7](0,0) = -1.0;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            } else if (param.gAmpl == "KFe2Se2_edoped_s++") {
                crystHarm = &MazinKFe2Se2_wZ;

                w.resize(10,MatrixType(1,1));

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0) =  1.0;
                // delta_1^in (#2 in Yan's notation)
                w[7](0,0) =  1.0;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            } else if (param.gAmpl == "KFe2Se2_edoped_d" || param.gAmpl == "KFeSe_overdoped_PhenD") { // Yan's parametrization of 3D d-wave RPA gap for KFe2Se2 electron doped n=6.2 case
                crystHarm = &dwaveRPAKFe2Se2_2;

                w.resize(10,MatrixType(6,4));

                // Note: Here we only setup the parameters for the open cylindrical electron pockets
                // The closed Z-pocket will be treated separately in dwaveRPAKFe2Se2_2 in CrystalHarmonics2D.h

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0)=-0.000376574254776 ; w[6](0,1)=0.00534803786464 ; w[6](0,2)=-0.00214944046528 ; w[6](0,3)=-0.00187465318503;
                w[6](1,0)=0.00378238017884   ; w[6](1,1)=-0.0439946250568 ; w[6](1,2)=0.0185213992803   ; w[6](1,3)=0.0288458107427;
                w[6](2,0)=0.00278525714697   ; w[6](2,1)=0.0384716073571  ; w[6](2,2)=-0.0420161072751  ; w[6](2,3)=-0.0251851927468;
                w[6](3,0)=-0.00245686417688  ; w[6](3,1)=0.0419046162243  ; w[6](3,2)=0.0430655044116   ; w[6](3,3)=-0.0284484516654;
                w[6](4,0)=-0.00419703475022  ; w[6](4,1)=-0.0461417984305 ; w[6](4,2)=-0.0197214817822  ; w[6](4,3)=0.0306719512114;
                w[6](5,0)=0.000467861780322  ; w[6](5,1)=0.00594295598114 ; w[6](5,2)=0.00232567634318  ; w[6](5,3)=-0.00231237440549;
                
                // delta_1^in (#2 in Yan's notation)
                w[7](0,0)=2.10636831722e-005 ; w[7](0,1)=0.00461572679488  ; w[7](0,2)=0.000814669581186 ; w[7](0,3)=-4.73567242256e-005;
                w[7](1,0)=-0.00335650333574  ; w[7](1,1)=0.00866267396061  ; w[7](1,2)=-0.0223425955457  ; w[7](1,3)=-0.0126869798646;
                w[7](2,0)=-0.00286959903879  ; w[7](2,1)=-0.0151165264446  ; w[7](2,2)=0.0565726121546   ; w[7](2,3)=0.0147252712602;
                w[7](3,0)=0.00263497545636   ; w[7](3,1)=-0.0193973893796  ; w[7](3,2)=-0.0590280563058  ; w[7](3,3)=0.0179297738751;
                w[7](4,0)=0.0037302837559    ; w[7](4,1)=0.0259317712993   ; w[7](4,2)=0.0264777293781   ; w[7](4,3)=-0.0215736713774;
                w[7](5,0)=-0.000119758571572 ; w[7](5,1)=-0.00425927235986 ; w[7](5,2)=-0.00267437899451 ; w[7](5,3)=0.00433047399819;
                
                // kappa^pipi (#3 in Yan's notation)
                // Note: The kappa pockets are generated from band 6 also so that we will take care of it in 
                // dwaveRPAKFe2Se2_2 in CrystalHarmonics2D.h

                kz.resize(10,VectorType(6,0.0));
                FieldType r1(0.0);
                FieldType nz(10.);
                for (size_t i=0;i<6;i++) {
                    r1 = 4.*param.pi_f/nz*float(i);
                    kz[6][i] = r1;
                    kz[7][i] = r1;
                }
                k0.resize(10,FieldType(0.0));
                r1 = param.pi_f/nz;
                k0[6] = r1; k0[7] = r1; 

            } else if (param.gAmpl == "KFe2Se2_edoped_s") { // Yan's parametrization of 3D s-wave RPA gap for KFe2Se2 electron doped n=6.2 case
                crystHarm = &swaveRPAKFe2Se2_elDoped;

                w.resize(10,MatrixType(6,4));

                // Note: Here we only setup the parameters for the open cylindrical electron pockets
                // The closed Z-pocket will be treated separately in dwaveRPAKFe2Se2_2 in CrystalHarmonics2D.h

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0)=-0.00541906931057 ; w[6](0,1)=-0.0608620450495 ; w[6](0,2)=-0.00218074698119  ; w[6](0,3)=0.000213793212534;
                w[6](1,0)=-0.159900745092   ; w[6](1,1)=-0.322024169464  ; w[6](1,2)=-0.362989325643    ; w[6](1,3)=-0.151198716254;
                w[6](2,0)=0.146948189182    ; w[6](2,1)=0.318316684241   ; w[6](2,2)=0.32275570624      ; w[6](2,3)=0.133743437695;
                w[6](3,0)=0.166711712337    ; w[6](3,1)=0.346215240894   ; w[6](3,2)=0.37324596359      ; w[6](3,3)=0.155210651093;
                w[6](4,0)=-0.176066400173   ; w[6](4,1)=-0.362064801518  ; w[6](4,2)=-0.402014776722    ; w[6](4,3)=-0.168112517042;
                w[6](5,0)=0.00376517756956  ; w[6](5,1)=-0.0321390865253 ; w[6](5,2)=0.0186857373038    ; w[6](5,3)=0.00931860794524;
                
                // delta_1^in (#2 in Yan's notation)
                w[7](0,0)=-0.0382853531086  ; w[7](0,1)=-0.137877987437  ; w[7](0,2)=-0.0423493896102   ; w[7](0,3)=-0.0142972247143;
                w[7](1,0)=0.092468430053    ; w[7](1,1)=0.298610737216   ; w[7](1,2)=0.0908793850882    ; w[7](1,3)=0.0337750152989;
                w[7](2,0)=-0.146241137909   ; w[7](2,1)=-0.421956480362  ; w[7](2,2)=-0.16954081477     ; w[7](2,3)=-0.0606662966785;
                w[7](3,0)=-0.168920311938   ; w[7](3,1)=-0.489603454608  ; w[7](3,2)=-0.202883654055    ; w[7](3,3)=-0.0724826887294;
                w[7](4,0)=0.0931291894257   ; w[7](4,1)=0.29137064612    ; w[7](4,2)=0.11297526081      ; w[7](4,3)=0.0434331608326;
                w[7](5,0)=-0.0305227829668  ; w[7](5,1)=-0.107572947419  ; w[7](5,2)=-0.0483881055599   ; w[7](5,3)=-0.0179197623905;
                
                // kappa^pipi (#3 in Yan's notation)
                // Note: The kappa pockets are generated from band 6 also so that we will take care of it in 
                // dwaveRPAKFe2Se2_2 in CrystalHarmonics2D.h

                kz.resize(10,VectorType(6,0.0));
                FieldType r1(0.0);
                FieldType nz(10.);
                for (size_t i=0;i<6;i++) {
                    r1 = 4.*param.pi_f/nz*float(i);
                    kz[6][i] = r1;
                    kz[7][i] = r1;
                }
                k0.resize(10,FieldType(0.0));
                r1 = param.pi_f/nz;
                k0[6] = r1; k0[7] = r1; 

            } else if (param.gAmpl == "KFe2Se2_edoped_SO0p05_d") { // Yan's parametrization of 3D d-wave RPA gap for KFe2Se2 electron doped n=6.2 case
                crystHarm = &dwaveRPAKFe2Se2_3;

                w.resize(10,MatrixType(6,4));

                // Note: Here we only setup the parameters for the open cylindrical electron pockets
                // The closed Z-pocket will be treated separately in dwaveRPAKFe2Se2_3 in CrystalHarmonics2D.h

                // delta_1^out (#1 in Yan's notation)
                w[6](0,0)= 0.00034125016618 ; w[6](0,1)=-0.000863895377164 ; w[6](0,2)=0.000708626805186 ; w[6](0,3)=0.00238960938588;
                w[6](1,0)= 0.00256725824231 ; w[6](1,1)=-0.00616225343566  ; w[6](1,2)=0.00177281138308  ; w[6](1,3)=0.00414999069821;
                w[6](2,0)= 0.00499553885784 ; w[6](2,1)= 0.00524222392915  ; w[6](2,2)=-0.0042122058828  ; w[6](2,3)=-0.00344259875465;
                w[6](3,0)=-0.00502393716805 ; w[6](3,1)= 0.0101684871774   ; w[6](3,2)=0.00371464189582  ; w[6](3,3)=-0.00788840161221;
                w[6](4,0)=-0.00253118403358 ; w[6](4,1)=-0.00676612530322  ; w[6](4,2)=-0.00148114885328 ; w[6](4,3)=0.00485880813237;
                w[6](5,0)=-0.00027646649292 ; w[6](5,1)=-0.000572812391484 ; w[6](5,2)=-0.000291979616997; w[6](5,3)=0.00233739565051;

                // delta_1^in (#2 in Yan's notation)
                w[7](0,0)=-0.000408835674636; w[7](0,1)=0.00875030000196 ; w[7](0,2)= -0.00177112891004 ; w[7](0,3)= -0.00404843722616;
                w[7](1,0)=-0.00305911423061 ; w[7](1,1)=-0.0152513427525 ; w[7](1,2)= -0.000116686887307; w[7](1,3)= 0.00402298980713;
                w[7](2,0)=-0.00307709554551 ; w[7](2,1)=0.0132677188101  ; w[7](2,2)= 0.0178416530675   ; w[7](2,3)= -0.00272195510584;
                w[7](3,0)= 0.00303427066314 ; w[7](3,1)=0.0167483821715  ; w[7](3,2)= -0.0170639035353  ; w[7](3,3)= -0.0035634917666;
                w[7](4,0)= 0.00273019016491 ; w[7](4,1)=-0.0166269114833 ; w[7](4,2)= -0.00105721250432 ; w[7](4,3)= 0.0050031551408;
                w[7](5,0)= 0.000614882564338; w[7](5,1)=0.00821397798803 ; w[7](5,2)= 0.00222879169371  ; w[7](5,3)= -0.00429517883495;

                // kappa^pipi (#3 in Yan's notation)
                // Note: The kappa pockets are generated from band 6 also so that we will take care of it in 
                // dwaveRPAKFe2Se2_3 in CrystalHarmonics2D.h

                kz.resize(10,VectorType(6,0.0));
                FieldType r1(0.0);
                FieldType nz(10.);
                for (size_t i=0;i<6;i++) {
                    r1 = 4.*param.pi_f/nz*float(i);
                    kz[6][i] = r1;
                    kz[7][i] = r1;
                }
                k0.resize(10,FieldType(0.0));
                r1 = param.pi_f/nz;
                k0[6] = r1; k0[7] = r1; 

            // } else if (param.gAmpl == "KFe2Se2_underdoped_PhenD") {
                // crystHarm = &dwavePhenomKFe2Se2;
            // } else if (param.gAmpl == "KFeSe_overdoped_PhenD") {
                // crystHarm2 = &dwavePhenomKFe2Se2_overdoped;
                
            } else if (param.gAmpl == "10-orbit_7Sheets_s+-") { // standard s+- gap for 10-orbital model with 7 sheets
                crystHarm = &MazinKFe2Se2_wZ;

                w.resize(10,MatrixType(1,1));

                // 1. hole pocket
                w[3](0,0) =  1.0;
                // 2. hole pocket
                w[4](0,0) =  1.0;
                // 3. hole pocket
                w[5](0,0) =  1.0;
                // outer electron pocket
                w[6](0,0) = -1.0;
                // inner electron pocket
                w[7](0,0) = -1.0;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            } else if (param.gAmpl == "10-orbit_7Sheets_s+-_aniso") { // standard s+- gap for 10-orbital model with 7 sheets and anisotropic gap on electron sheets
                crystHarm = &MazinKFe2Se2_wZ;

                w.resize(10,MatrixType(1,1));

                // 1. hole pocket
                w[3](0,0) =  1.0;
                // 2. hole pocket
                w[4](0,0) =  1.0;
                // 3. hole pocket
                w[5](0,0) =  1.0;
                // outer electron pocket
                w[6](0,0) = -1.0;
                // inner electron pocket
                w[7](0,0) = -0.5;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            } else if (param.gAmpl == "10-orbit_7Sheets_AntiPhase") { // anti-phase s+- gap for 10-orbital model with 5 sheets
                crystHarm = &MazinKFe2Se2_wZ;

                w.resize(10,MatrixType(1,1));

                // 1. hole pocket
                w[3](0,0) =  1.0;
                // 2. hole pocket
                w[4](0,0) =  1.0;
                // 3. hole pocket
                w[5](0,0) =  1.0;
                // outer electron pocket
                w[6](0,0) = -1.0;
                // inner electron pocket
                w[7](0,0) = -1.0;

                kz.resize(10,VectorType(0,0.0));
                k0.resize(10,FieldType(0.0));

            }  else {
                std::cout << "Gap not implemented! Bailing out.\n";
                exit(0);
            }
        }

    };
}


#endif
