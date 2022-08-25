#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include <cstdlib> // for atof and atoi
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Matrix.h"
#include "Range.h"
#include "momentumDomain.h"
#include "parameters.h"
#include "utilities.h"
// #include "tbFromFile.h"
// #include "SrRuO.h"
// #include "SrRuO_SO.h"
// #include "1band_wSpin.h"
// #include "BaFeAs_5orb.h"
// #include "KFe2Se2.h"
// #include "FourOrbital.h"
// #include "bilayer.h"
// #include "coupledLadders.h"
#include "gaps3D.h"
#include "model.h"

namespace rpa {

template <typename Field, template <typename> class MatrixTemplate,
          typename ModelType, typename ConcurrencyType>
class bandstructure {

private:
  typedef MatrixTemplate<Field> MatrixType;
  typedef std::complex<Field> ComplexType;
  typedef MatrixTemplate<ComplexType> ComplexMatrixType;
  typedef std::vector<Field> VectorType;
  typedef std::vector<ComplexType> ComplexVectorType;
  typedef Field FieldType;
  typedef momentumDomain<Field, MatrixTemplate, ConcurrencyType> kDomain;
  typedef PsimagLite::Range<ConcurrencyType> RangeType;

  rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &param;
  ModelType &model;
  ConcurrencyType &conc;
  size_t nbands;
  int nLines;
  VectorType dx, dy, dz, ht;
  // VectorType hti;
  std::vector<size_t> orb1, orb2;
  kDomain kmesh_;
  const kDomain &kmesh;
  bool caching_;
  std::vector<bool> cachedK;
  ComplexMatrixType Lm;

public:
  std::vector<VectorType> ek;
  std::vector<std::vector<ComplexType>> gapk;
  std::vector<ComplexMatrixType> ak;
  std::vector<ComplexMatrixType> Mk;
  std::vector<ComplexMatrixType> MkFF;
  std::vector<VectorType> ekq;
  std::vector<std::vector<ComplexType>> gapkq;
  std::vector<ComplexMatrixType> akq;
  std::vector<ComplexMatrixType> Mkq;
  std::vector<ComplexMatrixType> MkqFF;

  // model<FieldType, MatrixTemplate, ConcurrencyType> model;

  // #ifdef USE_SRRUO
  // 			SrRuO_SO<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_1BANDWSPIN
  // 			SingleBand_wSpin<FieldType,MatrixTemplate,ConcurrencyType>
  // s; #elif USE_BILAYER
  // 			orthoIIBilayer<FieldType,MatrixTemplate,ConcurrencyType>
  // s;
  // 			// bilayer<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_BILAYER_1BAND
  // 			bilayer<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_BSCCObilayer
  // 			BSCCObilayer<FieldType,MatrixTemplate,ConcurrencyType>
  // s; #elif USE_BILAYER_FESC
  // 			bilayerFESC<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_BAFEAS
  // 			BaFeAs<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_KFE2SE2
  // 			KFe2Se2<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_FOURORBITAL
  // 			FourOrbital<FieldType,MatrixTemplate,ConcurrencyType> s;
  // #elif USE_COUPLEDLADDERS
  // 			coupledLadders<FieldType,MatrixTemplate,ConcurrencyType>
  // s; #else
  // tbFromFile<FieldType,MatrixTemplate,ConcurrencyType> s; #endif

public:
  bandstructure(
      rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &parameters,
      ModelType &modelIn, ConcurrencyType &concurrency, const kDomain &kmeshIn,
      bool caching)
      : param(parameters), model(modelIn), conc(concurrency),
        nbands(param.nOrb), kmesh_(param, conc, 0, 2), // not being used
        kmesh(kmeshIn),
        // caching_((kmesh.nktot<=2e6)?caching:false),
        caching_(caching), cachedK(kmesh.nktot, false),
        Lm(2 * nbands, 2 * nbands),
        ek(caching_ ? kmesh.nktot : 1, VectorType(nbands)),
        gapk(caching_ ? kmesh.nktot : 1, ComplexVectorType(nbands)),
        ak(caching_ ? kmesh.nktot : 1, ComplexMatrixType(nbands, nbands)),
        Mk(caching_ ? kmesh.nktot : 1,
           ComplexMatrixType(nbands, nbands * nbands)),
        MkFF(caching_ ? kmesh.nktot : 1,
             ComplexMatrixType(nbands, nbands * nbands)),
        ekq(caching_ ? kmesh.nktot : 1, VectorType(nbands)),
        gapkq(caching_ ? kmesh.nktot : 1, ComplexVectorType(nbands)),
        akq(caching_ ? kmesh.nktot : 1, ComplexMatrixType(nbands, nbands)),
        Mkq(caching_ ? kmesh.nktot : 1,
            ComplexMatrixType(nbands * nbands, nbands)),
        MkqFF(caching_ ? kmesh.nktot : 1,
              ComplexMatrixType(nbands * nbands, nbands)) {
    // if (kmesh.nktot>=16384) caching_=false;
    // if (param.tbfile!="") readCSVFile();
    if (conc.rank() == 0)
      std::cout << "Caching=" << caching_ << "\n";
    // if (param.LS==1) setupLMatrix();
    // if (param.sublattice==1) fixdr();
  }

  // Constructor without kmesh input
  bandstructure(
      rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &parameters,
      ModelType &modelIn, ConcurrencyType &concurrency)
      : param(parameters), model(modelIn), conc(concurrency),
        nbands(param.nOrb), kmesh_(param, conc, 0, 2), kmesh(kmesh_),
        caching_(true), cachedK(0), Lm(2 * nbands, 2 * nbands),
        ek(caching_ ? kmesh.nktot : 1, VectorType(nbands)),
        gapk(caching_ ? kmesh.nktot : 1, ComplexVectorType(nbands)),
        ak(caching_ ? kmesh.nktot : 1, ComplexMatrixType(nbands, nbands)),
        Mk(caching_ ? kmesh.nktot : 1,
           ComplexMatrixType(nbands, nbands * nbands)),
        MkFF(caching_ ? kmesh.nktot : 1,
             ComplexMatrixType(nbands, nbands * nbands)),
        ekq(caching_ ? kmesh.nktot : 1, VectorType(nbands)),
        gapkq(caching_ ? kmesh.nktot : 1, ComplexVectorType(nbands)),
        akq(caching_ ? kmesh.nktot : 1, ComplexMatrixType(nbands, nbands)),
        Mkq(caching_ ? kmesh.nktot : 1,
            ComplexMatrixType(nbands * nbands, nbands)),
        MkqFF(caching_ ? kmesh.nktot : 1,
              ComplexMatrixType(nbands * nbands, nbands)) {
    // if (kmesh.nktot>=16384) caching_=false;
    // if (param.tbfile!="") readCSVFile();
    if (conc.rank() == 0)
      std::cout << "Caching=" << caching_ << "\n";
    // if (param.LS==1) setupLMatrix();
    // if (param.sublattice==1) fixdr();
  }

  void getEkAndAk(VectorType &kvec, VectorType &eigenvals,
                  ComplexMatrixType &eigenvects, int spin = 1) {
    // if(!caching_) {
    getBands(kvec, eigenvals, eigenvects,
             spin); // Now included in getBands routine
    // phaseFactor(kvec,eigenvects); // Now included in getBands routine
    return;
    // } else {
    // 	// NOTE: This scheme only works when tb input file has all atoms within
    // a unit cell on the same site, i.e.
    // 	// when the phase-factors associated with positions with a unit cell are
    // taken into account here. If the tb input
    // 	// file has these shifts already in it, the phase factors will already
    // be taken care of and are included in the
    // 	// eigenvectors. But in this case, since the caching maps the input
    // k-vector to the 1. BZ, the phase factors
    // 	// won't be correct. So for now, we will switch off mapping to 1. BZ
    // 	// Also, when dim=2 and q has a finite qz, there seems to be a problem.
    // Therefore caching is turned off! 	size_t ik; FieldType residue;
    // VectorType kOrg(kvec); 	kmesh.kToik(kvec,ik,residue,false); // false =
    // mapping
    // to 1.BZ switched off 	if (fabs(residue)<=1.0e-5) { // kvec falls on
    // k[ik] in kmesh 		if (cachedK[ik]) {
    // 			// std::cout << "in cached routine \n";
    // 			eigenvals = ev[ik];
    // 			eigenvects = ak[ik];
    // 			// phaseFactor(kOrg,eigenvects); // Now included in
    // getBands routine 			return; 		} else {
    // VectorType kik(3); kmesh.momenta.getRow(ik,kik);
    // 			getBands(kik,eigenvals,eigenvects,spin);
    // 			ev[ik] = eigenvals;
    // 			ak[ik] = eigenvects;
    // 			cachedK[ik] = true;
    // 			// phaseFactor(kOrg,eigenvects); // Now included in
    // getBands routine 			return;
    // 		}
    // 	} else {
    // 		getBands(kOrg,eigenvals,eigenvects,spin);
    // 		// std::cout << "ik="<<ik<<" not yet cached\n";
    // 		// phaseFactor(kOrg,eigenvects);
    // 		return;
    // 	}
    // }
  }

  inline void getBands(const VectorType &k, VectorType &eigenvals,
                       ComplexMatrixType &eigenvects, int spin = 1) {

    model.getBands(k, eigenvals, eigenvects);
    return;
  }

  void precalculate_ekak(bool calcGap = 0) {
    size_t nktot(kmesh.nktot);
    std::vector<FieldType> k(3);
    for (size_t ik = 0; ik < nktot; ik++) {
      kmesh.momenta.getRow(ik, k);
      calculateBandTensors(k, ek[ik], ak[ik], Mk[ik], gapk[ik], MkFF[ik],
                           calcGap, 0);
    }
  }

  void precalculate_ekqakq(const VectorType &q, bool calcGap = 0) {
    size_t nktot(kmesh.nktot);
    std::vector<FieldType> k(3);
    std::vector<FieldType> kq(3);
    for (size_t ik = 0; ik < nktot; ik++) {
      kmesh.momenta.getRow(ik, k);
      for (size_t i = 0; i < 3; ++i)
        kq[i] = k[i] + q[i];
      calculateBandTensors(kq, ekq[ik], akq[ik], Mkq[ik], gapkq[ik], MkqFF[ik],
                           calcGap, 1);
    }
  }

  void calculateBandTensors(VectorType &k_, VectorType &ek_,
                            ComplexMatrixType &ak_, ComplexMatrixType &MGG,
                            std::vector<ComplexType> &gap,
                            ComplexMatrixType &MFF, bool calcGap = 0,
                            bool qshift = 0) {

    getEkAndAk(k_, ek_, ak_);
    VectorType mk_(k_);
    VectorType emk_(ek_);
    ComplexMatrixType amk_(ak_);
    for (size_t i = 0; i < param.dimension; i++)
      mk_[i] = -mk_[i];
    getEkAndAk(mk_, emk_, amk_);

    for (size_t i = 0; i < nbands * nbands; i++) {
      size_t l1 = param.indexToOrb(i, 0);
      size_t l2 = param.indexToOrb(i, 1);
      for (size_t b = 0; b < nbands; b++) {
        if (qshift)
          MGG(i, b) = ak_(l1, b) * conj(ak_(l2, b));
        else
          MGG(b, i) = ak_(l1, b) * conj(ak_(l2, b));
        if (calcGap) {
          if (param.explicitSpin &&
              param.oppositeSpinPairing ==
                  1) { // spin explicitely taken into account and up-down gap
            // if (!qshift) MFF(b,i) = conj(ak_(l1,b)) *
            // conj(ak_(l2,(b+int(nbands/2))%nbands)); // Pseudospin up-down
            // pairing else MFF(i,b) = ak_(l1,b) *
            // ak_(l2,(b+int(nbands/2))%nbands);
            if (!qshift)
              MFF(b, i) =
                  conj(ak_(l1, b)) *
                  conj(amk_(l2, (b + int(nbands / 2)) %
                                    nbands)); // Pseudospin up-down pairing
            else
              MFF(i, b) = amk_(l1, b) * ak_(l2, (b + int(nbands / 2)) % nbands);
          } else {
            // if (!qshift) MFF(b,i) = conj(ak_(l1,b)) * conj(ak_(l2,b));
            // else MFF(i,b) = ak_(l1,b) * ak_(l2,b);
            if (!qshift)
              MFF(b, i) = conj(ak_(l1, b)) * conj(amk_(l2, b));
            else
              MFF(i, b) = amk_(l1, b) * ak_(l2, b);
          }
        }
      }
    }

    if (calcGap) {
      // kmesh.mapTo1BZ(k_); // make sure k is in 1. BZ because gap may not be
      // periodic in k
      gap3D<FieldType, psimag::Matrix, ModelType, ConcurrencyType> Delta(
          param, model, conc);
      for (size_t i = 0; i < nbands; i++) {
        gap[i] = Delta(k_, i, ak_);
        // gap[i] *=
        // std::pow(param.Omega0,2)/(std::pow(ek_[i],2)+std::pow(param.Omega0,2));
        // // Lorentzian cut-off
        gap[i] *= exp(-std::pow(ek_[i], 2) /
                      std::pow(param.Omega0, 2)); // Gaussian cut-off
      }
    }
  }

  void writeGap() {
    if (!caching_)
      return;

    size_t nktot(kmesh.nktot);
    std::vector<FieldType> k(3);

    std::string cstr = "gap_" + param.fileID + ".txt";
    const char *filename = cstr.c_str();
    if (conc.rank() == 0) {
      // const char *filename = file.c_str();
      std::ofstream os(filename);
      int precision = 5;
      os.precision(precision);
      for (size_t ik = 0; ik < nktot; ik++) {
        std::vector<FieldType> k(3);
        kmesh.momenta.getRow(ik, k);
        os << k[0] / param.pi_f << "," << k[1] / param.pi_f << ","
           << k[2] / param.pi_f;
        for (size_t i = 0; i < nbands; ++i)
          os << "," << real(gapk[ik][i]);
        os << "\n";
      }
      os.close();
    }
  }

  // void precalculate_ekqakq(const VectorType& q, bool calcGap=0) {
  // 	size_t nktot(kmesh.nktot);
  // 	for (size_t ik = 0; ik < nktot; ik++) {
  // 		std::vector<FieldType> k(3);
  // 		kmesh.momenta.getRow(ik,k);
  // 		std::vector<FieldType> kq(3);
  // 		for (size_t i = 0; i < 3; ++i) kq[i] = k[i] + q[i];
  // 		getEkAndAk(kq,ekq[ik],akq[ik]);
  // 		for (size_t i=0; i < nbands*nbands; i++) {
  // 			size_t l1 = param.indexToOrb(i,0); size_t l2 =
  // param.indexToOrb(i,1); 			for (size_t b=0; b < nbands;
  // b++) { Mkq[ik](i,b) = akq[ik](l1,b) * conj(akq[ik](l2,b));
  // 				// if (ik==0) std::cout << i << "," << b << ","
  // << Mkq[ik](i,b) << "\n"; 				if (calcGap) {
  // if (param.explicitSpin) { 						MkqFF[ik](i,b) =
  // akq[ik](l1,b)
  // * akq[ik](l2,(b+int(nbands/2))%nbands); // Pseudospin up-down pairing }
  // else
  // { 						MkqFF[ik](i,b) = akq[ik](l1,b) *
  // akq[ik](l2,b);
  // 					}
  // 				}
  // 			}
  // 		}
  // 		if (calcGap) {
  // 			kmesh.mapTo1BZ(kq); // make sure kq is in 1. BZ because
  // gap may not be periodic in kz
  // gap3D<FieldType,psimag::Matrix,ConcurrencyType> Delta(param,conc);
  // for (size_t i=0; i < nbands; i++) { gapkq[ik][i] = Delta(kq,i,akq[ik]);
  // gapkq[ik][i] *=
  // std::pow(param.Omega0,2)/(std::pow(ekq[ik][i],2)+std::pow(param.Omega0,2));
  // // Lorentzian cut-off
  // 			}
  // 		}
  // 	}
  // }

  void calcBandStructure(std::string file) {
    size_t nktot(kmesh.nktot);
    RangeType range(0, nktot, conc);
    std::vector<std::vector<FieldType>> ek(nktot, VectorType(nbands, 0));
    std::vector<MatrixType> weights(nktot, MatrixType(nbands, nbands));
    ComplexMatrixType ak(nbands, nbands);

    for (; !range.end(); range.next()) {
      size_t ik = range.index();
      std::vector<FieldType> k(3);
      kmesh.momenta.getRow(ik, k);
      getEkAndAk(k, ek[ik], ak);

      for (size_t iband = 0; iband < nbands; iband++)
        for (size_t iorb = 0; iorb < nbands; iorb++)
          weights[ik](iorb, iband) = pow(abs(ak(iorb, iband)), 2);
    }

    for (size_t ik = 0; ik < nktot; ik++) {
      conc.reduce(ek[ik]);
      conc.reduce(weights[ik]);
    }

    if (conc.rank() == 0) {
      const char *filename = file.c_str();
      std::ofstream os(filename);
      int precision = 5;
      os.precision(precision);
      for (size_t ik = 0; ik < nktot; ik++) {
        std::vector<FieldType> k(3);
        kmesh.momenta.getRow(ik, k);
        os << k[0] / param.pi_f << "," << k[1] / param.pi_f << ","
           << k[2] / param.pi_f;
        for (size_t i = 0; i < nbands; ++i)
          os << "," << ek[ik][i];

        for (size_t iorb = 0; iorb < nbands; iorb++)
          for (size_t iband = 0; iband < nbands; iband++)
            os << "," << weights[ik](iorb, iband);
        os << "\n";
      }
      os.close();
    }
  }

  FieldType calcFilling() {
    size_t nktot(kmesh.nktot);
    RangeType range(0, nktot, conc);
    VectorType ek(nbands, 0);
    ComplexMatrixType ak(nbands, nbands);
    FieldType occ(0.0);

    VectorType occupation(nktot, 0);
    for (; !range.end(); range.next()) {
      size_t ik = range.index();
      std::vector<FieldType> k(3);
      kmesh.momenta.getRow(ik, k);

      getEkAndAk(k, ek, ak);
      for (size_t i = 0; i < nbands; i++)
        occupation[ik] += fermi(ek[i], 1. / param.temperature);
    }

    conc.allReduce(occupation);
    // for (size_t ik=0;ik<nktot;ik++) conc.allReduce(occupation[ik]);

    // std::cout << "Rank:" << conc.rank() <<
    // "Occupation:"<<occupation[0]<<"\n";

    for (size_t ik = 0; ik < nktot; ik++)
      occ += occupation[ik];
    occ /= FieldType(nktot);
    // if (conc.rank()==0) std::cout << "\n\t\tFilling = " << 2 * occ  <<
    // "\n\n";
    return 2 * occ;
  }

  inline FieldType fermi(const FieldType &energy, const FieldType &invT) {
    FieldType xx = energy * invT;
    FieldType f;
    if (xx <= 0.)
      f = 1. / (1. + exp(xx));
    else if (xx > 0.)
      f = exp(-xx) / (1. + exp(-xx));
    else
      f = 0.5;
    return f;
  }

  void eigen(VectorType &eigenvals, ComplexMatrixType &matrix) const {
    int n = matrix.n_row();
    int lwork = 2 * n - 1;
    std::vector<ComplexType> work(lwork);
    std::vector<FieldType> rwork(3 * n - 2);
    int info;

    GEEV('V', 'U', n, matrix, n, eigenvals, work, lwork, rwork, &info);
    if (info != 0) {
      throw std::runtime_error("zheev: failed\n");
    }
  }
};

} // namespace rpa

#endif
