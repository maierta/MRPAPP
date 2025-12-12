//-*-C++-*-

#ifndef SUSCEPTIBILITY_H
#define SUSCEPTIBILITY_H

#include "Fermi.h"
#include "Matrix.h"
#include "bandstructure.h"
#include "chi0.h"
#include "ferminator.h"
#include "gaps2D.h"
#include "gaps3D.h"
#include "model.h"
#include "momentumDomain.h"
#include "parameters.h"
#include "rpa_CuO.h"
#include "utilities.h"
#include <fstream>
#include <string>
#include <vector>

namespace rpa {

template <typename Field, typename SuscType, typename BandsType,
          template <typename> class MatrixTemplate, typename ModelType,
          typename ConcurrencyType>
class susceptibility {
private:
  typedef Field FieldType;
  typedef std::complex<Field> ComplexType;
  typedef MatrixTemplate<Field> MatrixType;
  typedef MatrixTemplate<ComplexType> ComplexMatrixType;
  typedef std::vector<Field> VectorType;
  typedef std::vector<std::complex<Field>> ComplexVectorType;
  typedef PsimagLite::Range<ConcurrencyType> RangeType;
  typedef std::vector<SuscType> VectorSuscType;
#ifdef USE_SCGAP2D
  typedef rpa::gap2D<FieldType, psimag::Matrix, ConcurrencyType> GapType;
#elif USE_SCGAP3D
  // typedef rpa::gap3D<FieldType,psimag::Matrix,BandsType,ConcurrencyType>
  // GapType;
  typedef rpa::gap3D<FieldType, psimag::Matrix, ModelType, ConcurrencyType>
      GapType;
#else
  typedef rpa::gap2D<FieldType, psimag::Matrix, ConcurrencyType> GapType;
#endif

  rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &param;
  rpa::model<Field, MatrixTemplate, ConcurrencyType> &tbmodel;
  ConcurrencyType &conc;
  // momentumDomain<Field,psimag::Matrix> qMesh;
  size_t numberOfQ;
  size_t msize;
  size_t nq1;
  size_t nq2;
  size_t nq3;
  size_t nw;
  VectorType omega;
  std::vector<std::vector<FieldType>> QVec;
  std::vector<size_t> indexToiq;
  std::vector<size_t> indexToiw;
  FieldType wmin_, wmax_;
  bool writeFullChi0;
  bool kMap;

public:
  // typedef std::vector<SuscType> BaseType;

  susceptibility(
      rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &parameters,
      ModelType &modelIn, ConcurrencyType &concurrency, const FieldType &qxmin,
      const FieldType &qxmax, const size_t nq1In, const FieldType &qymin,
      const FieldType &qymax, const size_t nq2In, const FieldType &qzmin,
      const FieldType &qzmax, const size_t nq3In, const FieldType &wmin,
      const FieldType &wmax, const size_t nwIn)
      : param(parameters), tbmodel(modelIn), conc(concurrency),
        // qMesh(param,nqxIn,nqzIn,3),
        // numberOfQ(nq1In*nq2In*nq3In*nwIn),
        numberOfQ(1), msize(param.nOrb * param.nOrb), nq1(nq1In), nq2(nq2In),
        nq3(nq3In), nw(nwIn), omega(nw, 0),
        // QVec(numberOfQ,VectorType(4,0)),
        QVec(1, VectorType(4, 0)), indexToiq(numberOfQ, 0),
        indexToiw(numberOfQ, 0), wmin_(wmin), wmax_(wmax),
        writeFullChi0(param.writeFullChi0), kMap(0) {

    // Setup q-mesh based on nq's and q-limits
    setupQandOmegaMesh(nq1, nq2, nq3, numberOfQ, nw, qxmin, qxmax, qymin, qymax,
                       qzmin, qzmax, wmin, wmax, QVec);
    numberOfQ = QVec.size();

    if (param.scState == 1 && param.printGap == 1 && conc.rank() == 0) {
      if (conc.rank() == 0)
        std::cout << "Now writing gap.txt \n";
      printGap2();
      if (conc.rank() == 0)
        std::cout << "Done writing gap.txt \n";
    }

    if (param.Case == "Emery") {
      std::vector<ComplexMatrixType> chi0(numberOfQ, ComplexMatrixType(3, 3));
      std::vector<ComplexMatrixType> chi0_g(numberOfQ,
                                            ComplexMatrixType(19, 3));
      std::vector<ComplexMatrixType> chi0_gg(numberOfQ,
                                             ComplexMatrixType(19, 19));

      std::string filename("chi0Emery.txt");
      if (param.readChiForSus)
        readChi0Emery(QVec, chi0, chi0_g, chi0_gg, filename);
      else
        calcEmeryChi0(chi0, chi0_g, chi0_gg);
      calcEmeryRPA(chi0, chi0_g, chi0_gg);
    } else {
      VectorSuscType chi0Matrix(numberOfQ, SuscType(parameters, concurrency));
      calcElements(chi0Matrix);
      if (conc.rank() == 0) {
        std::cout << "Now printing out chiq \n";
        writeChiqTxt(chi0Matrix);
      }
    }
  }

  void calcElements(VectorSuscType &chi0Matrix) {
    // Setup k-mesh for chi0 calculation
    momentumDomain<Field, psimag::Matrix, ConcurrencyType> kmesh(
        param, conc, param.nkInt, param.nkIntz, param.dimension);
    kmesh.set_momenta(false);
    BandsType bands(
        param, tbmodel, conc, kmesh,
        param.cacheBands); // false = no Caching // true = Caching, needed here
                           // because we pre-calculate energies

    std::vector<FieldType> q(3);
    bool single_q = false;

    if (param.cacheBands) {
      // Pre-calculate band energies for k-mesh
      if (conc.rank() == 0)
        std::cout << "Pre-calculating e(k) \n";
      bands.precalculate_ekak(param.scState); // sets ek and ak
      bands.writeGap();
      if (nq1 * nq2 * nq3 == 1) { // only 1 q-vector --> Pre-calculate ekq,akq
        single_q = true;
        if (conc.rank() == 0)
          std::cout << "Pre-calculating e(k+q) \n";
        q[0] = QVec[0][0];
        q[1] = QVec[0][1];
        q[2] = QVec[0][2];
        bands.precalculate_ekqakq(q, param.scState);
      }
    }
    RangeType range(0, numberOfQ, conc);
    // SuscType chi0QW(param,conc);
    for (; !range.end(); range.next()) {

      size_t iQ = range.index();
      q[0] = QVec[iQ][0];
      q[1] = QVec[iQ][1];
      q[2] = QVec[iQ][2];
      if (!single_q && param.cacheBands) {
        bands.precalculate_ekqakq(
            q, param.scState); // resets ekq and akq for changing q vector
      }

      // if (param.scState==1 && !param.cacheBands) { // RPA/BCS calculation
      // with SC gap, band calculation on the fly 	GapType
      // Delta(param,conc);
      // 	calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
      //                calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],Delta,QVec[iQ][3],0);
      // } else {
      if (wmin_ == 0.0 && wmax_ == 0.0) { // only for zero frequency
        // if (param.cacheBands) { // Band energies and eigenvectors are
        // pre-calculated
        //    calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
        //        calcChi0(param,kmesh,bands,conc,chi0Matrix[iQ]);
        // } else {// Band energies and eigenvectors are calculated on the fly
        //    calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
        //        //
        //        calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],wmin_,kMap);
        //        calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],kMap);
        // }

        if (!param.scState) {
          calcChi0Matrix<FieldType, SuscType, BandsType, GapType,
                         MatrixTemplate, ConcurrencyType>
              calcChi0(param, kmesh, q, bands, conc, chi0Matrix[iQ],
                       param.cacheBands,
                       1); // Last parameter is kMap for calc. chi0k
        } else { // use finite w constructor for SC state calculation, but set
                 // w=0
          calcChi0Matrix<FieldType, SuscType, BandsType, GapType,
                         MatrixTemplate, ConcurrencyType>
              calcChi0(
                  param, kmesh, q, bands, conc, chi0Matrix[iQ], 0.0,
                  param.cacheBands); // Constructor for pre-calculated bands
        }
      } else { // Finite frequency calculation
        // if (param.cacheBands) { // Band energies and eigenvectors are
        // pre-calculated
        calcChi0Matrix<FieldType, SuscType, BandsType, GapType, MatrixTemplate,
                       ConcurrencyType>
            calcChi0(param, kmesh, q, bands, conc, chi0Matrix[iQ], QVec[iQ][3],
                     param.cacheBands); // Constructor for pre-calculated bands
        // } else { // Band energies and eigenvectors are calculated on the fly
        // calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
        //        calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],QVec[iQ][3],0);
        //        // Constructor for on the fly band diagonalization
        // }
      }
      // }

      if (conc.rank() == 0) {
        std::cout.precision(7);
        std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3]
                  << "  of "
                  << numberOfQ
                  // << " total. ChiPhys=" << chi0Matrix[iQ].calcSus()
                  << " total. ChiPhys=" << tbmodel.calcSus(chi0Matrix[iQ], "zz")
                  << tbmodel.calcSus(chi0Matrix[iQ], "+-")
                  << tbmodel.calcSus(chi0Matrix[iQ], "-+")
                  // << chi0Matrix[iQ].calcSus()
                  // << "chi0_{1133}" << chi0Matrix[iQ](0,18)
                  << "\n";
      }
    }

    for (size_t iq = 0; iq < numberOfQ; iq++)
      chi0Matrix[iq].allReduce();
  }

  void calcEmeryChi0(std::vector<ComplexMatrixType> &chi0,
                     std::vector<ComplexMatrixType> &chi0_g,
                     std::vector<ComplexMatrixType> &chi0_gg) {
    // Setup k-mesh for chi0 calculation
    momentumDomain<Field, psimag::Matrix, ConcurrencyType> kmesh(
        param, conc, param.nkInt, param.nkIntz, param.dimension);
    kmesh.set_momenta(false);
    BandsType bands(param, tbmodel, conc, kmesh, false); // false = no Caching
    RangeType range(0, numberOfQ, conc);

    for (; !range.end(); range.next()) {

      size_t iQ = range.index();
      std::vector<FieldType> q(3);
      q[0] = QVec[iQ][0];
      q[1] = QVec[iQ][1];
      q[2] = QVec[iQ][2];

      calcChi0Matrix<FieldType, SuscType, BandsType, GapType, MatrixTemplate,
                     ConcurrencyType>
          calcChi0(param, kmesh, bands, q, conc, chi0[iQ], chi0_g[iQ],
                   chi0_gg[iQ]);

      if (conc.rank() == 0)
        std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3]
                  << "  of " << numberOfQ << "\n";
    }

    for (size_t iq = 0; iq < numberOfQ; iq++) {
      conc.allReduce(chi0[iq]);
      conc.allReduce(chi0_g[iq]);
      conc.allReduce(chi0_gg[iq]);
    }

    if (conc.rank() == 0) {
      std::string filename("chi0Emery.txt");
      writeChi0Emery(QVec, chi0, chi0_g, chi0_gg, filename);
    }
  }

  void calcEmeryRPA(std::vector<ComplexMatrixType> &chi0,
                    std::vector<ComplexMatrixType> &chi0_g,
                    std::vector<ComplexMatrixType> &chi0_gg) {
    // Now calculate the RPA Chi
    interactionEmery<FieldType, psimag::Matrix, ConcurrencyType> rpaEmery(
        param);
    std::vector<ComplexMatrixType> chiS(numberOfQ, ComplexMatrixType(3, 3));
    std::vector<ComplexMatrixType> chiC(numberOfQ, ComplexMatrixType(3, 3));

    for (size_t iQ = 0; iQ < numberOfQ; iQ++) {
      // First calculate the effective interaction with chi0
      std::vector<FieldType> q(3);
      q[0] = QVec[iQ][0];
      q[1] = QVec[iQ][1];
      q[2] = QVec[iQ][2];
      ComplexMatrixType GammaS(19, 19);
      ComplexMatrixType GammaC(19, 19);
      ComplexMatrixType bareSpin(19, 19);
      ComplexMatrixType bareCharge(19, 19);
      ComplexMatrixType couplingSpin(19, 19);
      ComplexMatrixType couplingCharge(19, 19);
      rpaEmery.calcRPAResult(chi0_gg[iQ], 1, GammaS, bareSpin,
                             q); // renormalized spin interaction
      rpaEmery.calcRPAResult(chi0_gg[iQ], 0, GammaC, bareCharge,
                             q); // renormalized charge interaction
      // Now calculate the spin and charge susceptibilities
      ComplexMatrixType aGammaSa(3, 3), aGammaCa(3, 3);
      for (size_t l1 = 0; l1 < 3; l1++)
        for (size_t l2 = 0; l2 < 3; l2++) {
          for (size_t i = 0; i < 19; i++)
            for (size_t j = 0; j < 19; j++) {
              aGammaSa(l1, l2) +=
                  chi0_g[iQ](i, l1) * GammaS(i, j) * chi0_g[iQ](j, l2);
              aGammaCa(l1, l2) +=
                  chi0_g[iQ](i, l1) * GammaC(i, j) * chi0_g[iQ](j, l2);
            }
          chiS[iQ](l1, l2) = chi0[iQ](l1, l2) - aGammaSa(l1, l2);
          chiC[iQ](l1, l2) = chi0[iQ](l1, l2) - aGammaCa(l1, l2);
        }
      if (conc.rank() == 0) {
        std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3];
        ComplexType susS = calcChiSRPA(chiS[iQ]);
        ComplexType susN = calcChiNematicRPA(chiC[iQ]);
        std::cout << "spin sus.: " << susS << " , nematic sus.: " << susN
                  << "\n";
      }
    }
    writeChiSChiC(chiS, chiC);
  }

  ComplexType calcChiSRPA(ComplexMatrixType &chiS) {
    ComplexType sus(0.0);
    for (size_t i = 0; i < chiS.n_row(); i++)
      for (size_t j = 0; j < chiS.n_row(); j++) {
        sus += 0.5 * chiS(i, j);
      }
    return sus;
  }

  ComplexType calcChiNematicRPA(ComplexMatrixType &chiC) {
    ComplexType sus(0.0);
    sus = chiC(1, 1) - chiC(1, 2) + chiC(2, 2) - chiC(2, 1);
    return sus;
  }

  void setupQandOmegaMesh(size_t nq1, size_t nq2, size_t nq3, size_t numberOfQ,
                          size_t nw, const FieldType &qxmin,
                          const FieldType &qxmax, const FieldType &qymin,
                          const FieldType &qymax, const FieldType &qzmin,
                          const FieldType &qzmax, const FieldType &wmin,
                          const FieldType &wmax,
                          std::vector<std::vector<FieldType>> &QVec) {

    MatrixType momenta(0, 3);

    if (param.qGridType == "Path") {
      if (conc.rank() == 0)
        std::cout << "Setting up Q-Path along high sym. dir. \n";
      // Setup momenta along high-symmetry direction
      // std::string path("Path2");
      std::string path(param.momentumPath);
      size_t nkPath(nq1 * 3);
      if (param.momentumPath == "Path0_extended") {
        nkPath = nq1*5;
        }
      if (param.dimension == 3) {
        path = "Path3";
        nkPath = nq1 * 7;
      }
      rpa::momentumDomain<Field, psimag::Matrix, ConcurrencyType> kmesh2(
          param, conc, path, nkPath);
      numberOfQ = nkPath * nw;
      // std::cout << "nkPath: " << nkPath << "\n";
      momenta.resize(nkPath, 3);
      for (size_t ik = 0; ik < kmesh2.momenta.n_row(); ik++) {
        // if (conc.rank()==1) std::cout << kmesh2.momenta(ik,0) << "," <<
        // kmesh2.momenta(ik,1) << "," << kmesh2.momenta(ik,2) << "\n";
        momenta(ik, 0) = kmesh2.momenta(ik, 0);
        momenta(ik, 1) = kmesh2.momenta(ik, 1);
        momenta(ik, 2) = kmesh2.momenta(ik, 2);
      }

    } else {
      numberOfQ = nq1 * nq2 * nq3 * nw;
      if (conc.rank() == 0)
        std::cout << "Setting up regular Q-grid \n";
      momenta.resize(nq1 * nq2 * nq3, 3);
      if (nq1 == 1 && nq2 == 1 && nq3 == 1) { // Q is fixed
        momenta(0, 0) = qxmin;
        momenta(0, 1) = qymin;
        momenta(0, 2) = qzmin;
      } else if (nq1 > 1 && nq2 == 1 && nq3 == 1) { // 1D Q-Scan
        for (size_t iq1 = 0; iq1 < nq1; iq1++) {
          momenta(iq1, 0) =
              qxmin + float(iq1) / float(nq1 - 1) * (qxmax - qxmin);
          momenta(iq1, 1) =
              qymin + float(iq1) / float(nq1 - 1) * (qymax - qymin);
          momenta(iq1, 2) =
              qzmin + float(iq1) / float(nq1 - 1) * (qzmax - qzmin);
        }
      } else if (nq1 > 1 && nq2 > 1 && nq3 == 1) { // 2D Q-scan (in-plane)
        for (size_t iq1 = 0; iq1 < nq1; iq1++) {
          for (size_t iq2 = 0; iq2 < nq2; iq2++) {
            size_t index(iq2 + nq2 * iq1);
            momenta(index, 0) =
                qxmin + float(iq1) / float(nq1 - 1) * (qxmax - qxmin);
            momenta(index, 1) =
                qymin + float(iq2) / float(nq2 - 1) * (qymax - qymin);
            momenta(index, 2) = qzmin;
          }
        }
      } else if (nq1 > 1 && nq2 > 1 && nq3 > 1) { // 3D Q-scan
        for (size_t iq1 = 0; iq1 < nq1; iq1++) {
          for (size_t iq2 = 0; iq2 < nq2; iq2++) {
            for (size_t iq3 = 0; iq3 < nq3; iq3++) {
              size_t index(iq3 + nq3 * iq2 + nq3 * nq2 * iq1);
              momenta(index, 0) =
                  qxmin + float(iq1) / float(nq1 - 1) * (qxmax - qxmin);
              momenta(index, 1) =
                  qymin + float(iq2) / float(nq2 - 1) * (qymax - qymin);
              momenta(index, 2) =
                  qzmin + float(iq3) / float(nq3 - 1) * (qzmax - qzmin);
            }
          }
        }
      }
    }
    // Setup linear omega-mesh
    for (size_t i = 0; i < nw; i++)
      omega[i] =
          wmin + (wmax - wmin) * float(i) / fmax(float(nw - 1), float(1));
    // Now combine vectors
    size_t nQ = momenta.n_row();
    QVec.resize(numberOfQ, VectorType(4, 0));
    indexToiq.resize(numberOfQ);
    indexToiw.resize(numberOfQ);
    for (size_t iq = 0; iq < nQ; iq++)
      for (size_t iw = 0; iw < nw; iw++) {
        size_t index = iw + iq * nw;
        indexToiq[index] = iq;
        indexToiw[index] = iw;
        QVec[index][0] = momenta(iq, 0);
        QVec[index][1] = momenta(iq, 1);
        QVec[index][2] = momenta(iq, 2);
        QVec[index][3] = omega[iw];
      }
  }

  void writeChiqTxt(VectorSuscType &chi0Matrix) {
    int width(10);
    std::vector<FieldType> q(3);
    if (writeFullChi0) {
      std::string cstr = "chi0Full_" + param.fileID + ".txt";
      const char *filename = cstr.c_str();
      std::ofstream os(filename);
      os.precision(width);
      cstr = "chi1Full_" + param.fileID + ".txt";
      const char *filename1 = cstr.c_str();
      std::ofstream os1(filename1);
      os1.precision(width);
      os << std::fixed;
      os1 << std::fixed;
      os << "nq1,nq2,nq3,nw: \n";
      os1 << "nq1,nq2,nq3,nw: \n";
      os << nq1 << " , " << nq2 << " , " << nq3 << " , " << nw << "\n";
      os1 << nq1 << " , " << nq2 << " , " << nq3 << " , " << nw << "\n";
      SuscType chi1Matrix(param, conc);
      for (size_t iq = 0; iq < numberOfQ; iq++) {
        q[0] = QVec[iq][0];
        q[1] = QVec[iq][1];
        q[2] = QVec[iq][2];
        calcRPAResult(chi0Matrix[iq], tbmodel.spinMatrix, chi1Matrix, q);
        os << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3]
           << " , ";
        os1 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3]
           << " , ";
        for (size_t l1 = 0; l1 < msize; l1++)
          for (size_t l2 = 0; l2 < msize; l2++) {
            os << real(chi0Matrix[iq](l1, l2)) << " , ";
            os1 << real(chi1Matrix(l1, l2)) << " , ";
          }
        for (size_t l1 = 0; l1 < msize; l1++)
          for (size_t l2 = 0; l2 < msize; l2++) {
            os << imag(chi0Matrix[iq](l1, l2)) << " , ";
            os1 << imag(chi1Matrix(l1, l2)) << " , ";
          }
        ComplexType sus0(tbmodel.calcSus(chi0Matrix[iq], "zz"));
        ComplexType sus1(tbmodel.calcSus(chi0Matrix[iq], "+-"));
        ComplexType sus2(tbmodel.calcSus(chi0Matrix[iq], "-+"));
        os << real(sus0) << " , " << imag(sus0) << " , " << real(sus1) << " , "
           << imag(sus1) << " , " << real(sus2) << " , " << imag(sus2);
        os << "\n";
        sus0 = tbmodel.calcSus(chi1Matrix, "zz");
        sus1 = tbmodel.calcSus(chi1Matrix, "+-");
        sus2 = tbmodel.calcSus(chi1Matrix, "-+");
        os1 << real(sus0) << " , " << imag(sus0) << " , " << real(sus1) << " , "
           << imag(sus1) << " , " << real(sus2) << " , " << imag(sus2);
        os1 << "\n";
      }
    }
    std::string cstr = "chiRPA_" + param.fileID + ".txt";
    const char *filename = cstr.c_str();
    std::ofstream os2(filename);
    // interaction<FieldType,psimag::Matrix,ConcurrencyType> rpa(param,conc);
    os2.precision(width);
    os2 << std::fixed;
    SuscType chiRPA(param, conc);
//     os2 << "kx , ky, kz, w, ReChizz, ImChizz, ReChipm, ImChipm, ReChimp, ImChimp, ReChi0zz, ImChi0zz, ReChi0pm, ImChi0pm, "
//         << "ReChi0mp, ImChi0mp"
// #ifdef USE_1BANDALTERMAGNET
//         << ", ReChi"

    for (size_t iq = 0; iq < numberOfQ; iq++) {
      // std::cout << "iq:"<<iq<< " rank: " << conc.rank() << " sus: " <<
      // tbmodel.calcSus(chi0Matrix[0], "zz") << "\n";
      q[0] = QVec[iq][0];
      q[1] = QVec[iq][1];
      q[2] = QVec[iq][2];
      calcRPAResult(chi0Matrix[iq], tbmodel.spinMatrix, chiRPA, q);
      // std::cout << "RPA result: " << tbmodel.calcSus(chiRPA, "zz") << "\n";
      // std::cout << "spinMatrix: " << tbmodel.calcSus(tbmodel.spinMatrix,
      // "zz") << "\n"; std::cout << "chi0 result: " <<
      // tbmodel.calcSus(chi0Matrix[iq], "zz") << "\n";
      ComplexType susRzz(tbmodel.calcSus(chiRPA, "zz"));
      ComplexType susRpm(tbmodel.calcSus(chiRPA, "+-"));
      ComplexType susRmp(tbmodel.calcSus(chiRPA, "-+"));
      // ComplexType sus1(chi0Matrix[iq].calcSus());
      ComplexType sus1(tbmodel.calcSus(chi0Matrix[iq], "zz"));
      ComplexType sus2(tbmodel.calcSus(chi0Matrix[iq], "+-"));
      ComplexType sus3(tbmodel.calcSus(chi0Matrix[iq], "-+"));
      /* ComplexType sus3(tbmodel.calcSus(chi0Matrix[iq], "xx")); */
      /* ComplexType sus4(tbmodel.calcSus(chi0Matrix[iq], "yy")); */
#ifdef USE_1BANDALTERMAGNET
      ComplexType suspmAA(tbmodel.calcSusi1i2(chiRPA, "+-", 0, 0));
      ComplexType suspmBB(tbmodel.calcSusi1i2(chiRPA, "+-", 1, 1));
      ComplexType suspmAB(tbmodel.calcSusi1i2(chiRPA, "+-", 0, 1));
      ComplexType suspmBA(tbmodel.calcSusi1i2(chiRPA, "+-", 1, 0));
      ComplexType susmpAA(tbmodel.calcSusi1i2(chiRPA, "-+", 0, 0));
      ComplexType susmpBB(tbmodel.calcSusi1i2(chiRPA, "-+", 1, 1));
      ComplexType susmpAB(tbmodel.calcSusi1i2(chiRPA, "-+", 0, 1));
      ComplexType susmpBA(tbmodel.calcSusi1i2(chiRPA, "-+", 1, 0));
#endif
      os2 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3]
          << " , " << real(susRzz) << " , " << imag(susRzz)
          << " , " << real(susRpm) << " , " << imag(susRpm)
          << " , " << real(susRmp) << " , " << imag(susRmp)
          << " , " << real(sus1) << " , " << imag(sus1)
          << " , " << real(sus2) << " , " << imag(sus2)
          << " , " << real(sus3) << " , " << imag(sus3)
#ifdef USE_1BANDALTERMAGNET
          << " , " << real(suspmAA) << " , " << imag(suspmAA) // 16, 17
          << " , " << real(suspmBB) << " , " << imag(suspmBB) // 18, 19
          << " , " << real(suspmAB) << " , " << imag(suspmAB) // 20, 21
          << " , " << real(suspmBA) << " , " << imag(suspmBA) // 22, 23
          << " , " << real(susmpAA) << " , " << imag(susmpAA) // 24, 25
          << " , " << real(susmpBB) << " , " << imag(susmpBB) // 26, 27
          << " , " << real(susmpAB) << " , " << imag(susmpAB) // 28, 29
          << " , " << real(susmpBA) << " , " << imag(susmpBA) // 30, 31
#endif
          << "\n";
    }
    os2.close();
    std::cout << "File is written \n";
  }

  void writeChiSChiC(std::vector<ComplexMatrixType> &chiS,
                     std::vector<ComplexMatrixType> &chiN) {
    std::ofstream os1("chiSRPA.txt");
    std::ofstream os2("chiNRPA.txt");
    std::ofstream os3("chiCRPA.txt");
    int width(10);
    os1.precision(width);
    os2.precision(width);
    os3.precision(width);
    os1 << std::fixed;
    os2 << std::fixed;
    os3 << std::fixed;
    std::vector<FieldType> q(3);
    for (size_t iq = 0; iq < numberOfQ; iq++) {
      q[0] = QVec[iq][0];
      q[1] = QVec[iq][1];
      q[2] = QVec[iq][2];
      ComplexType chiRPAS = calcChiSRPA(chiS[iq]);
      ComplexType chiRPAN = calcChiNematicRPA(chiN[iq]);
      ComplexType chiRPAC = calcChiSRPA(chiN[iq]);
      os1 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3]
          << " , ";
      os1 << real(chiRPAS) << "," << imag(chiRPAS) << "\n";
      os2 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3]
          << " , ";
      os2 << real(chiRPAN) << "," << imag(chiRPAN) << "\n";
      os3 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3]
          << " , ";
      os3 << real(chiRPAC) << "," << imag(chiRPAC) << "\n";
    }
  }

  template <typename FieldType, typename ComplexMatrixType>
  void writeChi0Emery(const std::vector<std::vector<FieldType>> &qField,
                      const std::vector<ComplexMatrixType> &chi0,
                      const std::vector<ComplexMatrixType> &chi0_g,
                      const std::vector<ComplexMatrixType> &chi0_gg,
                      const std::string &filename) {
    if (conc.rank() == 0) {
      std::ofstream os(filename.c_str());
      int width(10);
      os.precision(width);
      os << std::fixed;
      size_t nq(qField.size());
      for (size_t iq = 0; iq < nq; iq++) {
        os << qField[iq][0] << " , " << qField[iq][1] << " , " << qField[iq][2]
           << " , " << qField[iq][3] << " , ";
        for (size_t i = 0; i < chi0[0].n_row(); i++)
          for (size_t j = 0; j < chi0[0].n_col(); j++)
            os << real(chi0[iq](i, j)) << " , " << imag(chi0[iq](i, j))
               << " , ";
        for (size_t i = 0; i < chi0_g[0].n_row(); i++)
          for (size_t j = 0; j < chi0_g[0].n_col(); j++)
            os << real(chi0_g[iq](i, j)) << " , " << imag(chi0_g[iq](i, j))
               << " , ";
        for (size_t i = 0; i < chi0_gg[0].n_row(); i++)
          for (size_t j = 0; j < chi0_gg[0].n_col(); j++) {
            os << real(chi0_gg[iq](i, j)) << " , " << imag(chi0_gg[iq](i, j));
            if (j < chi0_gg[0].n_col())
              os << " , ";
          }
        // os << "\n";
      }
    }
  }
  template <typename FieldType, typename ComplexMatrixType>
  void readChi0Emery(std::vector<std::vector<FieldType>> &qField,
                     std::vector<ComplexMatrixType> &chi0,
                     std::vector<ComplexMatrixType> &chi0_g,
                     std::vector<ComplexMatrixType> &chi0_gg,
                     std::string &filename) {
    size_t nq(0);
    if (conc.rank() == 0) {
      std::vector<FieldType> data;
      typedef std::complex<FieldType> ComplexType;
      loadVector(data, filename);
      size_t step(4 + 2 * (3 * 3 + 19 * 3 +
                           19 * 19)); // 4 qField (qx,qy,qz,w), 3x3 for chi0,
                                      // 19x3 for chi0_g, 19x19 for chi0_gg
      nq = data.size() / step;
      if (nq != numberOfQ) {
        qField.resize(nq);
        chi0.resize(nq);
        chi0_g.resize(nq);
        chi0_gg.resize(nq);
        for (size_t iq = 0; iq < nq; ++iq) {
          qField[iq] = std::vector<FieldType>(3, 0);
          chi0[iq] = ComplexMatrixType(3, 3);
          chi0_g[iq] = ComplexMatrixType(19, 3);
          chi0_gg[iq] = ComplexMatrixType(19, 19);
        }
      }
      // for (size_t iq=0;iq<nq;++iq) qField[iq] = std::vector<FieldType>(3,0);
      std::cout << "number of Q points from file: " << nq << "\n";
      for (size_t iq = 0; iq < nq; iq++) {
        qField[iq][0] = data[0 + iq * step];
        qField[iq][1] = data[1 + iq * step];
        qField[iq][2] = data[2 + iq * step];
        qField[iq][3] = data[3 + iq * step];
        size_t ishift(4 + iq * step);
        for (size_t i = 0; i < chi0[iq].n_row(); i++)
          for (size_t j = 0; j < chi0[iq].n_col(); j++) {
            FieldType r1 = data[ishift + 2 * (j + i * chi0[0].n_col())];
            FieldType i1 = data[ishift + 2 * (j + i * chi0[0].n_col()) + 1];
            chi0[iq](i, j) = ComplexType(r1, i1);
          }
        ishift = 4 + iq * step + 2 * chi0[0].n_row() * chi0[0].n_col();
        for (size_t i = 0; i < chi0_g[iq].n_row(); i++)
          for (size_t j = 0; j < chi0_g[iq].n_col(); j++) {
            FieldType r1 = data[ishift + 2 * (j + i * chi0_g[0].n_col())];
            FieldType i1 = data[ishift + 2 * (j + i * chi0_g[0].n_col()) + 1];
            chi0_g[iq](i, j) = ComplexType(r1, i1);
          }
        ishift = 4 + iq * step + 2 * chi0[0].n_row() * chi0[0].n_col() +
                 2 * chi0_g[0].n_row() * chi0_g[0].n_col();
        for (size_t i = 0; i < chi0_gg[iq].n_row(); i++)
          for (size_t j = 0; j < chi0_gg[iq].n_col(); j++) {
            FieldType r1 = data[ishift + 2 * (j + i * chi0_gg[0].n_col())];
            FieldType i1 = data[ishift + 2 * (j + i * chi0_gg[0].n_col()) + 1];
            chi0_gg[iq](i, j) = ComplexType(r1, i1);
          }
      }
    } else {
      conc.broadcast(nq);
      if (nq != numberOfQ) {
        qField.resize(nq);
        chi0.resize(nq);
        chi0_g.resize(nq);
        chi0_gg.resize(nq);
        for (size_t iq = 0; iq < nq; ++iq) {
          qField[iq] = std::vector<FieldType>(3, 0);
          chi0[iq] = ComplexMatrixType(3, 3);
          chi0_g[iq] = ComplexMatrixType(19, 3);
          chi0_gg[iq] = ComplexMatrixType(19, 19);
        }
      }
    }
    numberOfQ = nq;
    for (size_t iq = 0; iq < numberOfQ; iq++) {
      conc.broadcast(qField[iq]);
      conc.broadcast(chi0[iq]);
      conc.broadcast(chi0_g[iq]);
      conc.broadcast(chi0_gg[iq]);
    }
  }

  void printGap(const std::string &file = "gapRPA.txt") {
    VectorType data;
    loadVector(data, file);
    size_t step = 5;
    std::cout << "File=" << file << "\n";
    size_t nk(data.size() / step);
    std::cout << "number of k-points in gap file = " << nk << "\n";

    GapType Delta(param, conc);
    // momentumDomain<Field,psimag::Matrix>
    // kmesh(param,param.nkInt,param.nkIntz,param.dimension);
    // kmesh.set_momenta(false);
    // BandsType bands(param,conc,kmesh,true);
    // VectorType ek(param.nOrb,0), k(3);
    // ComplexMatrixType ak(param.nOrb,param.nOrb);
    std::ofstream os("gap.txt");
    int width(10);
    os.precision(width);
    os << std::fixed;
    for (size_t ik = 0; ik < nk; ++ik) {
      // kmesh.momenta.getRow(ik,k);
      // bands.getEkAndAk(k,ek,ak);
      VectorType k(3, 0);
      size_t iband;
      k[0] = data[ik * step];
      k[1] = data[ik * step + 1];
      k[2] = data[ik * step + 2];
      iband = size_t(data[ik * step + 3]);
      FieldType gapRPA(data[ik * step + 4]);
      // for (size_t iband=0;iband<param.nOrb;iband++) {
      ComplexType gap1 = Delta(k, iband);
      // gap1 *= pow(param.Omega0,2)/(pow(ek[iband],2)+pow(param.Omega0,2)); //
      // Lorentzian cut-off

      os << k[0] << " , " << k[1] << " , " << k[2] << " , " << iband << " , "
         << gapRPA << " , " << real(gap1) << "\n";
      std::cout << k[0] << " , " << k[1] << " , " << k[2] << " , " << iband
                << " , " << gapRPA << " , " << real(gap1) << "\n";
      // }
    }
  }

  void printGap2() {
    GapType Delta(param, tbmodel, conc);
    momentumDomain<Field, psimag::Matrix, ConcurrencyType> kmesh(
        param, conc, param.nkInt, param.nkIntz, param.dimension);
    kmesh.set_momenta(false);
    BandsType bands(param, tbmodel, conc, kmesh, true);
    VectorType ek(param.nOrb, 0), k(3);
    ComplexMatrixType ak(param.nOrb, param.nOrb);
    std::string file = "gapAll_" + param.fileID + ".txt";
    const char *filename = file.c_str();
    std::ofstream os(filename);
    int width(10);
    os.precision(width);
    os << std::fixed;
    for (size_t ik = 0; ik < kmesh.nktot; ++ik) {
      kmesh.momenta.getRow(ik, k);
      bands.getEkAndAk(k, ek, ak);
      std::vector<ComplexType> gap1(param.nOrb);
      for (size_t iband = 0; iband < param.nOrb; iband++) {
        gap1[iband] = Delta(k, iband, ak);
        // gap1[iband] *=
        // pow(param.Omega0,2)/(pow(ek[iband],2)+pow(param.Omega0,2)); //
        // Lorentzian cut-off
        gap1[iband] *= exp(-std::pow(ek[iband], 2) /
                           std::pow(param.Omega0, 2)); // Gaussian cut-off
      }

      os << k[0] << "  " << k[1] << "  " << k[2] << "  ";
      for (size_t ib = 0; ib < param.nOrb; ib++)
        os << real(gap1[ib]) << " " << imag(gap1[ib]) << " ";
      os << "\n";
    }
  }

  void printGap3() {
    ferminator<FieldType, BandsType, psimag::Matrix, ModelType, ConcurrencyType>
        FSpoints(param, tbmodel, conc, 1);
    size_t nk(FSpoints.nTotal);
    momentumDomain<Field, psimag::Matrix, ConcurrencyType> kmesh(
        param, conc, param.nkInt, param.nkIntz, param.dimension);
    BandsType bands(param, conc, kmesh, true);
    VectorType ek(param.nOrb, 0);
    ComplexMatrixType ak(param.nOrb, param.nOrb);
    GapType Delta(param, conc);
    std::ofstream os("gapOnFS.txt");
    int width(10);
    os.precision(width);
    os << std::fixed;
    for (size_t ik = 0; ik < nk; ++ik) {
      VectorType k(3, 0);
      size_t iband;
      k[0] = FSpoints.kFx[ik];
      k[1] = FSpoints.kFy[ik];
      k[2] = FSpoints.kFz[ik];
      iband = FSpoints.kFtoBand[ik];
      bands.getEkAndAk(k, ek, ak);
      ComplexType gap1 = Delta(k, iband, ak);
      // gap1 *= pow(param.Omega0,2)/(pow(ek[iband],2)+pow(param.Omega0,2)); //
      // Lorentzian cut-off
      os << k[0] << " , " << k[1] << " , " << k[2] << " , " << iband;
      for (size_t iorb = 0; iorb < param.nOrb; iorb++)
        os << " , " << FSpoints.weights[ik][iorb];
      os << " , " << real(gap1) << "\n";
    }
  }
};
} // namespace rpa

#endif
