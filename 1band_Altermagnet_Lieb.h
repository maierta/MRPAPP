// Model file for 1-band altermagnet Lieb model
// from PRL 134, 096703 (2025)
//
#ifndef SINGLEBAND_ALTERMAGNET_LIEB_H
#define SINGLEBAND_ALTERMAGNET_LIEB_H

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

namespace rpa {

template <typename Field, template <typename> class MatrixTemplate,
          typename ConcurrencyType>
class model {
private:
  typedef MatrixTemplate<Field> MatrixType;
  typedef std::complex<Field> ComplexType;
  typedef MatrixTemplate<ComplexType> ComplexMatrixType;
  typedef std::vector<Field> VectorType;
  typedef Field FieldType;
  const rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &param;
  ConcurrencyType &conc;
  size_t dim;

private:
  size_t msize;

public:
  FieldType nbands;
  ComplexMatrixType spinMatrix;
  ComplexMatrixType chargeMatrix;
  std::vector<size_t> orbOfEll;

  model(
      const rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &parameters,
      ConcurrencyType &concurrency)
      : param(parameters), conc(concurrency), dim(param.dimension),
        nbands(param.nOrb), spinMatrix(nbands * nbands, nbands * nbands),
        chargeMatrix(nbands * nbands, nbands * nbands), orbOfEll(nbands) {
    if (nbands != 2)
      std::cerr
          << "Number of orbitals should be 2 for this model! Exiting ...\n";

    msize = nbands * nbands;

    orbOfEll[0] = 0; // A-sublattice
    orbOfEll[1] = 1; // B-sublattice

    setupInteractionMatrix();
  }

  inline void getBands(const VectorType k, VectorType &eigenvals,
                       ComplexMatrixType &eigenvects) {
    FieldType t, ta, tb, cx, cy, cxy;

    t = param.hopping_t;   // inter-sub-lattice hopping
    ta = param.hopping_t1; // hopping ta
    tb = param.hopping_t2; // hopping tb

    cx = cos(k[0]);
    cy = cos(k[1]);
    cxy = cos(0.5 * k[0]) * cos(0.5 * k[1]);

    FieldType ekAB = -4 * t * cxy;
    FieldType ekAA = -2 * ta * cx - 2 * tb * cy;
    FieldType ekBB = -2 * tb * cx - 2 * ta * cy;

    // Write Hamiltonian into eigenvects
    // Basis is (A, B)

    ComplexMatrixType temp(2, 2);
    VectorType evals(2);

    for (size_t i = 0; i < nbands; i++)
      for (size_t j = 0; j < nbands; j++)
        eigenvects(i, j) = ComplexType(0., 0.);

    for (size_t i = 0; i < nbands; i++)
      for (size_t j = 0; j < nbands; j++)
        temp(i, j) = ComplexType(0., 0.);

    temp(0, 0) = -param.mu + ekAA;
    temp(1, 1) = -param.mu + ekBB;

    temp(0, 1) = ekAB;
    temp(1, 0) = ekAB;

    eigen(evals, temp);

    for (size_t b = 0; b < 2; b++) {
      eigenvals[b] = evals[b];
      for (size_t l = 0; l < 2; l++) {
        eigenvects(l, b) = temp(l, b);
      }
    }
    //
  }

  void setupInteractionMatrix() {
    FieldType U(param.U);
    int o1, o2, o3, o4, l1, l2, l3, l4;

    for (size_t ind1 = 0; ind1 < msize; ind1++) {
      for (size_t ind2 = 0; ind2 < msize; ind2++) {
        spinMatrix(ind1, ind2) = 0.0;
        l1 = param.indexToOrb(ind1, 0);
        l2 = param.indexToOrb(ind1, 1);
        l3 = param.indexToOrb(ind2, 0);
        l4 = param.indexToOrb(ind2, 1);
        o1 = orbOfEll[l1];
        o2 = orbOfEll[l2];
        o3 = orbOfEll[l3];
        o4 = orbOfEll[l4];

        // U-terms
        if (o1 == o2 && o1 == o3 && o1 == o4) {
          spinMatrix(ind1, ind2) = U;
        }
      }
    }
  }

  std::complex<Field> calcSus(const ComplexMatrixType &sus,
                              const std::string &component = "zz") const {
    std::complex<Field> chiPhys(0.0);
    int o1, o2, o3, o4, l1, l2, l3, l4;
    for (size_t ind1 = 0; ind1 < msize; ind1++) {
      for (size_t ind2 = 0; ind2 < msize; ind2++) {
        l1 = param.indexToOrb(ind1, 0);
        l2 = param.indexToOrb(ind1, 1);
        l3 = param.indexToOrb(ind2, 0);
        l4 = param.indexToOrb(ind2, 1);
        o1 = orbOfEll[l1];
        o2 = orbOfEll[l2];
        o3 = orbOfEll[l3];
        o4 = orbOfEll[l4];

        if (o1 == o2 && o3 == o4) {
          chiPhys += 0.25 * sus(ind1, ind2);
        }
      }
    }
    return chiPhys;
  }

  std::complex<Field> calcSCGap(
      VectorType &k, size_t band,
      ComplexMatrixType &Uk) { // dummy function; handled by gaps3g.h directly
    return ComplexType(0., 0.);
  }
};
} // namespace rpa

#endif
