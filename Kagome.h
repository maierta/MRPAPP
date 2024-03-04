// Model file for  single band model
#ifndef KAGOME_H
#define KAGOME_H

#include <cstdlib> // for atof and atoi
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Matrix.h"
#include "Range.h"
#include "math.h"
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
  const ComplexType ii;

public:
  FieldType nbands;
  ComplexMatrixType spinMatrix;
  ComplexMatrixType chargeMatrix;

  model(
      const rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &parameters,
      ConcurrencyType &concurrency)
      : param(parameters), conc(concurrency), dim(param.dimension),
        ii(ComplexType(0.0, 1.0)), nbands(param.nOrb),
        spinMatrix(nbands * nbands, nbands * nbands),
        chargeMatrix(nbands * nbands, nbands * nbands) {
    std::cout << "Kagome Hubbard Model \n";
    if (nbands != 3)
      std::cerr << "Number of orbitals should be 3! Exiting ...\n";

    msize = 9;
    setupInteractionMatrix();
  }

  inline void getBands(const VectorType k, VectorType &eigenvals,
                       ComplexMatrixType &eigenvects) {
    FieldType t = param.hopping_t;
    // const ComplexType ii = ComplexType(0.0,1.0);

    eigenvects(0, 0) = -param.mu;
    eigenvects(1, 1) = -param.mu;
    eigenvects(2, 2) = -param.mu;
    eigenvects(0, 1) = -2. * t * cos(0.5 * k[0]);
    eigenvects(0, 2) = -2. * t * cos(0.25 * k[0] + 0.25 * sqrt(3.) * k[1]);
    eigenvects(1, 2) = -2. * t * cos(0.25 * k[0] - 0.25 * sqrt(3.) * k[1]);
    // eigenvects(0,1) = -t * (1. + exp(-ii*k[0]));
    // eigenvects(0,2) = -t * (1. + exp(-ii*k[0]/2.)*exp(-ii*k[1]*sqrt(3.)/2.));
    // eigenvects(1,2) = -t * (1. + exp(+ii*k[0]/2.)*exp(-ii*k[1]*sqrt(3.)/2.));
    eigenvects(1, 0) = conj(eigenvects(0, 1));
    eigenvects(2, 0) = conj(eigenvects(0, 2));
    eigenvects(2, 1) = conj(eigenvects(1, 2));

    eigen(eigenvals, eigenvects);
  }

  void setupInteractionMatrix() {
    size_t nOrb(param.nOrb);
    FieldType U(param.U);

    for (size_t l1 = 0; l1 < nOrb; ++l1) {
      size_t ind1 = l1 + l1 * nOrb;
      spinMatrix(ind1, ind1) = U + param.deltaU[l1];
      chargeMatrix(ind1, ind1) = -U - param.deltaU[l1];
    }
  }

  std::complex<Field> calcSus(const ComplexMatrixType &sus,
                              const std::string &component = "zz") const {
    std::complex<Field> chiPhys(0.0);
    for (size_t l1 = 0; l1 < param.nOrb; ++l1) {
      for (size_t l2 = 0; l2 < param.nOrb; ++l2) {
        size_t ind1(l1 + l1 * param.nOrb);
        size_t ind2(l2 + l2 * param.nOrb);
        chiPhys += 0.5 * sus(ind1, ind2);
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
