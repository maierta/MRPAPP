// Model file for 1-band altermagnet Lieb model with explicit spin --> 4 bands
// total from PRL 134, 096703 (2025)
#ifndef SINGLEBAND_ALTERMAGNET_LIEB_WSPIN_H
#define SINGLEBAND_ALTERMAGNET_LIEB_WSPIN_H

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
  std::vector<int> spinOfEll;
  std::vector<size_t> orbOfEll;

  model(
      const rpa::parameters<Field, MatrixTemplate, ConcurrencyType> &parameters,
      ConcurrencyType &concurrency)
      : param(parameters), conc(concurrency), dim(param.dimension),
        nbands(param.nOrb), spinMatrix(nbands * nbands, nbands * nbands),
        chargeMatrix(nbands * nbands, nbands * nbands), spinOfEll(nbands),
        orbOfEll(nbands) {
    if (nbands != 4)
      std::cerr
          << "Number of orbitals should be 4 for this model! Exiting ...\n";

    msize = nbands * nbands;

    // Note that here the "orbitals" denote a combined orbital and spin index
    // The basis is (A, up; B, up; A, down; B, down)
    spinOfEll[0] = +1;
    orbOfEll[0] = 0; // spin up, A-sublattice
    spinOfEll[1] = +1;
    orbOfEll[1] = 1; // spin up, B-sublattice
    spinOfEll[2] = -1;
    orbOfEll[2] = 0; // spin dn, A-sublattice
    spinOfEll[3] = -1;
    orbOfEll[3] = 1; // spin dn, B-sublattice

    setupInteractionMatrix();
  }

  inline void getBands(const VectorType k, VectorType &eigenvals,
                       ComplexMatrixType &eigenvects) {
    FieldType t, t1, t2, h, cx, cy, cxy;

    t = param.hopping_t;
    t1 = param.hopping_t1;
    t2 = param.hopping_t2;

    h = param.zeemanField; // Zeeman field

    cx = cos(k[0]);
    cy = cos(k[1]);
    cxy = cos(0.5 * k[0]) * cos(0.5 * k[1]);

    FieldType ekAB = -4 * t * cxy;
    FieldType ekAA = -2 * t1 * cx - 2 * t2 * cy;
    FieldType ekBB = -2 * t2 * cx - 2 * t1 * cy;

    // Write Hamiltonian into eigenvects
    // Basis is (A, up; B, up; A, down; B, down)
    // Note that the Hamiltonian is block-diagonal in the spin

    ComplexMatrixType temp(4, 4);
    VectorType evals(4);

    for (size_t i = 0; i < nbands; i++)
      for (size_t j = 0; j < nbands; j++)
        eigenvects(i, j) = ComplexType(0., 0.);

    for (size_t i = 0; i < nbands; i++)
      for (size_t j = 0; j < nbands; j++)
        temp(i, j) = ComplexType(0., 0.);

    temp(0, 0) = -param.mu + h + ekAA;
    temp(1, 1) = -param.mu - h + ekBB;
    temp(2, 2) = -param.mu - h + ekAA;
    temp(3, 3) = -param.mu + h + ekBB;

    temp(0, 1) = ekAB;
    temp(1, 0) = ekAB;
    temp(2, 3) = ekAB;
    temp(3, 2) = ekAB;

    eigen(evals, temp);

    for (size_t b = 0; b < 4; b++) {
      eigenvals[b] =
          evals[b]; // 0,1 are for bonding (A+B), 2,3 for anti-bonding (A-B)
      for (size_t l = 0; l < 4; l++) {
        eigenvects(l, b) = temp(l, b);
      }
    }
    //

    // Only diagonalize spin up block
    // Spin down block can be generated with A, up = B, down and vice versa
    //
    /* ComplexMatrixType temp(2, 2); */
    /* VectorType evals(2); */
    /**/
    /* for (size_t i = 0; i < 2; i++) */
    /*   for (size_t j = 0; j < 2; j++) */
    /*     eigenvects(i, j) = ComplexType(0., 0.); */
    /**/
    /* for (size_t i = 0; i < 2; i++) */
    /*   for (size_t j = 0; j < 2; j++) */
    /*     temp(i, j) = ComplexType(0., 0.); */
    /**/
    /* temp(0, 0) = -param.mu + h + ekAA; */
    /* temp(1, 1) = -param.mu - h + ekAA; */
    /**/
    /* temp(0, 1) = ekAB; */
    /* temp(1, 0) = ekAB; */
    /**/
    /* eigen(evals, temp); */
    /**/
    /* for (size_t b = 0; b < 2; b++) { */
    /*   eigenvals[b] = evals[b]; // 0 is bonding (A+B), spin up, 1 is
     * anti-bonding */
    /*                            // (A-B), spin down */
    /*   eigenvals[b + 2] = */
    /*       evals[b]; // 2 is for bonding spin down, 3 for anti-bonding spin
     * down */
    /*   for (size_t l = 0; l < 2; l++) { */
    /*     eigenvects(l, b) = temp(l, b); */
    /*   } */
    /*   eigenvects(2, b + 2) = eigenvects(1, b); // A, down == B, up */
    /*   eigenvects(3, b + 2) = eigenvects(0, b); // B, down == A, up */
    /* } */
  }

  void setupInteractionMatrix() {
    FieldType U(param.U);
    int s1, s2, s3, s4, o1, o2, o3, o4, l1, l2, l3, l4;

    for (size_t ind1 = 0; ind1 < msize; ind1++) {
      for (size_t ind2 = 0; ind2 < msize; ind2++) {
        spinMatrix(ind1, ind2) = 0.0;
        l1 = param.indexToOrb(ind1, 0);
        l2 = param.indexToOrb(ind1, 1);
        l3 = param.indexToOrb(ind2, 0);
        l4 = param.indexToOrb(ind2, 1);
        s1 = spinOfEll[l1];
        s2 = spinOfEll[l2];
        s3 = spinOfEll[l3];
        s4 = spinOfEll[l4];
        o1 = orbOfEll[l1];
        o2 = orbOfEll[l2];
        o3 = orbOfEll[l3];
        o4 = orbOfEll[l4];

        // U-terms
        if (o1 == o2 && o1 == o3 && o1 == o4) {
          if (s1 == -s2 && s1 == s3 && s1 == -s4)
            spinMatrix(ind1, ind2) += U;
          if (s1 == s2 && s1 == -s3 && s1 == -s4)
            spinMatrix(ind1, ind2) -= U;
        }
      }
    }
  }

  std::complex<Field> calcSus(const ComplexMatrixType &sus,
                              const std::string &component = "zz") const {
    std::complex<Field> chiPhys(0.0);
    int s1, s2, s3, s4, o1, o2, o3, o4, l1, l2, l3, l4;
    for (size_t ind1 = 0; ind1 < msize; ind1++) {
      for (size_t ind2 = 0; ind2 < msize; ind2++) {
        l1 = param.indexToOrb(ind1, 0);
        l2 = param.indexToOrb(ind1, 1);
        l3 = param.indexToOrb(ind2, 0);
        l4 = param.indexToOrb(ind2, 1);
        s1 = spinOfEll[l1];
        s2 = spinOfEll[l2];
        s3 = spinOfEll[l3];
        s4 = spinOfEll[l4];
        o1 = orbOfEll[l1];
        o2 = orbOfEll[l2];
        o3 = orbOfEll[l3];
        o4 = orbOfEll[l4];

        if (o1 == o2 && o3 == o4) {
          if (component == "zz" && s1 == s2 && s3 == s4) {
            chiPhys +=
                sus(ind1, ind2) * ComplexType(s1, 0) * ComplexType(s3, 0);
          } else if (component == "+-" && s1 == +1 && s3 == +1 && s2 == -1 &&
                     s4 == -1) { // for +- - susc.
            chiPhys += sus(ind1, ind2);
          } else if (component == "-+" && s1 == -1 && s3 == -1 && s2 == +1 &&
                     s4 == +1) { // for -+ - susc.
            chiPhys += sus(ind1, ind2);
          }
        }
      }
    }
    return 0.125 * chiPhys;
  }

  std::complex<Field> calcSusi1i2(const ComplexMatrixType &sus,
                                  const std::string &component = "zz",
                                  size_t i1 = 0, size_t i2 = 0) const {
    std::complex<Field> chiPhys(0.0);
    int s1, s2, s3, s4, o1, o2, o3, o4, l1, l2, l3, l4;
    for (size_t ind1 = 0; ind1 < msize; ind1++) {
      for (size_t ind2 = 0; ind2 < msize; ind2++) {
        l1 = param.indexToOrb(ind1, 0);
        l2 = param.indexToOrb(ind1, 1);
        l3 = param.indexToOrb(ind2, 0);
        l4 = param.indexToOrb(ind2, 1);
        s1 = spinOfEll[l1];
        s2 = spinOfEll[l2];
        s3 = spinOfEll[l3];
        s4 = spinOfEll[l4];
        o1 = orbOfEll[l1];
        o2 = orbOfEll[l2];
        o3 = orbOfEll[l3];
        o4 = orbOfEll[l4];

        if (o1 == o2 && o3 == o4 && o1 == i1 && o3 == i2) {
          if (component == "zz" && s1 == s2 && s3 == s4) {
            chiPhys +=
                sus(ind1, ind2) * ComplexType(s1, 0) * ComplexType(s3, 0);
          } else if (component == "+-" && s1 == +1 && s3 == +1 && s2 == -1 &&
                     s4 == -1) { // for +- - susc.
            chiPhys += sus(ind1, ind2);
          } else if (component == "-+" && s1 == -1 && s3 == -1 && s2 == +1 &&
                     s4 == +1) { // for -+ - susc.
            chiPhys += sus(ind1, ind2);
          }
        }
      }
    }
    return 0.125 * chiPhys;
  }

  std::complex<Field> calcSCGap(
      VectorType &k, size_t band,
      ComplexMatrixType &Uk) { // dummy function; handled by gaps3g.h directly
    return ComplexType(0., 0.);
  }
};
} // namespace rpa

#endif
