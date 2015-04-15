//-*-C++-*-

#ifndef  Psimag_HermiteNormalForm
#define  Psimag_HermiteNormalForm

/** \ingroup extendedMatrices */
/*@{*/

/** \file HermiteNormalForm.h
 *
 *  Contains: 
 *
 *    - a class for representing the integer Hermite Normal Form (HNF)
 *      of a given matrix.
 *
 *    - a class for solving (integer) equations represented by a HNF
 *      and a given vector.
 */
 
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

#include "Tag.h"
#include "Mat.h"
#include "Vec.h"
#include "rational.h"
#include "PSIMAGAssert.h"

namespace psimag {
  
  /** \ingroup extendedMatrices 
   *
   *\brief Objects of this class are matrices in hermite normal forms
   *       of the matrix they were constructed from.
   *
   * HermiteNormal Form: (google Hermite Normal Form)
   *
   * The code is a restructured version of the code used in the cctbx:
   *
   *     cctbx_sources/scitbx/include/scitbx/matrix/row_echelon.h
   * 
   * This code in turn was derived from the CrystGAP package (GAP
   * Version 3.4.4).
   *
   * <see>
   *     B. Eick, F. Ga"hler and W. Nickel
   *     Computing Maximal Subgroups and Wyckoff Positions of Space Groups
   *     Acta Cryst. (1997). A53, 467 - 474
   * </see>
   * 
   */
  template <typename Field, size_t NROW, size_t NCOL=NROW, typename Traits=ColMajorTraits<Field,NROW,NCOL> > 
  class HermiteNormalForm: public Mat< Field, NROW, NCOL, Traits> {

  public:

    typedef HermiteNormalForm<Field, NROW, NCOL, Traits> ThisType;
    typedef Mat<Field, NROW, NCOL, Traits>               BaseType;
    typedef Mat<Field, NCOL, NROW, Traits>               TransformType;
    
  protected:
    
    size_t _rank;
    
  public:

    size_t rank() const { return _rank;}
    
    /**
     * \brief A copy of the matrix this object was
     *        constructed from.
     */
    BaseType origional;

    /**
     * \brief The assocaited transformation matrix. T*origional = this.
     */
    TransformType T;
    
    /**
     * \brief Default constructor.
     *
     *        T*origional = this
     */
    HermiteNormalForm():
      BaseType(),
      _rank(0),
      origional()
    {MakeIdentity(T);}

 //    /**
//      * \brief Copy Construct
//      */
//     HermiteNormalForm(const HermiteNormalForm& other):
//       BaseType(other),
//       _rank(other.rank)
//       //      origional(other.origional),
//       //      T(other.T)
//     {
//       Copy(origional, other.origional);
//       copy(T, other.T);
//     }

    /**
     * \brief Construct a Mat< Field, NROW, NCOL > which is the
     *        hermite normal form (also called the row echelon form)
     *        of the given origional matrix, also construct and
     *        contain the assocaited transformation matrix, T.
     *
     *        T*origional = this
     */
    HermiteNormalForm(const BaseType& origional):
      BaseType(origional),
      _rank(0),
      origional(origional)
    {
      MakeIdentity(T);
    }

    /** 
     * \brief Exchange the values in row r1 with those in row r2.
     */
    void swap2Rows(size_t r1, size_t r2) {

      this->swapRows(r1,r2);
      T.swapRows(r1,r2);

      //display();
    }

    /** 
     * \brief Make the values in the specified row the negative of
     *        their origional values.
     */
    void negateRow(size_t row) {

      for(size_t col=0; col < NCOL; col++) 
	(*this)(row,col) *= -1;

      for(size_t col=0; col < NCOL; col++) 
	T(row,col) *= -1;

      //display();
    }

    void display() {
      check();
    }

    void check() {

      bool ok = true; 

      if (!(abs(Det(T)) == 1 )) {
	ok = false;
      }

      Mat<Field,NROW,NCOL> temp;
      Multiply(T, origional, temp);
      if (!(*this == temp  )) {
	ok = false;
      }
      
    }
    
    /** 
     * \brief Make row r1 equal to that row minus a times row 2.
     *
     */
    void combineRows(Field a, size_t r1, size_t r2) {

      if (a == 0) return;

      for(size_t c=0; c < NCOL; c++) 
	(*this)(r1,c) -= a * (*this)(r2,c);

      for(size_t c=0; c < NCOL; c++) 
	T(r1,c) -= a * T(r2,c);

      //display();
    }
    
    /** 
     * \brief Compute the values of the HNF matrix from the values of
     *        the matrix this object was constructed from.
     */
    void normalize() {
      
      size_t count = 0;
      size_t row = 0;
      size_t col = 0;
      
      while( row < NROW && col < NCOL ) {

	size_t row2 = row; 
	
	// First row2 with non-zero val in col
	while (row2 < NROW && (*this)(row2,col) == 0) 
	  row2++;

	// all rows below row had 0 col: go to next col
	if (row2 == NROW) {

	  col++;
	  continue;
	}

	// row had 0 for col but row2 did not swap em
	if (row2 != row) {

	  this->swap2Rows(row,row2);
	}

	// row now has a non-zero in col, row2 may have a zero
	for (row2++; row2 < NROW; row2++) {

	  if ((*this)(row2, col) == 0) continue;
	  
	  if (abs((*this)(row2, col)) < abs((*this)(row,col))) {

	    this->swap2Rows(row,row2);
	  }
	}

	if ((*this)(row,col) < 0) {
	  negateRow(row);
	}

	bool cleared = true;
	for (row2 = row+1; row2 < NROW; row2++) {

	  Field a = (*this)(row2,col) / (*this)(row,col);

	  combineRows(a, row2, row);

	  if ((*this)(row2,col) != 0) 
	    cleared = false;
	}

	if (cleared) {
	  row++; col++; 
	}
      
      }
      //M = scitbx::mat_ref<Field>(M.begin(), i, MC);
      _rank = row;
      //      check();
    }

    /**
     * \brief Return the index of the first non-zero column in the
     *        given row.
     */
    size_t pivotColumn(size_t row) const {
      for(size_t col=0; col<NCOL; col++)
	if ((*this)(row,col)) 
	  return col; 
      throw -1;
    }
  };

  /** \ingroup extendedMatrices 
   *
   *\brief Objects of this class are used to find the solution to
   *       equations represented a given integer HNF matrix and a
   *       given integer vector.
   *
   * If scaled_integer_solution = true is used in the construction of
   * BackSubstitution, the solution will be an integer vector which is
   * scaled by the object's scaleFactor member (also an integer).
   *
   * <see> "Algorithms for deriving crystallographic space-group
   *   information", R.W. Grosse-Kunstleve, March 1999.
   * </see>
   * 
   * \param SolutionType: The type of the solution vector.  
   *
   * The code is a restructured version of the code used in the cctbx:
   *
   *     cctbx_sources/scitbx/include/scitbx/matrix/row_echelon.h
   * 
   * This code in turn was derived from the CrystGAP package (GAP
   * Version 3.4.4).
   *
   */
  template <typename Field, 
	    size_t NROW, size_t NCOL,
	    typename SolutionType=Field> 
  class BackSubstitution {
    
  private:
    
    void initFlags() {
      if (flag_indep[0]) 
	for(size_t c=0; c < NCOL; c++) 
	  flag_indep[c] = true;
    }
    
  public:
    
    const HermiteNormalForm<Field,NROW,NCOL>& m;
    const Vec<SolutionType, NROW>&              v;

    Vec<SolutionType, NCOL>        sol;
    Vec<bool,         NCOL>        flag_indep;
    Field                        scaleFactor;

    bool computeSolution;
    bool scaled_integer_solution;
    bool computeIndependence;
    
    /**
     * \brief Objects of this class represent the solution to a system
     *        of equations specified by the given HermiteNormalForm
     *        object and a Vec of values.
     *
     * \param solution: compute the solution, turn this off if
     *        you are just interested in the independence by row of
     *        the hnf matrix.
     *
     * \param scaled: perform any scaling necessary
     *        so that the solution vector is an iteger vector
     *        (whatever it's SolutionType). The scaling is kept in the
     *        scaleFactor member.
     *
     * \param independence: Turn this off if you are just
     *        interested in the independence by row of the hnf matrix.
     *
     */
    BackSubstitution(HermiteNormalForm<Field,NROW,NCOL> hnf,
		     Vec< Field, NCOL > values, 
		     bool solution         = true,
		     bool scaled           = false,
		     bool independence     = false):
      computeSolution(solution),
      scaled_integer_solution(scaled),
      computeIndependence(independence),
      m(hnf),
      v(values), 
      scaleFactor(1)
    {
      if (computeIndependence)
	initFlags();
      
      // Iterate over the rows in reverse order
      for (size_t r = NROW; r > 0;) {
	r--;
	try {
	  size_t pivotC = m.pivotColumn(r);

	  if (computeIndependence) 
	    flag_indep[pivotC] = false;
	  if (computeSolution)     
	    set_sol(r, pivotC);
	  continue;
	}
	catch (int) {
	  if (v[r] != 0) 
	    throw std::range_error("BackSubstitution: row with no pivot was non-zero");

	  if (computeSolution)
	    sol[r] = 1;  // arbitrary, make it small but non-zero
	}
      }
    }
    
    /**
     * \brief Scale up the current scaleFactor by the given factor,
     *        adjusting the solution vector as appropriate.
     *
     */
    void rescale(Field newFactor) {

      if (newFactor == Field(1)) return;
      scaleFactor *= newFactor;
      for (size_t i=0; i< NCOL; i++)
	sol[i] *= newFactor;
    }

    /**
     * \brief Compute the solution value for the row r, given the
     *        pivot column of that row.
     *
     * If scaled_integer_solution is true, compute a new scale factor
     * and rescale if necessary. 
     *
     * \note This assumes that the solution for all rows greater than
     *       the given row have already been computed.
     *
     */
    void set_sol(size_t  r, size_t  pivotC) {
      
      Field sum = Field(0);
      for (size_t jc= pivotC; jc < NCOL; jc++) 
	sum += m(r,jc) * sol[jc];
	  
      Field num = (v[r]  - sum);
      Field den = m(r,pivotC);

      if (scaled_integer_solution) {

	Field gcd_ = gcd(num, den);
	if (den < 0) 
	  gcd_ *= Field(-1);
	Field newfactor = den / gcd_;
	rescale(newfactor);
	sol[r] = num/gcd_;
	return;
      }

      sol[r] = (SolutionType) num / (SolutionType) den;
    }

  };

  //====================================================================== 
  /** \ingroup XML
   *
   * SpaceGroup  XML output
   **/
  template<typename Field, size_t NROW, size_t NCOL>
  Tag toXML(const HermiteNormalForm<Field, NROW, NCOL>& hnf,
	    std::string name="HermiteNormalForm")
  {
    typedef HermiteNormalForm<Field,NROW,NCOL>       HermiteNormalFormType;
    typedef typename HermiteNormalFormType::BaseType MatType;
    typedef typename MatType::Traits                 MatTraitsType;

    Tag result(name);

    result["rank"] = hnf.rank();

    Tag transformTag("Transformation");
    std::ostringstream tbuff;
    MAT_PRINT<MatType,MatTraitsType>::JUST_NUMBERS( hnf.T, tbuff);
    transformTag.content << tbuff.str();
    
    result.add(transformTag);
    
    std::ostringstream buff;
    MAT_PRINT<MatType,MatTraitsType>::JUST_NUMBERS( hnf, buff);
    result.content << buff.str();
    
    return result;
  }
  
} /* namespace psimag */

#endif //Psimag_HermiteNormalForm
  /*@}*/
