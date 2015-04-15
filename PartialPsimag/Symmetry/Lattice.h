//-*-C++-*-

#ifndef PSIMAG_Lattice_H
#define PSIMAG_Lattice_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file Lattice.h  
 *
 *  Containe the base class definition for cell objects. */ 

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include "Vec.h"
#include "Mat.h"
#include "SeitzMatrix.h"
#include "Real.h"

#include "CellPosition.h"
#include "LatticeCoordinates.h"
#include "CellParameters.h"
#include "MetricTensor.h"
#include "Reciprocal.h"

#include "CellTranslation.h"
#include "Pattern.h"

#include "PSIMAGAssert.h"
#include "Matrix.h"

namespace psimag {

  template <typename Field, size_t DIM> class MetricTensor;
  template <typename Field, size_t DIM, typename Algorithms> class CellParameters;

  /** \ingroup latticeConcepts
   *
   * \brief The Lattice class inherits from SeitzMatrix<Field,DIM> the
   *        columns of this matrix are the Lattice's basis vectors.
   *
   * Apart from orientation a Lattice can be defined by \b CellParameters,
   * which include the length of the basis vectors and the angles
   * between them. 
   *  
   * The Lattice as a SeitzMatrix<Field, DIM> is an operator that:
   *  -  transforms from CellPositions to cartesian coordinates and 
   * Lattice->inverse()
   *  -  transforms from cartesian to (unNormalized) CellPositions.
   *
   * Cartesian positions can be 'normalized' or brought into the cell
   * formed by this lattice's basis vectors by:
   *
   * - transforming them to cell positions (multiply by *this)

   *
   * - normalizing the cell positions (modulus operation on the
   *   cellpositions components)
   * 
   * - transforming the 'normalized' cell positions back to caresian
   *   coordinates.
   *
   * The translation part of the Seitz matrix is generally not used in Lattices,
   * altough it is convenient to use the Sietz Matrix.
   *
   * \param Field: The scalar type used in the representation of the
   *               cell. Legal values include double, rational and
   *               sqrtExtendedRational.
   *
   * \param DIM: Dimensionality of the lattice being
   *             represented. Legal values are 1,2,3.
   */
  template <typename Field,size_t DIM,typename Algorithms>
  class Lattice: public SeitzMatrix<Field,DIM> {

  public:
    
    typedef CellParameters<Field,DIM,Algorithms>          CellParamType;
    typedef MetricTensor<Field,DIM>                       MetricType;
    typedef Lattice<Field,DIM,Algorithms>                 ThisType;
    typedef SeitzMatrix<Field,DIM>                        MatType;
    typedef typename MatType::Traits                      MatTraitsType;
    typedef CartesianTranslation<Field,DIM>               CartesianTranslationType;
    typedef CartesianPosition<Field,DIM>                  CartesianPositionType;
    typedef std::vector<CartesianTranslationType>         BasisVectors;
    typedef SymmetryOperation<Field,DIM,Algorithms>       SymmetryOperationType;
    typedef CellPosition<Field,DIM,Algorithms>            CellPositionType;
    typedef CellTranslation<Field,DIM>                    CellTranslationType;
    typedef LatticeTransformation<Field,DIM>              LatticeTransformationType;
    typedef LatticeCoordinates<DIM>                       LatticeCoordType;
    typedef std::vector<LatticeCoordType>                 LatticeCoordsType;

  public:

    static std::string typeName() { return "InputLattice";   }

    /** Default constructor. (Not for general use.)        **/
    Lattice(): 
      MatType(Field(0)),
      reciprocalSet(false),
      det(Field(0))
    {}

//     /**   
//      *  \brief It is not necessary to transform the pattern object to relect this lattice.
//      *
//      * \note This is here as an 'abstract' fulfillment of a
//      *       contract. Other Lattice classes will make transformations.
//      */
//     template<typename Occupant, size_t NUMPOS>
//     void updatePattern(Pattern<Field,DIM,NUMPOS,Occupant,Algorithms>& pattern) {
//       return;
//     }

  private:
    
    /** The inverse matrix of this Lattice.                    **/
    MatType                     inverseBasis;     

    /** Indicates whether or not the inverseBasis member is ready to use. **/
    bool                        inverseBasisSet;  

    /** A Lattice which is the reciprocal of this cell. **/
    //ThisType*                   reciprocal;       

  public:
    
    /**
     * \brief The transformation which was used to convert an origional Lattice to this lattice.
     */
    LatticeTransformationType   transform;
    
    /** Indicates whether or not the inverseBasis member is ready to use. **/
    bool                        reciprocalSet; // Should init to false.

    /** The metric tensor of the coordinate system defined by the basis of this Lattice. **/
    MetricType                  metric;         

    /**< The cell parameters (also called the lattice parameters) that correspond to this Lattice. **/  
    CellParamType               parameters;     

    /**< The determinent of this cell's matrix. **/
    Field                       det;            

    /** This Vec maps the Lattice's basis vector into a canonical order, (see 2002 ITC, pg 750) **/
    Vec<size_t, DIM>            order;           


    /**
     * Copy Construct a Lattice.
     */
    Lattice(const Lattice& lat):
      MatType(lat),
      inverseBasis(lat.inverseBasis),
      inverseBasisSet(lat.inverseBasisSet),
      //reciprocal(lat.reciprocal),
      transform(lat.transform),
      //reciprocalSet(lat.reciprocalSet),
      metric(lat.metric),
      parameters(lat.parameters),
      det(lat.det),
      order(lat.order)
    {}

    /**
     * The destructor deletes any reciprocal cell embedded in this cell and visa versa.
     */
    ~Lattice() {
//       if (reciprocalSet) {
// 	reciprocal->reciprocalSet = false; // avoid recursion
// 	delete reciprocal;
//       }
    }

    /**
     * Constructs a Lattice from CellParameters.
     */
    Lattice(const CellParamType& cellParameters): 
      MatType(Field(0)),
      //reciprocalSet(false),
      parameters(cellParameters),
      det(Field(0))
    {
      cellParameters.makeBasisFor(*this);
      bool skipParameters = true;
      update(skipParameters);
    }

    /**
     * Constructs a Lattice from a std::vector of CartesianTranslation vectors.
     */
    Lattice(const BasisVectors& translations):
      //reciprocalSet(false),
      det(Field(0))
    {
      assert(translations.size() == DIM);
      
      for (size_t r=0; r< DIM; r++) 
	for (size_t c=0; c< DIM; c++) 
	  (*this)(r,c) = translations[r][c];

      // This does not work because of dimensional mis match, 
      // *this is a 3x3 Seitz Matrix, translations is a 2x3 collection
      // Using Copy would require a different (specialized) traits object (later) &*&*&*
      //       COPY<MatType, 
      // 	BasisVectors, 
      // 	typename MatType::Traits, 
      // 	DoubleIndexTraits<Field,DIM,DIM> >::EXEC(*this, translations);
      
      update();
    }

    void getBasisVectors(std::vector<std::vector<Field> >& basisVectors) {

      basisVectors.resize(DIM);
      for (size_t r=0; r< DIM; r++) {
	basisVectors[r].resize(DIM);
	for (size_t c=0; c< DIM; c++) 
	  basisVectors[r][c] = (*this)(r,c);
      }
    }

    void getBasisVectors(Matrix<Field>& basisVectors) {
      
      basisVectors.resize(DIM,DIM);
      for (size_t r=0; r< DIM; r++) {
	for (size_t c=0; c< DIM; c++) 
	  basisVectors(r,c) = (*this)(r,c);
      }
    }

    /**
     * \brief Update the cells metric and parameters based on the
     *        cell's current set of basis vectors.  This method is
     *        called during initialization and during reduction.
     */
    void update(bool skipParameters=false) {

      setDet();

      metric = metricTensor(*this);

      order  = axisOrder(metric);

      if (!skipParameters)
	parameters = CellParamType(metric);
    }

    /**
     * Compute the determinent of this Lattice and throw an error if the value is zero.
     * This function is called within each Lattice constructor.
     */
    void setDet() {
      det = Det((MatType)*this);
      if ( det == Field(0) ) {
	std::ostringstream buff;
	buff << "ERROR Lattice(translations) resulted in a Det of 0" << std::endl;
	buff << ((MatType) *this) << std::endl;
	throw std::range_error(buff.str());
      }
    }

    /**
     * Return the ith basis vector of this cell, by canonical order.
     * (see 2002 ITC, pg 750)
     */
    CartesianTranslationType inOrder(size_t i) {
      assert(i < DIM);
      return (*this)[order[i]];
    }

    /**
     * Return 0 if the given cartesian symmetry operation maps the lattice onto
     * itself, otherwise return the "distance" of the mapped lattice
     * from this lattice.
     */
    Field symmetryDistance(const SymmetryOperation<Field,DIM,Algorithms>& cartSymmetryOperation) const {
      
      typedef typename MatType::Traits               MatTriats;
      typedef typename SymmetryOperationType::Traits SymTriats;
      
      MatType      mappedCartesianBasisVectors;
      MatType      mappedCellPositions;
      //      static bool trace(false);
      
      // Use the sym op and the (cartesian) basis vectors contained in
      // this lattice to generate mapped (cartesian) basis
      // vectors.
      Multiply(cartSymmetryOperation,(*this),mappedCartesianBasisVectors);
      
      // Convert the mapped basis vectors into the corresponding
      // (un-normalized) cell positions.
      Multiply(this->inverse(),mappedCartesianBasisVectors,mappedCellPositions);
      
      // Normalize the cell positions
      // Any basis vector that satifies the given symmetry operation will
      // normalize to zero.
      for (size_t pos=0; pos < DIM; pos++) 
	for (size_t i=0; i< DIM; i++) 
	  mappedCellPositions(i,pos) = Algorithms::modulus(mappedCellPositions(i,pos),Field(1));
      
      // Regenerate new cartesian positions that are within this
      // lattices basis.  Because we want the euclidean
      // distance. (Should we use the MetricTensor here instead?)
      Multiply(*this,mappedCellPositions,mappedCartesianBasisVectors);
      
      // Compute the 'distance' between this lattice and the mappedLattice.
      // If the lattice satisfies the symetry
      Field maxd(0);
      
      for (size_t pos=0; pos < DIM; pos++) {
	Field  d = 0; //distanceFromZero(pos); 
	for (size_t i=0; i< DIM; i++) {
	  Field di = mappedCartesianBasisVectors(i,pos);
	  di = di*di;
	  d += di;
	}
	if (d > maxd) maxd = d;
      }
      
      //return the distance
      return maxd;
    }

//     /* 
//      * Return a reference to a recipricol of this cell.  
//      *
//      * The recipricol is created by lazy evaluation.
//      *
//      *\note The cell and it's reciprocal are memory managed
//      *      together. When a cell is destroyed its reciprocal (if it
//      *      has one) is also destroyed.
//      */
//     ThisType& getReciprocal() {
      
//       if (reciprocalSet)
// 	return *reciprocal;

//       reciprocal = Reciprocal(this);
//       reciprocalSet = true;

//       reciprocal->reciprocal = this;
//       reciprocal->reciprocal.reciprocalSet = true;

//       return *reciprocal;
//     }

   /* 
     * Return the volume of the cell.
     * The cell's volume is given by it's determinant.
     */
    Field volume() const { return abs( det); }

   /* 
    * Return the CartesianPosition corresponding to the given CellPosition.
    */
    template<typename OtherField>
    CartesianPositionType cartesian(const CellPosition<OtherField,DIM,Algorithms>& pos) const {
      
      CartesianPosition< Field, DIM> result;
      Multiply((*this), pos, result);
      return result;
    }
    
   /* 
    * Return the CartesianPosition corresponding to the given CellPosition.
    */
    template<typename OtherField>
    CartesianPositionType cartesianPosition(const CellTranslation<OtherField,DIM>& trns) const {
      
      CartesianPosition< Field, DIM> result;
      Multiply((*this), trns, result);
      return result;
    }
    
   /* 
    * Return the CartesianPosition corresponding to the given CellPosition.
    */
    template<typename OtherField>
    CartesianPositionType cartesianPosition(const CellPosition<OtherField,DIM,Algorithms>& pos) const {
      
      CartesianPosition< Field, DIM> result;
      Multiply((*this), pos, result);
      return result;
    }
    
    /* 
    * Return the CartesianPosition corresponding to the given CellPosition.
    */
    CartesianTranslationType cartesianTranslation(const CellTranslationType& trns) const {
      CartesianTranslationType result;
      Multiply((*this), trns, result);
      return result;
    }

    /* 
    * Return the CartesianPosition corresponding to the given CellPosition.
    */
    CartesianTranslationType cartesianTranslation(const CellPositionType& pos) const {
      CartesianTranslationType result;
      Multiply((*this), pos, result);
      result[DIM] = Field(0);
      return result;
    }
    
   /* 
     * Return the CartesianTranslation corresponding to the given LatticeCoordinates.
     */
    CartesianTranslationType cartesianTranslation(const LatticeCoordinates<DIM>& latCoord) const {
      CartesianTranslation<Field,DIM> result;
      for (size_t i=0; i<DIM; i++) 
	result += (*this)[i] * Field(latCoord[i]);
      return result;
    }

    /* 
     * Return the CellTranslation corresponding to the given LatticeCoordinates and factor.
     */
    CellTranslationType cellTranslation(const LatticeCoordinates<DIM>& latCoord, Field factor) const {

      // Lattice Coordinates are CellPositions whose values are restricted to integers.
      CellTranslation<Field,DIM> result = latCoord;
      result *= factor;
      return result;
    }

    /* 
     * Return the CellTranslation corresponding to the given Cartesian Translation.
     */
    CellTranslationType cellTranslation(const CartesianTranslationType& t) const {
      
      CellTranslationType result;
      Multiply<Field>(this->inverse(), t, result);
      return result;
    }

    /* 
     *
     */
    bool cellContains(const CellPositionType& cellPos) const {
      for (size_t i = 0; i< DIM; i++) {


	if (cellPos[i] >= Field(1)) return false;

	if (Algorithms::close(cellPos[i], Field(0)) ) {

	  continue;
	}
	if (cellPos[i] <  Field(0)) return false;
      }

      return true;
    }

    /* 
     *
     */
    bool cellContains(const CartesianPositionType& cartPos) const {
      return cellContains(cellPosition(cartPos));
    }

   /* 
    * Return the CellPosition corresponding to the given CartesianPosition.
    */
    CellPosition<Field,DIM,Algorithms> cellPosition(const CartesianPositionType& pos) const {
      
      CellPosition< Field, DIM ,Algorithms> result;
      Multiply(this->inverse(), pos, result);
      return result;
    }

    /* 
     * Returns CartesianTranslation corresponding to the i'th vector of this Cell's basis.
     *
     */
    CartesianTranslationType operator [] (size_t j) const {
      
      CartesianTranslationType trans;
      for(size_t i=0; i < DIM; i++)
	trans[i] = (*this)(i,j);
      return trans;
    }


    ThisType subLattice(const LatticeCoordsType& coords) const {
      
      LatticeTransformation<Field,DIM>  aTransform;
      for(size_t i=0; i < DIM; i++) 
	for(size_t j=0; j < DIM; j++) 
	  aTransform(j,i) = coords[i][j];
      aTransform.finalize();

      ThisType result;
      aTransform(*this, result);
      return result;
    }

//     template<typename OtherField>
//     CartesianTranslationType translationTo( const CellPosition<OtherField,DIM,Algorithms>& coord ) const {
      
//       CartesianTranslation< Field, DIM > result;
//       Multiply<Field>((*this), coord, result);
//       return result;
//     }
 



    //---------------------------------------------------------------------- SymmetryOperationType

    /** 
     * \brief Transform the given cartesian Symmetry operation into a cellbased symmetry operation. 
     *        This is a Similarity transform. cell-> cart->cart->cell.
     *
     */
    void cellOperation(const SymmetryOperation<Field,DIM,Algorithms>& cartesianSymop,  
		       SymmetryOperation<Field,DIM,Algorithms>& cellSymop) const
    {
      typedef typename SymmetryOperation<Field,DIM,Algorithms>::MatType MatType;
      MatType temp;
      Multiply(cartesianSymop, (*this), temp);
      Multiply(this->inverse(), temp, cellSymop);
    }

  };

  //====================================================================== 
  
  /** \ingroup xml
   *
   * XML Output function for 2D CellParameters.
   */
  template<typename Field, size_t DIM,typename Algorithms>
  Tag toXML(const Lattice<Field,DIM,Algorithms>& lat,
	    std::string name="Lattice") {
      
    typedef typename Lattice<Field,DIM,Algorithms>::MatType MatType; 
    typedef typename MatType::Traits                        MatTraitsType;

    Tag result(name);
    result["det"]    = lat.det;
    result["order"]  = lat.order;
    result["type"]   = lat.parameters.typeString();
    result.add(toXML(lat.parameters));
    result.add(toXML(lat.metric));
    result.add(toXML(lat.transform));
    
    std::ostringstream buff;
    MAT_PRINT<MatType,MatTraitsType>::JUST_NUMBERS( lat, buff);
    result.content << buff.str();
    
    return result;
  }
  
  //====================================================================== 
  /** \ingroup ostream
   * Lattice output stream operator 
   **/
  template<typename Field, size_t DIM,typename Algorithms>
  std::ostream& operator << (std::ostream& os, const Lattice<Field,DIM,Algorithms>& cell) {

    typedef typename Lattice<Field,DIM,Algorithms>::MatType LatticeMatType;
    
    os.setf(std::ios_base::fixed, std::ios_base::floatfield); os.precision(6);
    
    os << "Input Lattice { det: "<< cell.det << ", "
       << "order: [" << cell.order << "], "
       << "params : "<< cell.parameters << ", " 
       << "type: " << cell.parameters.typeString() << "}" << std::endl;
    os << " ---------------------------------- Basis:"  << std::endl;
    os << ( (LatticeMatType) cell ) << std::endl;
    os << " ---------------------------------- End Basis:"  << std::endl;
    os << cell.metric << std::endl;
    os << "----------------------------------- End Lattice " << std::endl;
    
    return os;
  }

} /* namespace psimag */

#endif

/*@}*/

