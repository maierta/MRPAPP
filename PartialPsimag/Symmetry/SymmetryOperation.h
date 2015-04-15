//-*-C++-*-

#ifndef  Psimag_SYMMETRY_OPERATION
#define  Psimag_SYMMETRY_OPERATION

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SymmetryOperation.h
 *
 *  Contains a class for implementing symmetry operations.
 *
 *  \author Mike Summers
 *
 */
 
#include <iostream>
#include <sstream>
#include <string>

#include "Vec.h"
#include "Mat.h"
#include "PSIMAGAssert.h"

#include "FieldConvert.h"

#include "SeitzMatrix.h"
#include "SeitzVectors.h"
#include "SeitzPosition.h"
#include "SeitzTranslation.h"
#include "CellPosition.h"
#include "LatticeCoordinates.h"
#include "CellTranslation.h"
#include "CellRotation.h"
#include "SymmetryElement.h"
#include "SymmetryElementName.h"
#include "MatMagnitude.h"
#include "HermiteNormalForm.h"

namespace psimag {
  

  //template<typename Field, size_t DIM> class SeitzMatrix;
  template<typename Field, size_t DIM> class CellRotation;

  /** \ingroup symmetryConcepts
   *
   * \brief A class for implementing symmetry operations.
   *
   * If the dimension of the space is DIM the SymmetryOperations are
   * DIM+1 by DIM+1 Seitz (or augmented) matricies. (c.f. 2002 ITC page
   * 721). These can be used to represent rotations, reflections, inversions,
   * rotoinversions, translations, screw rotations, and glide reflections
   *
   * These matricies work together with seitz positions and translation vectors.
   *
   * \notes DIM template parameter added to provide support for
   *        3D and 2D periodicity. Requires representation 
   *        of Zero and One and also the negative operator (-)
   *        from class Field.  
   */

  template<typename Field, size_t DIM, typename Algorithms>
  class SymmetryOperation: public SeitzMatrix< Field, DIM > {
    
  public:

    typedef SeitzMatrix< Field, DIM >               SuperClass;   
    typedef SuperClass                              MatType;
    typedef typename MatType::RotationType          RotationType;

    typedef SymmetryOperation<Field,DIM,Algorithms> SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>   SymmetryElementType;
    typedef SeitzTranslation<Field,DIM>             SeitzTranslationType;
    typedef SeitzPosition<Field,DIM>                SeitzPositionType;
    typedef LatticeCoordinates<DIM>                 LatticeCoordinatesType;
    typedef HermiteNormalForm<Field, DIM+1, DIM+1>  HermiteNormalFormType;

    std::string             name;

    //====================================================================== Constructors
    /**
     *\brief  Change this symmetryOperation to the inversion operation
     */
    void setToInversion() {
      size_t i,j;
      for(i=0;i<DIM;++i) 
	for(j=0;j<DIM;++j) 
	  if (i==j)
	    if (i==DIM-1)
	      (*this)(i,j) = Field(1); 
	    else
	      (*this)(i,j) = Field(-1); 
	  else
	    (*this)(i,j) = Field(0); 
      //name = "inversion";
    }

    //====================================================================== Constructors
    /**
     * \brief  The default constructor sets to identity
     */
    SymmetryOperation(): 
      SuperClass(),
      name("un-named")
    {}

    /**
     * \brief  The copy constructor 
     */
    SymmetryOperation(const SymmetryOperationType& other): 
      SuperClass(other), 
      name(other.name)
    {}

    /**
     * \brief Construct from a SeitzTranslation.
     */
    template<typename IN_TYPE>
    SymmetryOperation(const SeitzVector<IN_TYPE,DIM,0>& t): 
      SuperClass(t),
      name("un-named(Translation)")
    {}

    /**
     * \brief  Construct from a DIM*Dim length single array.
     *
     * This sets only the rotation part of this SietzMatrix.
     *
     * Convert data types as necessary.
     */
    template<typename IN_TYPE>
    SymmetryOperation(const IN_TYPE data [DIM*DIM]) : 
      SuperClass(data),
      name("un-named")
    {check();}

    //template<typename T> class TComp {public: typedef const T ArrayType[DIM*DIM]; typedef const T VecType[DIM]; };

    /**
     * \brief  Construct from a DIM*Dim length single array.
     *
     * This sets only the rotation part of this SietzMatrix.
     *
     * Convert data types as necessary.
     */
    template<typename IN_TYPE>
    SymmetryOperation(const IN_TYPE data [DIM*DIM], const IN_TYPE tdata [DIM]) : 
      SuperClass(data, tdata),
      name("un-named")
    {check();}


    HermiteNormalFormType hermiteNormalForm() const {
      const MatType& thisMat = *this;
      Mat<Field,DIM+1,DIM+1> mat;
      Copy(mat,thisMat);
      HermiteNormalFormType result(mat);
      return result;
    }

    /**
     * \brief  Make sure we are constructing a unimodular matrix.
     * 
     * \throw Throw a range_error if the determinant is not +/- 1.
     */
    void check() {
      Field det = Det((*this).rotation);
      if ( det == Field(0) ) {
	//if ( det != Field(1) and det != Field(-1) ) {

	std::string message("ERROR SymmetryOperation: ");
	message += " Tried to construct a SymmetryOperation with a Det != +/- 1";
	throw std::range_error(message);
      }
      return;
      Field tr  = Trace((*this).rotation);
      if ( tr < Field(-3) or tr > Field(3) ) {

	std::string message("ERROR SymmetryOperation: ");
	message += " Tried to construct a SymmetryOperation with a Trace out of range";
	throw std::range_error(message);
      }
      
    }

    SymmetryOperation& operator = (const SymmetryOperation& op) {

      SeitzMatrix<Field, DIM>&       tgtMatrix = (*this);
      const SeitzMatrix<Field, DIM>& srcMatrix = op;

      name            = op.name;
      tgtMatrix       = srcMatrix;

      return (*this);
    }
    
    //====================================================================== Classification

    void normalize()  {
      for (size_t i=0; i< DIM; i++) 
	(*this)(i,DIM) = Algorithms::modulus((*this)(i,DIM),Field(1));
    }
    
    /**
     * \brief Return an integer (Rotation::RotationType) representing
     *        the type of rotation contained in this SymmetryOperation.
     *
     */
    typename CellRotation<Field, DIM>::RotationType type() {
      return (*this).rotation.type();
    }

    /**
     * \brief Return the a (canonical) point of this operation or throw
     *        an error.
     *
     *   The fixed point for:
     * 
     *    - the identity -> the origin
     *
     *    - a mirror     -> the a (or b) axis intercept
     *
     *    - a glide      -< the corresponding mirror 'fixed'point' 
     *
     * \note only 2D rotations have fixed points This should be else where
     */
    SeitzPositionType fixedPoint() {

      static SeitzPosition<Field,DIM>  zero;
      static SymmetryOperationType     ident;
      SeitzPositionType                fixedPoint;

      MatType      symOp_ident         = (*this)-ident;
      Field        detSymOp_ident      = Det(symOp_ident);

      if (Algorithms::close(detSymOp_ident,Field(0))) {
	throw std::range_error("SymmetryOperation Fixed Point");
	//later return a fixed point which is the cannonical offset
	// for the mirror or glide.
	//	return fixedPoint;
      }
      MatType      symOp_ident_inverse = symOp_ident.inverse();
      Multiply(symOp_ident_inverse, zero, fixedPoint);
      return fixedPoint;
    }

    /**
     * \brief Sets the translation part so that the operation is
     *        has the given point as a fixed point.
     *
     * \note This does not apply to glides.
     */
    void setFixedPoint(const SeitzPositionType& pos) {
      static RotationType identity;
      MakeIdentity(identity);
      RotationType i_r = identity - this->rotation;
      Multiply(i_r,pos,this->translation);
      //SeitzPositionType fPoint = fixedPoint();
      //      if (! (fPoint == pos) ) {
      // 	std::ostringstream buff;
      // 	buff << std::endl << "setFixedPoint(" << pos << "): " << std::endl;
      // 	buff << "      identity           = " << identity << std::endl;
      // 	buff << "      rotation           = " << this->rotation << std::endl;
      // 	buff << "      identity -rotation = " << i_r << std::endl;
      // 	buff << "      trans              = " << this->translation << std::endl;
      // 	buff << "      fPoint = " << fPoint << std::endl;
      // 	throw std::logic_error(buff.str());
      //       }
    }
    
    /**
     * \brief Return the trace of this operation
     *
     */
    Field trace() const {
      return Trace(this->rotation);
    }

    /**
     * \brief Return the glide vector for this operation.
     *
     * \note This assumes that the operation is either 
     *       a translation, a mirror or a glide.
     */
    SeitzTranslationType translationPart() const {
      static Field two(2);
      static Field zero(0);
      MatType sq;
      Multiply(*this,*this,sq);
      SeitzTranslationType result( sq.translation / two);
      for (size_t i=0; i< DIM; i++)
	if (Algorithms::close(result[i],zero)) result[i] = zero;
      return result;
    }
    
    /**
     * \brief Return the sign of the angle that this operation rotates
     *        about the givenAxis at the fixed point.
     *
     */
    Field sinTheta(const SeitzPositionType&        fixedPoint,
		   const SeitzTranslation<Field,3> givenAxis) const {

      SeitzPositionType testPoint(fixedPoint);
      testPoint[0] = testPoint[0] + Field(1);

      SeitzTranslationType testDir = testPoint - fixedPoint;  

      SeitzPositionType movedPoint;
      Multiply((*this), testPoint, movedPoint);

      SeitzTranslationType movedPointDir = movedPoint - fixedPoint;

      SeitzTranslation<Field,3> cross = testDir % movedPointDir;

      Field cnorm = cross * givenAxis;

      Field result = cnorm/(testDir.l2norm() * movedPointDir.l2norm());

      return result;
    }

    /**
     * \brief Return the trace of this operation
     *
     */
    int sense(const SeitzPositionType& fixedPoint,
	      const SeitzTranslation<Field,3> givenAxis) const {

      Field st = sinTheta(fixedPoint,givenAxis);
      int result = Algorithms::sign(st);

      return result;
    }

  };

  //======================================================================
  /**
   * \brief Set values in the given symmetry element to reflect the given operation.
   *
   * \note This assumes that the operation type is either mirror or
   *       glide and that the translation part has already been set.
   */
  template<typename Field, typename LatticeType, typename Algorithms>
  void setRotationElements(SymmetryElement<Field,2,Algorithms>&         element,
			   const SymmetryOperation<Field,2,Algorithms>& cartesianOp,
			   const LatticeType& lat) {
    
    enum {DIM=2};

    typedef SymmetryOperation<Field,DIM,Algorithms> SymmetryOperationType;
    typedef CellPosition<Field,DIM,Algorithms>      CellPositionType;
    typedef CartesianPosition<Field,DIM>            CartesianPositionType;

    static Field                     zero(0);
    static Field                     minusTwo(-2);
    static Field                     minusOne(-1);
    static Field                     one(1);

    SymmetryOperationType            ident;
    CellPositionType                 zeroPos;
    static CellTranslation<Field,3>  zaxis     = cellTranslation<Field>(0,0,1);
    Field                            trace     = element.trace;
    //    SymmetryOperationType            cellOp;
    CartesianPositionType            cartesianPosition;

    // Use the lattice to translate to cell coordinates!
    //lat.cellOperation(cartesianOp, cellOp);

    // Set the cell position to the fixed-point
    Multiply((cartesianOp-ident).inverse(), zeroPos, cartesianPosition);
    element.cellPosition = lat.cellPosition(cartesianPosition);
    
    if (Algorithms::close(trace,minusTwo)) {
      element.type = "twoFold";
      //      element.rotationAngle = Field(180);
    }
    if (Algorithms::close(trace , minusOne)) {
      element.type = "threeFold";
      //      element.rotationAngle = Field(120);
    }
    if (Algorithms::close(trace , zero)) {
      element.type = "fourFold";
      //      element.rotationAngle = Field(90);
    }
    if (Algorithms::close(trace , one)) {
      element.type = "sixFold";
      //      element.rotationAngle = Field(60);
    }

    Field sinTheta_       = cartesianOp.sinTheta(cartesianPosition,zaxis);
    element.rotationAngle = degrees(asin(sinTheta_)); 
    element.sense         = Algorithms::sign(sinTheta_);

    if (element.sense == -1)
      element.type = element.type + "N";

  }

  //======================================================================
  /**
   * \brief Uses the given cell operation to set the element's
   *        position and netDirection.
   *
   * \note This assumes that the operation type is either mirror or
   *       glide.
   */
  template<typename Field, typename Algorithms>
  void setCanonicalPositionAndNetDirection(SymmetryElement<Field,2,Algorithms>&        element,
					   const SymmetryOperation<Field,2,Algorithms> cellOp,
					   bool doNetDirection=true)
  {
    enum        {DIM=2};
    typedef CellPosition<Field,DIM,Algorithms>  CellPositionType;

    CellPositionType a = seitzPosition<Field,Field>(1,0);
    CellPositionType movedA;
    Multiply(cellOp,a,movedA);
    
    try {    // Make some data points so we can find the slope and intercepts

      CellPositionType origin;
      CellPositionType movedOrigin;
      Multiply(cellOp,origin,movedOrigin);
      
      element.cellPosition = mirrorLineIntercept(origin,movedOrigin,a,movedA);
      
      if (doNetDirection)
	element.netDirection = getNetDirectionFor(origin,movedOrigin,a,movedA);
    }
    catch (std::logic_error e) {  // try some other points
      
      CellPositionType b = seitzPosition<Field,Field>(0,1);
      CellPositionType movedB;
      Multiply(cellOp,b,movedB);
      
      element.cellPosition = mirrorLineIntercept(a,movedA,b,movedB);
      
      if (doNetDirection)
	element.netDirection = getNetDirectionFor(a,movedA,b,movedB);

      // If this catch fails again a logic_error will be thrown
    }
  }

  //======================================================================
  /**
   * \brief Return the x intercept of the mirror line assocaited with
   *        the cell operation that produced the given set of points.
   *
   * \note This assumes that the operation type producing the points
   *       is either a 2D mirror or a 2D glide.
   */
  template<typename Field, typename Algorithms>
  bool acceptableInterceptPt(const CellPosition<Field,2,Algorithms>&  p)
  {
    static Field zero(0);
    // On the x axis
    if (Algorithms::close(p[1],zero)) return true;
    // Not on the x axis but on the y axis
    if (Algorithms::close(p[0],zero)) return true;
    return false;
  }
  
  //======================================================================
  /**
   * \brief Return the x intercept of the mirror line assocaited with
   *        the cell operation that produced the given set of points.
   *
   * \note This assumes that the operation type producing the points
   *       is either a 2D mirror or a 2D glide.
   */
  template<typename Field, typename Algorithms>
  const CellPosition<Field,2,Algorithms> 
  mirrorLineIntercept(const CellPosition<Field,2,Algorithms>&  pt1,
		      const CellPosition<Field,2,Algorithms>&  movedPt1,
		      const CellPosition<Field,2,Algorithms>&  pt2,
		      const CellPosition<Field,2,Algorithms>&  movedPt2) 
  {
    enum {DIM=2};
    typedef CellPosition<Field,DIM,Algorithms> CellPositionType;
    
    // If pt1 did not move then it is on the mirrorLine 
    if (pt1 == movedPt1 && acceptableInterceptPt(pt1)) return pt1;
    
    // The points just swapped
    if (movedPt1 == pt2 && movedPt2 == pt1) {
      const CellPositionType& biSect = bisector(pt1, pt2);
      if ( acceptableInterceptPt(biSect) )
	return biSect;
    }
    
    // get the intercept from the bisectors of the given points
    const CellPositionType p1(bisector(pt1, movedPt1 ));
    const CellPositionType p2(bisector(pt2, movedPt2 ));
    
    // If the bisectors are the same point
    if (p1.closeTo(p2)) {
      if (acceptableInterceptPt(p1))
	return p1;
      // The bisectors cannot be used to find the intercept
      throw std::logic_error("mirrorLineIntercept could not find a position");
    }
    
    return intercept<Field,Algorithms>(p1,p2);
  }

  //======================================================================
  /**
   * \brief Return the cell slope of the mirror line assocaited with
   *        the given cell operation.
   *
   * \note This assumes that the operation type is either mirror or
   *       glide.
   */
  template<typename Field, typename Algorithms>
  const LatticeCoordinates<2> getNetDirectionFor(const CellPosition<Field,2,Algorithms>&  pt1,
						 const CellPosition<Field,2,Algorithms>&  movedPt1,
						 const CellPosition<Field,2,Algorithms>&  pt2,
						 const CellPosition<Field,2,Algorithms>&  movedPt2) 
  {
    enum {DIM=2};
    typedef CellPosition<Field,DIM,Algorithms> CellPositionType;
    
    CellPositionType p1(bisector(pt1, movedPt1 ));
    CellPositionType p2(bisector(pt2, movedPt2 ));

    // If the bisectors are the same point
    if (p1.closeTo(p2)) 
      throw std::logic_error("mirrorLineXIntercept could not find a position");
    
    try {
      Field slp = slope<Field,Algorithms>(p2,p1);
      return latticeCoordinate<Field,Algorithms>(slp);
    }
    catch (std::range_error& e) {
      // vertical
      return latticeCoordinate(0,1);
    }
  }
  
  //======================================================================
  /**
   * \brief Set values in the given symmetry element to reflect the given operation.
   *
   * \note This assumes that the operation type is either mirror or
   *       glide and that the translation part has already been set.
   *
   * \note We are copying cartesianOp in so we can modify a temp of it.
   */
  template<typename Field, typename LatticeType, typename Algorithms>
  void setMirrorElements(SymmetryElement<Field,2,Algorithms>&  element,
			 SymmetryOperation<Field,2,Algorithms> cartesianOp,
			 const LatticeType& lat) {
    

    enum {DIM=2};

    typedef SymmetryOperation<Field,DIM,Algorithms>                  SymmetryOperationType;
    typedef typename SymmetryOperationType::MatType                  MatType;
    typedef typename MatType::TranslationType                        TranslationType;

    typedef CellPosition<Field,DIM,Algorithms>                       CellPositionType;
    typedef CellTranslation<Field,DIM>                               CellTranslationType;

    static Field               two(2);
    static Field               zero(0);
    CellTranslation<Field,DIM> zeroDiff;
    SymmetryOperationType      cellOp;

    element.sense = 0;

    // Use the lattice to translate to cell coordinates!
    lat.cellOperation(cartesianOp, cellOp);

    element.translationPart   = cellOp.translationPart();

    // Set the element type
    bool zeroTranslation = closeTo<Field,DIM,0,Algorithms>(element.translationPart,zeroDiff);

    if(element.trace == two) {
      // Set identity
      if (zeroTranslation) {
	element.type          = "identity";
	element.glideFraction = zero;
	return;
      }
      else {
	element.type          = "translation";
	element.netDirection  = latticeCoordinate<Field,Algorithms>(element.translationPart);
	element.setGlideFraction();
	return;
      }
    }

    if (zeroTranslation)  {
      element.type = "mirror";
      setCanonicalPositionAndNetDirection<Field,Algorithms>(element,cellOp);
      element.setGlideFraction();
      return;
    }
    else { 
      element.type = "glide";

      // Just use the mirror part of the glide operation to find the position
      COPY<TranslationType,
	CellTranslationType, 
	typename TranslationType::Traits, 
	typename CellTranslationType::Traits>::EXEC_MINUS(cellOp.translation, element.translationPart);
      
      setCanonicalPositionAndNetDirection<Field,Algorithms>(element,cellOp,false);

      element.netDirection  = latticeCoordinate<Field,Algorithms>(element.translationPart);
      element.setGlideFraction();
      return;
    }
  }

  //======================================================================
  /**
   * \brief Return a Symmetry Element representing the given cartesian
   *        symmetry operation.
   *
   */
  template<typename Field, typename LatticeType, typename Algorithms>
  SymmetryElement<Field,2,Algorithms> symmetryElement(const SymmetryOperation<Field,2,Algorithms> cartesianOp,
						      const LatticeType& lat) {
    
    enum {DIM=2};
    typedef SymmetryOperation<Field,DIM,Algorithms> SymmetryOperationType;
    typedef SymmetryElement<Field,DIM,Algorithms>   SymmetryElementType;
    
    static Field                  zero(0);
    SymmetryOperationType         ident;
    SymmetryElementType           result;
    
    result.trace             = cartesianOp.trace();
    
    if (!Algorithms::close(Det(cartesianOp - ident), zero))  {  // has a fixed point
      setRotationElements(result, cartesianOp, lat);
    }
    else { // no fixed point => ident mirror glide or translation
      setMirrorElements(result, cartesianOp, lat);
    }
    result.unNormalizedName = SymmetryElementName<DIM>(result);
    normalize(result);
    result.name = SymmetryElementName<DIM>(result);
    return result;
  }

  //====================================================================== 
  //======================================================================
 
  template<typename Field, typename IN_TYPE, typename Algorithms>
  SymmetryOperation<Field,2,Algorithms> symmetryOperation(IN_TYPE x1, IN_TYPE x2,
							  IN_TYPE x3, IN_TYPE x4 ) {
    SymmetryOperation<Field,2,Algorithms> result;
    result(0,0) = convert<Field>(x1);
    result(0,1) = convert<Field>(x2);
    result(1,0) = convert<Field>(x3);
    result(1,1) = convert<Field>(x4);

    return result;
  }

  template<typename Field, typename IN_TYPE, typename Algorithms>
  SymmetryOperation<Field,2,Algorithms> symmetryOperation(IN_TYPE x1, IN_TYPE x2, IN_TYPE x3,
							  IN_TYPE x4, IN_TYPE x5, IN_TYPE x6 ) {
    SymmetryOperation<Field,2,Algorithms> result;
    result(0,0) = convert<Field>(x1);
    result(0,1) = convert<Field>(x2);
    result(0,2) = convert<Field>(x3);
    result(1,0) = convert<Field>(x4);
    result(1,1) = convert<Field>(x5);
    result(1,2) = convert<Field>(x6);

    return result;
  }

  template<typename Field, typename IN_TYPE, typename Algorithms>
  SymmetryOperation<Field,3,Algorithms> symmetryOperation(IN_TYPE x1, IN_TYPE x2, IN_TYPE x3,
							  IN_TYPE x4, IN_TYPE x5, IN_TYPE x6,
							  IN_TYPE x7, IN_TYPE x8, IN_TYPE x9) {
    SymmetryOperation<Field,3,Algorithms> result;
    result(0,0) = convert<Field>(x1);
    result(0,1) = convert<Field>(x2);
    result(0,2) = convert<Field>(x3);
    result(1,0) = convert<Field>(x4);
    result(1,1) = convert<Field>(x5);
    result(1,2) = convert<Field>(x6);
    result(2,0) = convert<Field>(x7);
    result(2,1) = convert<Field>(x8);
    result(2,2) = convert<Field>(x9);
    return result;
  }

  //====================================================================== 
  /** \ingroup XML
   *
   * SpaceGroup  XML output
   **/
  template<typename Field, size_t DIM, typename LatticeType, typename Algorithms> 
  Tag toXML(const SymmetryOperation<Field,DIM,Algorithms>& symOp,
	    const LatticeType& lat,
	    std::string        name = "SymmetryOperation",
	    std::string        type = "Cartesian" ) {

    typedef SymmetryOperation<double,2,Algorithms>  SymmetryOperationType;
    typedef SeitzPosition<double,2>                 SeitzPositionType;
    typedef CartesianPosition<double,2>             CartesianPositionType;
    typedef SeitzTranslation<double,2>              SeitzTranslationType;
    typedef HermiteNormalForm<Field, DIM+1, DIM+1>  HermiteNormalFormType;
    typedef SeitzMatrix< Field, DIM >               MatType;
    typedef typename MatType::Traits                MatTraitsType;

    Tag result(name);

    result["name"]             = symOp.name;
    result["DIM"]              = DIM;
    result["type"]             = type;
    //result.add(toXML(symmetryElement(symOp, lat)));

    //     Tag hermitNormalFormTag("HermitNormalForm");
    //     Mat<Field,DIM+1,DIM+1> mat;
    //     Copy(mat,b);
    //     result.add(toXML(HermiteNormalFormType(mat)));

    std::ostringstream buff;
    MAT_PRINT<MatType,MatTraitsType>::JUST_NUMBERS( symOp, buff);
    result.content << buff.str();

    return result;
  }

  //   //====================================================================== 
  //   /** \ingroup XML
  //    *
  //    * SpaceGroup  XML output
  //    **/
  //   template<typename Field, size_t DIM, typename Algorithms, typename LatticeType> 
  //   Tag toXML(const SymmetryOperation<Field,DIM,Algorithms>& symOp, 
  // 	    const LatticeType& lat,
  // 	    std::string name="SymmetryOperation") {

  //     Tag result = toXML(symOp,name);

  //     result.add(toXML(lat.cartesianTranslation(symOp.element.netDirection),"CartesianDirection"));
  //     result.add(toXML(lat.cartesianPosition   (symOp.element.cellOffset)  ,"CartesianPosition"));

  //     return result;
  //   }

  //====================================================================== 
  /** \ingroup ostream
   *
   * \brief  SymmetryOperation output stream operator.
   *
   */
  template<typename Field,size_t DIM, typename Algorithms>
  std::ostream& operator << (std::ostream& os, const SymmetryOperation<Field,DIM,Algorithms>& S) {
      
    os << "Symmetry Operation:  " << S.name << "------------------------------" << std::endl;
    os << static_cast<SeitzMatrix< Field, DIM > >(S) << std::endl;
    os << "End Symmetry Operation ------------------------------" << std::endl;
    return os;
  }

} /* namespace psimag */

#endif   //Psimag_SYMMETRY_OPERATION

/*@}*/
