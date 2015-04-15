//-*-C++-*-

#ifndef PSIMAG_LatticeTransformationBase_H
#define PSIMAG_LatticeTransformationBase_H

/** \ingroup latticeConcepts */
/*@{*/

/** \file Since templates and virtuals don't mix a fragment of code to acheive the same effect. */

    //====================================================================== Common base code for LatticeTransformation and 
// ======================================================================    Inverse Lattice Transforamtion
    //======================================================================

    //---------------------------------------------------------------------- LatticeTemplate
    /** 
     * \brief Return a transformed version of the given lattice.
     *
     */
    template<typename Algorithms, 
	     template<typename,size_t,typename> class InLatticeTemplate,
	     template<typename,size_t,typename> class OutLatticeTemplate>
    OutLatticeTemplate<Field,DIM,Algorithms>
    operator() (const InLatticeTemplate<Field,DIM,Algorithms>& srcLattice) const 
    {
      OutLatticeTemplate<Field,DIM,Algorithms> result;
      (*this)(srcLattice, result);
      return result;
    }

    //---------------------------------------------------------------------- SymmetryOperationType
    /** 
     * \brief Return a transformed version of the given SymmetryOperation.
     *
     */
    template<typename Algorithms>
    SymmetryOperation<Field,DIM,Algorithms>  
    operator()( const SymmetryOperation<Field,DIM,Algorithms>& symop) const {   
      SymmetryOperation<Field,DIM,Algorithms> result;
      (*this)(symop, result);
      return result;
    }

    //---------------------------------------------------------------------- PatternType
    /** 
     * \brief Transform the given PatternType.
     *
     */
    template<typename Occupant, 
	     typename Algorithms>
    void operator()
    (const Pattern<Field, DIM, Occupant, Algorithms>& srcPattern,
     Pattern<Field, DIM, Occupant, Algorithms>& tgtPattern) const {
      
      tgtPattern.setNumPos(srcPattern.NUMPOS);

      tgtPattern.occupants            = srcPattern.occupants;
      tgtPattern.occupant             = srcPattern.occupant;
      tgtPattern.occupantStart        = srcPattern.occupantStart;
      tgtPattern.occupantNumPositions = srcPattern.occupantNumPositions;
      
      // Transform each cell position
      for (size_t i=0; i<srcPattern.NUMPOS; i++)
	tgtPattern.cellPositions[i] = (*this)(srcPattern.cellPositions[i]);
      
      tgtPattern.normalizeCellPositions();
    }
    
    //---------------------------------------------------------------------- LatticeWithPatternType
    /** 
     * \brief Transform the given LatticeWithPatternType.
     *
     */
    template<typename Occupant, 
	     typename SrcLatticeType,
	     typename TgtLatticeType,
	     typename Algorithms>
    void operator()
      (const LatticeWithPattern<Field, DIM, Occupant, SrcLatticeType, Algorithms>& srcLatticeWithPattern,
       LatticeWithPattern<Field, DIM, Occupant, TgtLatticeType, Algorithms>& tgtLatticeWithPattern) const {

      // Transform Lattice Part
      const SrcLatticeType&  srcLattice = srcLatticeWithPattern;
      TgtLatticeType&        tgtLattice = tgtLatticeWithPattern;
      (*this)(srcLattice, tgtLattice);

      typedef Pattern<Field, DIM, Occupant, Algorithms>  PatternType;

      // Transform Pattern Part
      const PatternType& srcPattern = srcLatticeWithPattern.pattern;
      PatternType& tgtPattern = tgtLatticeWithPattern.pattern;
      (*this)(srcPattern, tgtPattern);

      tgtLatticeWithPattern.generateCartesianPositions();
    }
    
    //---------------------------------------------------------------------- SymmetryElementType
    /** 
     * \brief Transform the given symmmetry element into the target
     *        symmetry element.
     */
    template<typename Algorithms>
    void operator()(const SymmetryElement<Field,DIM,Algorithms>& el,  
		    SymmetryElement<Field,DIM,Algorithms>& tgtEl) const {   

      tgtEl.type              = el.type;
      tgtEl.name              = el.name;
      tgtEl.unNormalizedName  = el.unNormalizedName;
      tgtEl.netDirection      = (*this)(el.netDirection);
      tgtEl.translationPart   = (*this)(el.translationPart);
      tgtEl.glideFraction     = el.glideFraction;
      tgtEl.cellPosition      = (*this)(el.cellPosition);
      tgtEl.id                = el.id;
      tgtEl.trace             = el.trace;
      tgtEl.rotationAngle     = el.rotationAngle;
      tgtEl.sense             = el.sense;
    }


    //---------------------------------------------------------------------- AppliedSymmetryElementType
    /** 
     * \brief Transform the given AppliedSymmetryElementType.
     *
     */
    template<typename Occupant, 
	     typename SrcLatticeType,
	     typename TgtLatticeType,
	     typename Algorithms>
    void operator()
    (const AppliedSymmetryElement<Field,DIM,Occupant,SrcLatticeType,Algorithms>& srcAppliedSymmetryElement,
     AppliedSymmetryElement<Field,DIM,Occupant,TgtLatticeType,Algorithms>& tgtAppliedSymmetryElement) const {

      typedef SymmetryElement<Field, DIM, Algorithms> SymmetryElementType;

      // Transform SymmetryElement part
      (*this)(srcAppliedSymmetryElement.element, tgtAppliedSymmetryElement.element);

      // Cartesian parts don't need transformation (optionally doublechecking them would be nice!) &*&*&*
      tgtAppliedSymmetryElement.cartesianTranslation = srcAppliedSymmetryElement.cartesianTranslation;
      tgtAppliedSymmetryElement.cartesianPosition    = srcAppliedSymmetryElement.cartesianPosition;

      // These scalar quantites don't need to be transformed (optionally doublechecking them would be nice!) &*&*&*
      tgtAppliedSymmetryElement.latticeDistance =  srcAppliedSymmetryElement.latticeDistance;
      tgtAppliedSymmetryElement.patternDistance =  srcAppliedSymmetryElement.patternDistance;
      tgtAppliedSymmetryElement.distance        =  srcAppliedSymmetryElement.distance;
      tgtAppliedSymmetryElement.id              =  srcAppliedSymmetryElement.id;

      // Transform the associated Symmetryoperation
      (*this)(srcAppliedSymmetryElement.operation,  tgtAppliedSymmetryElement.operation);

      // Transform the associated testPattern (transform the cell positions and regenerate cartesian positions)
      (*this)(srcAppliedSymmetryElement.testPattern, tgtAppliedSymmetryElement.testPattern);
    }
    

#endif //PSIMAG_LatticeTransformationBase_H

/*@}*/

