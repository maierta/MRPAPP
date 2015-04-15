//-*-C++-*-

#ifndef  Psimag_GroupAction
#define  Psimag_GroupAction

/** \ingroup symmetryConcepts */
/*@{*/

/** \file GroupAction.h
 *
 */  

#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>

#include "SymmetryOperation.h"
#include "AppliedSymmetryElement.h"
#include "Pattern.h"
#include "LatticeWithPattern.h"
#include "TestPattern.h"
#include "CellPosition.h"
#include "CartesianPosition.h"
#include "SymmetryElement.h"
#include "AppliedSymmetryElement.h"
#include "Matrix.h"
#include "SymmetryGroup.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class SymmetryGroup;

  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   typename LatticeType,
	   typename Algorithms>
  class AppliedSymmetryElement;

  /** \ingroup symmetryConcepts 
   *
   * \brief The Group class
   *
   * \ note http://en.wikipedia.org/wiki/Group_action
   */
  template<typename Field, size_t DIM, typename Occupant,
	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  class GroupAction: 
    public Matrix<int>
  {
  public:

    typedef LatticeTemplate<Field,DIM,Algorithms>                             LatticeType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                            PatternType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>     LatticeWithPatternType;

    typedef CartesianPosition<Field,DIM>                                      CartesianPositionType;
    typedef CellPosition<Field,DIM, Algorithms>                               CellPositionType;

    typedef SymmetryGroup<Field,DIM,Occupant,LatticeTemplate,Algorithms>      ThisType;

    typedef AppliedSymmetryElement<Field,DIM,Occupant,LatticeType,Algorithms> AppliedSymmetryElementType;
    typedef std::vector<AppliedSymmetryElementType*>                          AppliedSymmetryElementsType;
    typedef std::pair<int,AppliedSymmetryElementType>                         AppliedSymmetryElementPairType;
    typedef std::map<std::string,AppliedSymmetryElementPairType>              GroupType;
    typedef SymmetryElement<Field,DIM,Algorithms>                             SymmetryElementType;
    typedef SymmetryOperation<Field,DIM,Algorithms>                           SymmetryOperationType;
    typedef TestPattern<Field,DIM,Occupant,LatticeType,Algorithms>            TestPatternType;
    typedef const std::pair<int,int>                                          IntPairType;
    typedef std::map<IntPairType,int>                                         MultiplicationTableMapType;
    typedef std::vector<int>                                                  PointPermutationType;


    size_t numPatternPos;
    size_t numOperations;

    const AppliedSymmetryElementsType&       appliedElements; // an index into group
    const PatternType&                       pattern; 

    /** \ingroup symmetryConcepts 
     *
     * \brief The SymmetryGroup class
     */
    GroupAction(const AppliedSymmetryElementsType&  applied,
		const PatternType&                  pat):
      numPatternPos(pat.NUMPOS),
      numOperations(0),
      appliedElements(applied),
      pattern(pat)
    {}

    /** \ingroup symmetryConcepts 
     *
     * \brief The SymmetryGroup class
     */
    void init()
    {
      numOperations = appliedElements.size();
      
      (*this).resize(numOperations, numPatternPos);
      
      // Initialize the group action table
      for (size_t o = 0; o < numOperations; o++) 
 	for (size_t p = 0; p < numPatternPos; p++) 
	  (*this)(o,p) = -1;    // This means that the index is undefined, i.e. it would map an occupant into a different occupant
      
      // Fill the table up
      for (size_t o = 0; o < numOperations; o++) {
	
 	const AppliedSymmetryElementType& op           = *appliedElements[o];
	const TestPatternType&            testPattern  = op.testPattern;

	for (size_t p = 0; p < numPatternPos; p++) {
	  try {
	    (*this)(o,p) = testPattern.pattern.indexOf(pattern.cellPositions[p]);
	  }
	  catch (std::out_of_range& e) {	
	    const AppliedSymmetryElementType& op           = *appliedElements[o];
	    //const TestPatternType&            testPattern  = op.testPattern;
	    std::ostringstream msg;
	    msg << "====================================================================== Group Action Error" << std::endl;
	    msg << "pattern is " << std::endl;
	    msg << " " << pattern << std::endl;
	    msg << "for operation[" << o << "]" << op.element.name << "(" << p << "=" << pattern.cellPositions[p] << ")" << std::endl;
	    msg << "   " << e.what() << std::endl;
	    msg << "operation is" << std::endl;
	    msg << " " << toXML(op) << std::endl;
	    //msg << toXML(testPattern) << std::endl;
	    throw std::logic_error(msg.str());
	    (*this)(o,p) = -2;
	  }
	}
      }
    }
    
    void print(std::ostream& os, std::vector<int> subsetOfOps=std::vector<int>(0)) const {

      bool   doAll  = (subsetOfOps.size()==0);
      size_t numOps = doAll?this->n_row():subsetOfOps.size();

      os << std::endl;
      for (size_t i = 0; i < numOps; i++) {
	int r = doAll?i:subsetOfOps[i];
	AppliedSymmetryElementType& appliedElement = *appliedElements[r];
	for (size_t c = 0; c < numPatternPos; c++) 
	  os << std::setw(3) << right << (*this)(r,c) << " ";
	os << "          /* " 
	   << std::setw(3) << left << r << "    "
	   << std::setw(30) << left << appliedElement.element.name << "*/";
	os << std::endl;
      }
    }
    
    bool hasSameAction(int op1, int op2) const {

      size_t numPos = (*this).n_col();
      for(size_t i=0; i<numPos;i++) {
	if ((*this)(op1,i) != (*this)(op2,i) )
	  return false;
      }
      return true;
    }

    bool isRedundantAction(int op, std::vector<int>& ops) const {
      
      // We already have collected op in ops
      if (count(ops.begin(),ops.end(),op) > 0) return true;
      
      for (size_t i=0; i< ops.size(); i++)
	if(hasSameAction(op,ops[i]))
	  return true;
      
      return false;
    }

    std::vector<int> uniquePermutationOps() const {

      size_t numOps         = (*this).n_row();
      std::vector<int>                  result(0);

      // Iterate over all ops
      for (size_t op=0; op < numOps; op++) 
	if(!isRedundantAction(op,result))
	  result.push_back(op);
      
      return result;
    }

    //============================================================ 
    /** \ingroup symmetryConcepts 
     *
     * \brief Return a matrix containing a unique set of 0group action permutations
     *
     * \ note buildGroupAction(const SymmetryType& symmetry) must be called first
     */
    void uniqueActions(Matrix<int>& permutations) const {
      
      std::vector<int> unique = uniquePermutationOps();
      size_t       numPositions  = (*this).n_col();
      size_t       numOperations = unique.size();
      
      permutations.resize(numOperations, numPositions);

      for (size_t i=0; i< numOperations; i++) 
	for(size_t p=0; p< numPositions; p++)
	  permutations(i,p) = (*this)(unique[i],p);
    }

  };

  template<typename Field, size_t DIM, typename Occupant,
 	   template<typename, size_t, typename> class LatticeTemplate,
	   typename Algorithms>
  Tag toXML(const GroupAction<Field,DIM,Occupant,LatticeTemplate,Algorithms>& groupAction,
	    std::string name="GroupAction",
	    const std::vector<int> selectedOPs=std::vector<int>(0) ) 
  {
    
    int numSelected = selectedOPs.size();
    Tag result(name);
    result["rows"] = (numSelected == 0)?groupAction.n_row():numSelected;
    result["cols"] = groupAction.n_col();
      
    std::ostringstream buff;
    groupAction.print(buff,selectedOPs);
    result.content << std::endl << std::endl <<buff.str();
    
    return result;
  }
  
} /** namespace spimag **/

#include "Symmetry.h"

#endif

/*@}*/
