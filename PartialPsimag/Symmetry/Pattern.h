//-*-C++-*-

#ifndef PATTERN_H
#define PATTERN_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file Pattern.h
 *  Contains the class definition for Pattern objects.
 */

#include <iostream>                                  
#include <ostream>
#include <sstream>                                  
#include <iomanip>

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <functional> 

#include "Vec.h"
#include "Tag.h"
#include "Mat.h"
#include "Matrix.h"

#include "PSIMAGAssert.h"
#include "CellPosition.h"
#include "CartesianPosition.h"
#include "CellTranslation.h"
#include "CartesianTranslation.h"

#include "LatticeTransformation.h"
#include "PatternData.h"

namespace psimag {

  /** \ingroup patternConcepts
   *
   * \brief Template class for Patternbase, inherited by PatternWithLattice and TestPattern
   *
   * \note Occupants must be comparable
   *  
   */
  template<typename Field, 
	   size_t   DIM,
	   typename Occupant, 
	   typename Algorithms>
  class Pattern  
  {
    //============================================================ typedefs

  public:

    enum { NROWS=DIM+1};
    
    typedef Pattern<Field,DIM,Occupant,Algorithms>     PatternType;
    typedef PatternData<Field,DIM,Occupant,Algorithms> PatternDataType;
    
    typedef CellPosition<Field,DIM,Algorithms>         CellPositionType;
    typedef CellTranslation<Field,DIM>                 CellTranslationType;
    typedef CartesianTranslation<Field,DIM>            CartesianTranslationType;
    typedef CartesianPosition<Field,DIM>               CartesianPositionType;

    typedef std::vector<CellPositionType>              CellPositionsType;
    typedef std::vector<CartesianPositionType>         CartPositionsType;
    typedef LatticeTransformation<Field,DIM>           LatticeTransformationType;


    //============================================================ members
    
  public:

    /** \brief  The number of positions contained in this pattern. */
    size_t                  NUMPOS;

    /** \brief The list of occupants contained in this
     *         pattern. Indexes into occupants are refered to as
     *         occupantIndicies.
     * 
     * \note Occupants are supposed to be light weight objects.
     *       (Which can however contain pointers or iterators refering
     *       to heavier objects.)
     */    
    std::vector<Occupant>   occupants; 
   
    /** \brief Map the occupantIndex to it's starting posIndex     */
    std::vector<int>        occupantStart;   
    
    /** \brief The number of positions held by the occupant refered to by a occupant index. */
    std::vector<int>        occupantNumPositions;

    /** \brief Map the positionIndex to a Ref to the Occupant The
     *         first DIM Occupants positions refer to the basis
     *         vectors, so the Occupants for these positions are
     *         ignored.
     */
    std::vector<int>        occupant; 

    /**
     * Each column is a cellPosition of an occupant in the pattern,
     * corresponding cellPositions and cartesian positions are at the
     * same index. This same index is used to look up the occupant in
     * occupant, which gives the ocupant index corresponding to the
     * pattern position index.
     */
    CellPositionsType       cellPositions;
    CartPositionsType       cartesianPositions;
    
  public:

    //============================================================ Constructors
    //============================================================

    /** \brief  Construct an empty patternbase Patternbase */
    Pattern(): 
      NUMPOS(0)
    {

    }

   /** \brief  Construct a Patternbase from the given PatternData object. */
    Pattern(const PatternDataType& patternData): 
      NUMPOS(0)
    {
      patternData.loadCellPositionsInto(*this);
    }

    /** \brief  Construct a trivial Patternbase with a single occupant at the origin. */
    Pattern(const Occupant& occupant): 
      NUMPOS(1),
      occupants(1,occupant),
      occupantStart(1,0),
      occupantNumPositions(1,1),
      occupant(1,0),
      cellPositions(1,CellPositionType()),
      cartesianPositions(1,CartesianPositionType())
    {

    }

    /** \brief  Copy Construct a Pattern. */
    Pattern(const Pattern& pat): 
      NUMPOS(pat.NUMPOS),
      occupants(pat.occupants),
      occupantStart(pat.occupantStart),
      occupantNumPositions(pat.occupantNumPositions),
      occupant(pat.occupant),
      cellPositions(pat.cellPositions),
      cartesianPositions(pat.cartesianPositions)
    {

    }
    
    //============================================================ Simple


    /** \brief  Set the number of positions in the pattern. */
    void setNumPos(size_t numpos) {
      NUMPOS = numpos;
      occupant.resize(numpos,-1); // Initialize occupant to -1 which means unused.
      cellPositions.resize(numpos);
      cartesianPositions.resize(numpos);
    }
    
    /** \brief  Set the number of occupants. */
    void setNumOccupants(size_t numOccupants) {
      occupants.resize(numOccupants); // Initialize occupant to -1 which means unused.
      occupantStart.resize(numOccupants);
      occupantNumPositions.resize(numOccupants);
    }

    //============================================================

    /**  Return the ith occupant in this Pattern. */
    const Occupant& getOccupant(int posIndex) const { 
      return occupants[occupant[posIndex]];
    }

    /**  Return the ith occupant in this Pattern. */
    void set(size_t posIndex, size_t occupantIndex, const CellPositionType& cellPos) { 
      occupant[posIndex]      = occupantIndex;
      cellPositions[posIndex] = cellPos;
    }
    
//     CellPositionType cellPosition(size_t posIndex) const {
//       return cellPositions[posIndex];
//     }

//     /**
//      * \brief Return a CartesianPosition object which corresponds to
//      *        the given position index.
//      */
//     CartesianPositionType cartesianPosition(size_t posIndex) const {
//       return cartesianPositions[posIndex];
//     };
  
    /**
     * \brief Return a CartesianTranslation object which corresponds
     *        to the difference between the origin and the cell
     *        position for the given index.
     */
    CartesianTranslationType cartesianTranslation(size_t posIndex) const {
      CartesianTranslationType result;
      for (size_t i=0; i< DIM; i++) {
	result[i] = cartesianPositions[posIndex][i];
      }
      return result;
    };
  
    /**
     * \brief Return a CellTranslation object which corresponds
     *        to the difference between the origin and the cell
     *        position for the given index.
     */
    CellTranslationType cellTranslation(size_t posIndex) const {
      CellTranslationType result;
      for (size_t i=0; i< DIM; i++) {
	result[i] = cellPositions[posIndex][i];
      }
      return result;
    };
  
    /**
     * \brief Returns the position index of the given CellPosition within this
     *        Pattern.
     */
    size_t indexOf(const CellPositionType& position) const {
      
      for(size_t i=0; i < cellPositions.size(); i++) 
	if (cellPositions[i].closeTo(position)) 
	  return i;
      std::ostringstream buff;
      buff << "Pattern.indexOf(" << position << ") could not find an index for the given position!" << std::endl;
      buff << "The Pattern is: " << (*this) << std::endl;
      throw std::out_of_range(buff.str());
    }
    

    //======================================================================

    /** \brief  Fille the given matrix with the cartesian positions. */
    void loadCartesian(Matrix<Field>& mat) {
      mat.resize(NUMPOS,DIM); 
      for (size_t i=0; i< NUMPOS; i++) 
	for (size_t j=0; j< DIM; j++)
	  mat(i,j) = cartesianPositions[i][j];
    }    

    /** \brief  Fille the given matrix with the cell positions. */
    void loadCellPositions(Matrix<Field>& mat) {
      mat.resize(NUMPOS,DIM); 
      for (size_t i=0; i< NUMPOS; i++) 
	for (size_t j=0; j< DIM; j++)
	  mat(i,j) = cellPositions[i][j];
    }    

    //======================================================================
    //======================================================================

    /** 
     * \brief Move each position in this pattern by the given originShift.
     **/
    void shiftPattern(const CartesianTranslationType&  originShift) {
      for(size_t pos=0; pos<NUMPOS; pos++)
	cartesianPositions[pos] = cartesianPositions[pos] + originShift;
    }
    
    //============================================================

    /** 
     * Review and if neccesary change all of the positions so that
     * their coefficients are in the range [0,1).
     *
     **/
    void normalizeCellPositions() {
      for (size_t pos=0; pos < NUMPOS; pos++) 
	cellPositions[pos].normalize();
    }
    
    //============================================================
    
    void buildPlusIndex(Matrix<int>& plusIndex) const {

      plusIndex.resize(NUMPOS, NUMPOS);

      for (size_t p1 = 0; p1 < NUMPOS; p1++) 
	for (size_t  p2 = 0; p2 < NUMPOS; p2++) 
	  plusIndex(p1,p2) = -1;    // This means that the index is undefined, i.e. it would map an occupant into a different occupant

      for (size_t o = 0; o < occupantStart.size(); o++) { 

	size_t oStart = occupantStart[o];
	size_t oEnd   = occupantStart[o] + occupantNumPositions[o];

	for (size_t p1 = oStart; p1 < oEnd; p1++) 
	  for (size_t p2 = oStart; p2 < oEnd; p2++) {

	    CellPositionType newPos;

	    for (size_t i=0; i<DIM; i++) 
	      newPos[i] =  cellPositions[p1][i] + cellPositions[p2][i];
	    
	    newPos.normalize();
	    
	    try {
	      plusIndex(p1,p2) = indexOf(newPos);
	    }
	    catch (std::out_of_range& ) {	  
	      plusIndex(p1,p2) = -2;
	    }
	  }
      }
    }

    //============================================================
    
    void setCartesianSites(Matrix<Field>& sites) const {
      
      sites.resize(cartesianPositions.size(),DIM);
      for (size_t i=0;i<cartesianPositions.size();i++) {
	for (size_t j=0;j<DIM;j++) 
	  sites(i,j)=cartesianPositions[i][j];
      }
    }

    //============================================================
    
    void buildDiffIndex(Matrix<int>& diffIndex) const {
      
      diffIndex.resize(NUMPOS, NUMPOS);

      for (size_t p1 = 0; p1 < NUMPOS; p1++) 
	for (size_t p2 = 0; p2 < NUMPOS; p2++) 
	  diffIndex(p1,p2) = -1;    // This means that the index is undefined, i.e. it would map an occupant into a different occupant
      
      for (size_t o = 0; o < occupantStart.size(); o++) { 
	
	size_t oStart = occupantStart[o];
	size_t oEnd   = occupantStart[o] + occupantNumPositions[o];

	for (size_t p1 = oStart; p1 < oEnd; p1++) 
	  for (size_t p2 = oStart; p2 < oEnd; p2++) {
	    
	    CellPositionType newPos;
	    
	    for (size_t i=0; i<DIM; i++) 
	      newPos[i] =  cellPositions[p1][i] - cellPositions[p2][i];
	    
	    newPos.normalize();
	    try {
	      diffIndex(p1,p2) = indexOf(newPos);
	    }
	    catch (std::out_of_range&) {	  
	      diffIndex(p1,p2) = -2;
	    }
	  }
      }
    }
    //============================================================
  };


  //======================================================================
  
  /** \ingroup XML
   *
   * \brief the Pattern XML outputer
   *
   */
  template<typename Field, typename Occupant, size_t DIM, typename Algorithms>
  Tag addSiteTableXML(const Pattern<Field,DIM,Occupant,Algorithms>& pat,
		      std::string name="SiteTable") {
      
    Tag siteTable(name);
    siteTable["type"] = "Add";
    siteTable["DIM"] = pat.NUMPOS;
    Matrix<int> table;
    pat.buildPlusIndex(table);
    std::ostringstream buff;
    table.print(buff);
    siteTable.content << std::endl << buff.str();
    return siteTable;
  }
    
  //======================================================================
  /** \ingroup XML
   *
   * \brief the Pattern XML outputer
   *
   */
  template<typename Field, typename Occupant, size_t DIM, typename Algorithms>
  Tag diffSiteTableXML(const Pattern<Field,DIM,Occupant,Algorithms>& pat,
		       std::string name="SiteTable") {

    Tag siteTable(name);
    siteTable["type"] = "Diff";
    siteTable["DIM"] = pat.NUMPOS;
    Matrix<int> table;
    pat.buildDiffIndex(table);
    std::ostringstream buff;
    table.print(buff);
    siteTable.content << std::endl << buff.str();
    return siteTable;
  }
    
  //============================================================

  /** \ingroup XML
   *
   * \brief the Pattern XML outputer
   *
   */
  template<typename Field, typename Occupant, size_t DIM, typename Algorithms>
  Tag toXML(const Pattern<Field,DIM,Occupant,Algorithms>& pat,
	    std::string name="Pattern",bool doTables=true) {
    
      Tag result(name);

      //============================================================ Occupants

      Tag occsTag("Occupants");
      for(size_t occ=0; occ<pat.occupants.size(); occ++) {
	Tag occTag("Occupant");
	occTag["patternIndex"] = occ;
	occTag["name"]         = pat.occupants[occ].name;
	occTag["color"]        = pat.occupants[occ].color;
	occTag["startIndex"]   = pat.occupantStart[occ];
	occTag["numPositions"] = pat.occupantNumPositions[occ];
	occsTag.add(occTag);
      }

      result.add(occsTag);

      //============================================================ Occupants

      Tag sitesT("Sites");
      sitesT["numSites"] = pat.NUMPOS;
      for(size_t pos=0; pos<pat.NUMPOS; pos++) {

	Tag siteTag("Site");
	siteTag["index"]=pos;
	siteTag["occupantIndex"]=pat.occupant[pos];
	siteTag["occupantName"] =pat.getOccupant(pos).name;
	siteTag.add(toXML(pat.cellPositions[pos]));
	siteTag.add(toXML(pat.cartesianPositions[pos]));
	sitesT.add(siteTag);
      }
      result.add(sitesT);

      //============================================================ Site Map Tables

      if (doTables) {
	result.add(diffSiteTableXML(pat));
	result.add(addSiteTableXML(pat));
      }
	  
      //============================================================ For Easy Reading

      result.content <<  pat;

      return result;
  }

  //============================================================

  /** \ingroup ostream
   *
   * \brief the Pattern ostream operator
   *
   */
  template<typename Field, typename Occupant, size_t DIM, typename Algorithms>
  std::ostream& operator << (std::ostream& os, 
			     const Pattern<Field,DIM,Occupant,Algorithms>& pat) {

    size_t NUMPOS = pat.NUMPOS;

    os << "Pattern ================================================== Pattern" << std::endl;
    os << "==================================================================" << std::endl;
    os << "| index | atom   | cell position         | cartesian pos          |" << std::endl;
    os << "==================================================================" << std::endl;

    for(size_t pos=0; pos<NUMPOS; pos++) {
      os << "| " <<  std::setw(5) << pos 
	 << " | " <<  std::setw(6) << pat.occupants[pat.occupant[pos]] 
	 << " | " ;
      for (size_t row=0; row<DIM; row++) 
	os <<  " " << std::setw(10) << std::right << pat.cellPositions[pos][row]; 
      os << "|" ;
      for (size_t row=0; row<DIM; row++) 
	os <<  " " << std::setw(10) << std::right << pat.cartesianPositions[pos][row];
      os << " |" ;
      os << std::endl;
    }
    os << "==================================================================" << std::endl;

    return os;
  }
}  

#endif

/*@}*/
