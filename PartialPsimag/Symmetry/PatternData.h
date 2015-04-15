//-*-C++-*-

#ifndef TEMP_PATTERN_H
#define TEMP_PATTERN_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file Pattern.h
 *  Contains the class definition for Temporary Pattern objects.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <ostream>

#include "PSIMAGAssert.h"
#include "CellPosition.h"
#include "Pattern.h"

namespace psimag {

  template<typename Field, size_t DIM, typename Occupant, typename Algorithms>
  class Pattern;

  /** \ingroup patternConcepts
   *
   * \brief Template class for Pattern.
   *
   * \note Occupants must be comparable
   *  
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms>
  class PatternData  {
    
  public:
    
    typedef Pattern<Field,DIM,Occupant,Algorithms>       PatternType;
    typedef CellPosition<Field,DIM, Algorithms>          CellPositionType;
    
    typedef std::set<CellPositionType>                   CellPositions;
    typedef std::set<Occupant>                           Occupants;
    typedef std::map<Occupant,CellPositions>             OccupantMap;
    
    typedef typename Occupants::const_iterator           OccItr;
    typedef typename OccupantMap::const_iterator         OccMapItr;
    typedef typename CellPositions::const_iterator       CellPosItr;

  public: 
    
    /** A map of occupants to the set of cell positions where they are located. */
    OccupantMap      occupantPositions;
    
    /** The set of positions that are contained in the this pattern. */
    CellPositions    positions;

    /** The set of Occupants that are contained in this pattern. */
    Occupants        occupants;
    
  public:
    
    /** \brief  Construct an empty pattern PatternData */
    PatternData() {}
    
    /** \brief Insert the given cellPosition into the CellPositions
     *         associated with the given occupant.
    */
    void insertRaw(const Occupant& occupant, const CellPositionType& cellPosition) {

      bool insertSucceded = positions.insert(cellPosition).second;
      if (!insertSucceded) {
	std::ostringstream msgBuffer;
	msgBuffer << "In Pattern.positions: unsuccesfull attempt to add position:" << std::endl;
	msgBuffer << "    cellPosition: " << cellPosition << std::endl;
	throw std::length_error(msgBuffer.str());
      }
      
      insertSucceded = occupantPositions[occupant].insert(cellPosition).second;
      if (!insertSucceded) {
	std::ostringstream msg;
	msg << "In Pattern[oc]=pos unsuccesfull attempt to add occupant and position. \n";
	throw std::length_error(msg.str());
      }
	
      occupants.insert(occupant);
    }

    /** \brief Insert the given cellPosition into the CellPositions
     *         associated with the given occupant.
    */
    void insert(const Occupant& occupant, const CellPositionType& cellPosition) {

      CellPositionType normalized = cellPosition;
      normalized.normalize();
      insertRaw(occupant,normalized);
    }

    /** \brief Insert the given cellPosition into the CellPositions
     *         associated with the given occupant.
    */
    void insertSkipDup(const Occupant& occupant, const CellPositionType& cellPosition) {

      CellPositionType normalized = cellPosition;
      normalized.normalize();

      if (positions.count(normalized) > 0) {
	return;
      }
      insertRaw(occupant, normalized);
    }

//     /** \brief Insert the given cellPosition into the CellPositions
//      *         associated with the given occupant.
//     */
//     void insert(Occupant occupant, CellPositionType cellPosition) {
//       insert(occupant,cellPosition);
//     }

    //============================================================

    /**
     * Use this patternData to setup the given pattern.
     *
     */
    void loadCellPositionsInto (PatternType& pat)  const {

      size_t     posIndex = pat.cellPositions.size();
      size_t     occIndex = pat.occupant.size();
      OccItr     o;
      CellPosItr ocp;
      
      for(o  = occupants.begin(); 
	  o != occupants.end(); 
	  o++) {

	const Occupant& occupant = *o;
	pat.occupants.push_back(occupant);
	pat.occupantStart.push_back(posIndex);
	const CellPositions&  dataPositions   = occupantPositions.find(occupant)->second; 

	size_t         numPos           = dataPositions.size();
	pat.occupantNumPositions.push_back(numPos);
	
	for(ocp=dataPositions.begin(); ocp != dataPositions.end(); ocp++) {
	  pat.cellPositions.push_back(*ocp);
	  pat.occupant.push_back(occIndex);
	  posIndex++;
	}
	occIndex++;
      }
      pat.NUMPOS = pat.cellPositions.size();
      pat.cartesianPositions.resize(pat.NUMPOS);

    }

    //============================================================

    /** 
     *  Return the Occupant object for the given position.
     */
    const Occupant& getOccupant(CellPositionType& position) const { 
      
      OccMapItr occPositions;
      for(occPositions = occupantPositions.begin();
	  occPositions != occupantPositions.end(); 
	  occPositions++) 
	if (occPositions.second.count(position) == 1)
	  return occPositions.first;
      
      throw std::range_error("PatternData.getOccupant(position) not found");
    }
  };


  //============================================================

  /** \ingroup ostream
   *
   * \brief the Pattern ostream operator
   *
   */
  template<typename Field, size_t DIM, typename Occupant, typename Algorithms>
  std::ostream& operator << (std::ostream& os, const PatternData<Field,DIM,Occupant, Algorithms>& patData) {

    typedef PatternData<Field,DIM,Occupant,Algorithms>        PatternDataType;

    typedef typename PatternDataType::CellPositionType        CellPositionType;
    typedef typename PatternDataType::OccupantMap             OccupantMap;
    typedef typename PatternDataType::CellPositions           CellPositions;
    typedef typename PatternDataType::OccItr                  OccItr;
    typedef typename PatternDataType::OccMapItr               OccMapItr;
    typedef typename PatternDataType::CellPosItr              CellPosItr;

    //os.setf(std::ios_base::fixed, std::ios_base::floatfield); os.precision(6);

    os << " Positions -----------------------" << std::endl;

    CellPosItr cellPos;
    for(cellPos= patData.positions.begin();
	cellPos != patData.positions.end();
	cellPos++) 
      os << "      " << (*cellPos) << std::endl;

    os << " Occupants -----------------------" << std::endl;

    OccItr occ;
    for(occ = patData.occupants.begin();
	occ != patData.occupants.end();
	occ++) 
      os << "      " << (*occ) << std::endl;


    os << " Occupant to Positions -----------------------" << std::endl;

    OccMapItr occPositions;
    for(occPositions = patData.occupantPositions.begin();
	occPositions != patData.occupantPositions.end(); 
	occPositions++) {

      os << "  " << occPositions->first 
	 << "  " << occPositions->second.size() << " positions: " << std::endl;

      const CellPositions positions = occPositions->second;

      if (positions.size() == 1) {
	os << " { " << (*(positions.begin())) <<  "}" << std::endl;
	continue;
      }

      os << std::endl << " { ";
      bool first = true;
      for(cellPos = positions.begin();
	  cellPos != positions.end(); 
	  cellPos++) {
	
	if (first) 
	  first = false;
	else
	  os << "   ";
	os << (*cellPos) << std::endl;
      }
      os << " }" << std::endl;
    }
    return os;
  }
  
}  

#endif

/*@}*/
