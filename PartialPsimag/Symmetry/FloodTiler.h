//-*-C++-*-

/** \ingroup crystallography */
/*@{*/

/**  \file SuperCrystal.h  
 *
 *  Contains the SuperCrystal class.
 */

#ifndef Psimag_FloodTiler_H
#define Psimag_FloodTiler_H

#include <vector>
#include <ostream>
#include <map>
#include "Vec.h"
#include "PSIMAGAssert.h"
#include "Matrix.h"

#include "Lattice.h"
#include "LatticeWithPattern.h"
#include "CartesianPosition.h"
#include "CartesianTranslation.h"
#include "LatticeCoordinates.h"
#include "Crystal.h"
#include "STLUtil.h"

namespace psimag {


  /** \ingroup crystallography
   *  
   *\brief FloodTiler class .
   *
   */ 
  template<typename Field, size_t DIM, 
	   typename Occupant, 
	   template<typename, size_t, typename> class LatticeTemplate, 
	   typename Algorithms>
  class FloodTiler {
  public:

    typedef Crystal<Field,DIM,Occupant,Algorithms>                            CrystalType;
    typedef CrystalBase<Field,DIM,Occupant,Lattice,Algorithms>                CrystalBaseType;

    typedef LatticeTemplate<Field,DIM,Algorithms>                             LatticeType;
    typedef LatticeWithPattern<Field,DIM,Occupant,LatticeType,Algorithms>     LatticeWithPatternType;
    typedef Pattern<Field,DIM,Occupant,Algorithms>                            PatternType;
    typedef PatternData<Field,DIM,Occupant,Algorithms>                        PatternDataType;
    typedef CellPosition<Field,DIM,Algorithms>                                CellPositionType;
    typedef CartesianPosition<Field,DIM>                                      CartesianPositionType;
    typedef CartesianTranslation<Field,DIM>                                   CartesianTranslationType;

    typedef LatticeCoordinates<DIM>                                           CoordType;
    typedef std::vector<CoordType>                                            CoordsType;

    std::map<CoordType, bool>       inside;
    
    const PatternType&              tilePattern;
    const LatticeType&              tileLattice;
    
    LatticeType                     targetLattice;
    PatternDataType                 tiledPatternData;
    CoordType                       origin;
    
    FloodTiler(const CoordsType&             subLatticeVectors,
	       const LatticeWithPatternType& tileLatticeWithPattern):
      tilePattern(tileLatticeWithPattern.pattern), 
      tileLattice(tileLatticeWithPattern),
      targetLattice(tileLattice.subLattice(subLatticeVectors))
    {
      placeTileAt(origin);
      placeTilesAdjacentTo(origin);
    }
    
    FloodTiler(const LatticeType&            subLat,
	       const LatticeWithPatternType& tileLatticeWithPattern):
      tilePattern(tileLatticeWithPattern.pattern), 
      tileLattice(tileLatticeWithPattern),
      targetLattice(subLat)
    {
      placeTileAt(origin);
      placeTilesAdjacentTo(origin);  // the tile placed at the origin could miss the target lattice
    }
    
    LatticeWithPatternType getTiledLatticeWithPattern() {
      LatticeWithPatternType result(targetLattice, tiledPatternData);
      return result;
    }

    void placeTileAt(const CoordType& coord) {

      static int maxCoord = 1000;
      
      if (inside.count(coord) > 0) {
	return;
      }

      for (size_t i = 0; i< DIM; i++) {
	if(coord[i] >maxCoord) {
	  std::ostringstream message;
	  message << "FloodTiler::placeTileAt: coord[" << i << "] = " << coord[i] << "is too large" << std::endl;
	  throw std::range_error(message.str());
	}
	if(coord[i] < -maxCoord) {
	  std::ostringstream message;
	  message << "FloodTiler::placeTileAt: coord[" << i << "] = " << coord[i] << "is too small" << std::endl;
	  throw std::range_error(message.str());
	}
      }

      CartesianTranslation<Field,DIM> trans = tileLattice.cartesianTranslation(coord);

      bool inTargetLattice = tileInTargetLatticeFor(trans);
      inside[coord]      = inTargetLattice;

      if(inTargetLattice) {
     
	addTileOccupants(trans);
	placeTilesAdjacentTo(coord);
      }

    }

    void placeTilesAdjacentTo(const CoordType& coord) {

      placeTilesAdjacentToInt(coord,DIM-1);
    }
    
    void placeTilesAdjacentToInt(const CoordType& coord, 
				 size_t dim_index) {


      if (dim_index >= 1 ) placeTilesAdjacentToInt(coord,dim_index-1);

      CoordType forwardCoord(coord);
      forwardCoord[dim_index] += 1;
      placeTileAt(forwardCoord);
      if (dim_index >= 1 )  placeTilesAdjacentToInt(forwardCoord,dim_index-1);
      
      CoordType backwardCoord(coord);
      backwardCoord[dim_index] -= 1;
      placeTileAt(backwardCoord);
      if (dim_index >= 1 )  placeTilesAdjacentToInt(backwardCoord,dim_index-1);
      
    }
	
    void addTileOccupants(const CartesianTranslationType& trans) {

      // for each point in the tile pattern
      for(size_t p = 0; p< tilePattern.NUMPOS; p++) {
	  
	// compute its translated cartesian position
	CartesianPositionType newCartPos(trans + tilePattern.cartesianPositions[p]);
	
	// compute its translatedall  position
	CellPositionType newPos = targetLattice.cellPosition(newCartPos);
	
	if (targetLattice.cellContains(newPos)) 
	  tiledPatternData.insertSkipDup(tilePattern.getOccupant(p),newPos);  
      }
    }

    bool tileInTargetLatticeFor(const CartesianTranslationType& trans) {

      CartesianPositionType newTileOrigin = trans.plusOrigin();

      if ( targetLattice.cellContains(newTileOrigin) ) {
	return true;
      }
      for (size_t i = 0; i<DIM; i++) {
	CartesianPositionType corner(tileLattice[i] + newTileOrigin);
	if ( targetLattice.cellContains(corner) ) {
	  return true;
	}
      }
      return false;
    }
    
  };


} /** namespace psimag */

#endif
/*@}*/
