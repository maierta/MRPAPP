//-*-C++-*-

#ifndef PSIMAG_Occupant_H
#define PSIMAG_Occupant_H

/** \ingroup patternConcepts **/
/*@{*/

/** \file Occupant.h
 *  Contains the class definition for generic Occupant objects.
 */

#include <string>
#include <vector>

namespace psimag {

  /** \ingroup patternConcepts
   *
   * \brief Class definition for  Simple Occupant. Usefull for testing.
   */
  class Occupant {
    
  public:
    
    Occupant(): 
      name("None"), 
      color("red") 
    {};
    
    Occupant(char* cString): 
      name(cString), 
      color("red") 
    {};

    Occupant(const Occupant& aOccupant): 
      name(aOccupant.name), 
      color(aOccupant.color) 
    {};

    Occupant(std::string aName): 
      name(aName), 
      color("red") 
    {};

    Occupant(std::string aName,std::string color_): 
      name(aName), 
      color(color_)
    {};
    
    bool operator==(const Occupant &occupant) const {
      if(name==occupant.name) return true;
      return false;
    }

    bool operator< (const Occupant& other) const {
      return (this->name < other.name);
    }
 
    Occupant& operator=(const Occupant& occupant) {
      name  = occupant.name;
      color = occupant.color;
      return *this;
    }
  
    std::string name;
    std::string id;
    std::string color;
  
  };

  /** \ingroup ostream
   *
   * \brief The Occupant ostream operator.
   *
   */
inline  std::ostream& operator << (std::ostream& os, const Occupant& occupant) {
    os << occupant.name;
    return os;
  }

}

#endif

/*@}*/
