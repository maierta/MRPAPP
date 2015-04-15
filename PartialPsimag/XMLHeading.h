//-*-C++-*-

#ifndef XML_HEADING_H
#define XML_HEADING_H

/** \ingroup ostream **/
/*@{*/

/** \file Tag.h
 *  Contains the class definition for XML Tag and related objects.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <ostream>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <functional>

#include "PSIMAGAssert.h"

namespace psimag {

  class XMLHeading {
  public:
    typedef enum { XHTML, SVG} DocType;
    std::string version;
    std::string standAlone;
    DocType     docType;

    XMLHeading(DocType docType_=XMLHeading::XHTML):
      version("1.0"),
      standAlone("no"),
      docType(docType_)
    {}
    
    const std::string docTypeString() const {
      switch (docType) { 
      case SVG:   return std::string("svg  PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\" ");
      case XHTML: return std::string("html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"DTD/xhtml1-strict.dtd\" ");
      default:    return std::string("html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"DTD/xhtml1-transitional.dtd\" ");
      }
    }
  };

/* ATTENTION/FIXME: Header files should not contain implementations
        (because it causes multiple definitions when using more than one compilation unit)
        quick fix: I added the keyword inline (Gonzalo) */

  /** \ingroup ostream
   *   XML Heading ostream operator
   **/
  inline std::ostream& operator << (std::ostream& os, const XMLHeading& h ) {
    os << "<?xml version=\"" << h.version << "\" standalone=\""<< h.standAlone << "\"?>"  << std::endl;
    os << "<!DOCTYPE " << h.docTypeString() << ">" << std::endl;
    os << std::endl;
    return os;
  }

}  

#endif

/*@}*/
