// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides some general type mapping tools
 */

#ifndef MRPAPP_TYPE_MAPPING_H
#define MRPAPP_TYPE_MAPPING_H

#include <type_traits>

namespace mrpapp {
template <typename V1, typename V2, typename T>
struct OnTypesEqual : std::bool_constant<std::is_same<V1, V2>::value> {
  using type = T;
};

template <typename T>
struct default_type : std::true_type {
  using type = T;
};

}  // namespace mrpapp

#endif
