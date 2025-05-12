// Copyright (C) 2025 ETH Zurich
// Copyright (C) 2025 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This file provides the explicit instantiation Vector object

#include "vector.hpp"

namespace mrpapp {

template class Vector<double, DeviceType::CPU>;

}


