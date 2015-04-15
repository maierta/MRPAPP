//-*-C++-*-

/** \file NullType.h 
 *  Code in this file is taken from the book: A. Alexandrescu, Modern C++ Design
 *  since the code is similar to that contained in the Loki library, it can 
 *  also be considered to be a modified version of NullType.h in Loki.
 *  Loki comes with the copyright statement in the source code.
 */

///////////////////////////////////////////////////////////////////////////////
// The Loki Library
// Copyright (c) 2001 by Andrei Alexandrescu
// This code accompanies the book:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design 
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
// Permission to use, copy, modify, distribute and sell this software for any 
//     purpose is hereby granted without fee, provided that the above copyright 
//     notice appear in all copies and that both that copyright notice and this 
//     permission notice appear in supporting documentation.
// The author or Addison-Welsey Longman make no representations about the 
//     suitability of this software for any purpose. It is provided "as is" 
//     without express or implied warranty.
///////////////////////////////////////////////////////////////////////////////

#ifndef PSIMAG_NullType_H
#define PSIMAG_NullType_H

namespace psimag {

  /** \brief Used for TypeList termination.
   *  \author Loki
   *  \ingroup Loki
   */
  class NullType {};

}

#endif // PSIMAG_NULLType_H
