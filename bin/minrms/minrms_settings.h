//
// Copyright (c) 2002 The Regents of the University of California.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions, and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions, and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//   3. Redistributions must acknowledge that this software was
//      originally developed by the UCSF Computer Graphics Laboratory
//      under support by the NIH National Center for Research Resources,
//      grant P41-RR01081.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _MINRMS_SETTINGS_H
#define _MINRMS_SETTINGS_H

//
// This file contains declarations for the "MinrmsSettings" class.
// MinrmsSettings contains all of the initial settings and
// input data that will be used throughout the entire program.
//
//      Minrms is divided into several classes:
// an orientation generator, an alignment generator
// (which is a dynamic programming algorithm),
// and an orientation searcher ("SearchOr").
// (-The orientation searcher is the glue between
//   the first two classes.)
//
//      Similarly, the code that handles storing and reading the
// settings necessary for these classes to operate, is organized
// into separate classes as well.  Because some of the data
// necessary for an OrientationGenerator class, for example,
// is also the same data required for a "SearchOr" class,
// I organized this data in a hierarchy.
//      The "MinrmsSettings" class multiply-inherits from
// the all of the classes which handle the settings
// for all of the complenents that make up minrms.
//
//
// MinrmsSettings class Hierarchy.
//
//            PairAlignSettings
//              /          `-._ 
//             /               `-._
//            /                    `-._
//           /                 SearchOrSettings
//          /                 /                \
// OrientationGenerator  SearchOrNW::    SearchOrDyn3d::
//    ::Settings         Settings          Settings
//           \           |                  _.-'
//            \          |              _.-'
//             \         |          _.-'
//              \        |      _.-'
//               \       |  _.-'
//               MinrmsSettings
//
//  (Also see "minrms_parser.h".)
//

#include "search_or.h"
#include "search_or_dyn3d.h"
#include "search_or_nw.h"

namespace minrms
{


struct MinrmsSettings: public virtual OrientationGenerator::Settings,
                       public virtual SearchOrNW::Settings,
                       public virtual SearchOrDyn3d::Settings
{
  enum WhichAlgorithm{ALGO_DYN3D, ALGO_NEEDLEMAN_WUNSCH};
  WhichAlgorithm which_algorithm;
}; //struct //SearchOrSettings



} //namespace minrms

#endif //#indef _MINRMS_SETTINGS_H

