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

#include "search_or.h"

using namespace minrms;

//const int   SearchOr::ITERATE_UNTIL_CONVERGENCE = -1;
//const int   SearchOr::DISSABLE_BACKUP = -1;
const Real SearchOr::Settings::ALTERNATES_NOT_LIMITED_BY_RMSD = -1.0f;


#ifdef COMPILER_BUG2
SearchOr::Settings::Settings()
 :PairAlignSettings(),
  num_solutions(1),
  alt_method(ALT_METHOD_3D),
  alt_min_3d_difference(0.0),
  alt_max_pairs_common(0),
  alt_rmsd_tolerance(ALTERNATES_NOT_LIMITED_BY_RMSD),
  refine_method(REFINE_BEST),
  refine_max_iters(ITERATE_UNTIL_CONVERGENCE),
  refine_while_searching(true)
{
  //cerr << "SearchOr::Settings, this = " << this << endl;
  //cerr << "default constructor called." << endl;
}

    //copy constructor:
    //(The alpha compiler produces buggy
    // binaries unless I define one explicitely)

SearchOr::Settings::Settings(Settings const& s):
  PairAlignSettings(s),
  num_solutions(s.num_solutions),
  alt_method(s.alt_method),
  alt_min_3d_difference(s.alt_min_3d_difference),
  alt_max_pairs_common(s.alt_max_pairs_common),
  alt_rmsd_tolerance(s.alt_rmsd_tolerance),
  refine_method(s.refine_method),
  refine_max_iters(s.refine_max_iters),
  refine_while_searching(s.refine_while_searching)
{
  // ----- hack ------
  //For some reason, the parent class's copy constructor
  //        PairAlignSettings(Settings const& s)
  //is not geting called.  Instead, its default constructor is invoked.
  //Don't know what's wrong, but I got to get code working sometime.
  //Rather than rely on the parent's copy constructor,
  //I just copy what I need here.
  //This code is presently only compiled on the alpha anyway.
  //(COMPILE_BUG2 is only defined when compiling on the alpha.)

  //These are all members of the parent class:
  vPdbFileNames = s.vPdbFileNames;
  aMol = s.aMol;    //copy by (const) reference, not by value
  aMol_f = s.aMol_f;
  n_min = s.n_min;
  vMsfOutputLabels = s.vMsfOutputLabels;
  msf_use_lower_case = s.msf_use_lower_case;
  // -----------------

  //cerr << "SearchOr::Settings, this = " << this << endl;
  //cerr << "s.aMol = " << s.aMol << ";  aMol = " << aMol << endl;
  //cerr << "s.aMol_f = " << s.aMol_f << ";  aMol_f = " << aMol_f << endl;
}
#endif //#ifdef COMPILER_BUG2


