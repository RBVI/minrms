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

#include "pair_align_settings.h"


namespace minrms
{
PairAlignSettings::PairAlignSettings():
  vPdbFileNames(2),
  aMol(NULL),
  aMol_f(NULL),
  n_min(0),
  vMsfOutputLabels(2),
  msf_use_lower_case(false)
{
  //cerr << "PairAlignSettings(), this = " << this << endl;
  //cerr << "...default constructor invoked" << endl;
}
  
#ifdef COMPILER_BUG2
// Copy constructor
// The alpha compiler produced buggy code
// unless I explicitely define one.
PairAlignSettings::PairAlignSettings(PairAlignSettings const &s):
  vPdbFileNames(s.vPdbFileNames),
  n_min(s.n_min),
  vMsfOutputLabels(s.vMsfOutputLabels),
  msf_use_lower_case(s.msf_use_lower_case)
{
  aMol = s.aMol;    //copy by (const) reference, not by value
  aMol_f = s.aMol_f;
  //cerr << "PairAlignSettings(), this = " << this << endl;
  //cerr << "s.aMol = " << s.aMol << ";  aMol = " << aMol << endl;
  //cerr << "s.aMol_f = " << s.aMol_f << ";  aMol_f = " << aMol_f << endl;
}
#endif //#ifdef COMPILER_BUG2

} //namespace minrms

