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


#include <cassert> //defines assert()
#include <cstring> //defines strcmp()
using namespace std;

#include "apply_trans.h"

namespace minrms
{

void ApplyTransform(ConstMatrix3x4 transform,
                    Biopolymer const& source,
                    Biopolymer &dest)
{
  assert(source.size() == dest.size());
  Biopolymer::const_iterator pSsource = source.begin();
  Biopolymer::iterator pSdest   =   dest.begin();
  while(pSdest != dest.end())
  {
    assert(pSsource != source.end());
    Biopolymer::Residue::const_iterator pAsource = pSsource->begin();
    Biopolymer::Residue::iterator pAdest = pSdest->begin();
    assert(pSsource->size() == pSdest->size());
    while(pAdest != pSdest->end())
    {
      assert(pAsource != pSsource->end());
      assert(pAdest   != pSdest->end());
      assert((*pAsource).second.altLoc == (*pAdest).second.altLoc);
      assert(strcmp((*pAsource).second.name, (*pAdest).second.name) == 0);
      Mult_Mat3x4_by_Vect3(transform,
                           (*pAsource).second.xyz,
                           (*pAdest).second.xyz);
      ++pAsource;
      ++pAdest;
    }
    ++pSsource;
    ++pSdest;
  }
} //ApplyTransform()


void ApplyTransformBackbone(ConstMatrix3x4 transform,
                            Biopolymer const& source,
                            Biopolymer& dest)
{
  assert(source.size() == dest.size());
  Biopolymer::const_iterator pSsource = source.begin();
  Biopolymer::iterator pSdest   =   dest.begin();
  while(pSdest != dest.end())
  {
    assert(pSsource != source.end());
    #ifdef DEBUG
      Real const *pVs;
      Real       *pVd;
    #endif //#ifdef DEBUG

    for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
    {
      Biopolymer::Residue::const_iterator pAsource;
      Biopolymer::Residue::iterator pAdest;
      pAsource = pSsource->GetBackboneAtom(q);
      pAdest   = pSdest->GetBackboneAtom(q);
      assert(pAsource != pSsource->end());
      assert(pAdest   != pSdest->end());
      assert((*pAsource).second.altLoc == (*pAdest).second.altLoc);
      assert(strcmp((*pAsource).second.name, (*pAdest).second.name) == 0);
      #ifdef DEBUG
        pVs = (*pAsource).second.xyz;
        pVd = (*pAdest).second.xyz;
      #endif
      Mult_Mat3x4_by_Vect3(transform,
                           (*pAsource).second.xyz,
                           (*pAdest).second.xyz);
    }
    ++pSsource;
    ++pSdest;
  }
} //ApplyTransformBackbone()


} //namespace minrms
