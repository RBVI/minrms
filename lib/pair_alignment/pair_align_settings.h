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

#ifndef _PAIR_ALIGN_SETTINGS_H
#define _PAIR_ALIGN_SETTINGS_H


#include <vector>
#include <string>
using namespace std;

#include <biopolymer.h>

//Settings that will be useful for most algorithms that aligns
//two structures together (using intermolecular criteria.)
namespace minrms
{


struct PairAlignSettings
{
  vector<string> vPdbFileNames; //Names of the PDB files containing
                                //the structures being aligned

  Biopolymer const *aMol; //An array containing the two structures being
                          //aligned (exactly as read from the PDB file)

  Biopolymer const *aMol_f; //A subset of the structures being aligned
                            //storing only the subset of atoms/residues
                            //that will be used in the computation.
                            //Note: A separate class is responsible for
                            //      allocating and filling aMol and aMol_f
                            //      (called PairAlignSettingsParser).
  int n_min; //The minimum number of matches made in an alignment
             //
             //(Details:   (can skip)
             // This information is useful for both the alignment generator
             // (using the lower-bound optimization)
             // and for the orientation generator,
             // for determining how to orient the
             // structures when sumperimposing them together.
             // Since n_min is shared by several different classes,
             // I put it here, where they can all access it.)

  vector<string> vMsfOutputLabels;//The labels that will be used to identify
                                     //which structure corresponds to which
                                     //line of each MSF-file generated
                                     //by the program.

  bool msf_use_lower_case; //Use Conrad's compression trick?

  PairAlignSettings();

#ifdef COMPILER_BUG2
  PairAlignSettings(PairAlignSettings const &s);
#endif //#ifdef COMPILER_BUG2

}; //struct PairAlignSettings


} //namespace minrms

#endif //#ifndef _PAIR_ALIGN_SETTINGS_H

