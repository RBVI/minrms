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



#ifndef _APPLY_TRANS_H
#define _APPLY_TRANS_H
#include <vect_3d.h>
#include "biopolymer.h"

namespace minrms {

//Multiplies all of the atoms in molecule "source" by a 3x4 matrix "transform"
//(this matrix stores an affine-transformation).
//The resulting positions are stored in the "dest" molecule,
//The "dest" molecule must be pre-allocated and
//(apart from the position of the atoms) identical
//to the "source" molecule.
//Note:  It's okay if source and dest are the same.
void ApplyTransform(ConstMatrix3x4 transform,
                    Biopolymer const& source,
                    Biopolymer &dest);

//Same as above, only it saves some time, perhaps by only transforming the
//"backbone-atoms", and ignoring all the other atoms.
//If all the calculations involve "backbone-atoms" exclusively, then you should
//use this version instead of "ApplyTransform()"
//Note:  Again, it's okay if source and dest are the same.
void ApplyTransformBackbone(ConstMatrix3x4 transform,
                            Biopolymer const& source,
                            Biopolymer& dest);

} //namespace minrms

#endif //#ifndef _APPLY_TRANS_H
