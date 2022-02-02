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

#ifndef _TWO_DIM_TABLE_H
#define _TWO_DIM_TABLE_H

#include <global_utils.h>

namespace minrms {

// The TwoDimTable class is actually a three-dimensional array
//  of width = size_x, height = size_y, and depth = 2.
// It is meant to be used as a part of the Dyn3d class.
//
// What it does:
// In order to reduce the size of the dynamic programming table
// stored in the Dyn3d class as much as possible, only the pointers
// which point to the previous cells are needed later on,
// so only they must be stored in a fully three dimensional table.
// Pointers in the three dimensional table
// can be represented by as few as two bits.
//   (however, we chose to use one
//    byte for each pointer instead).
// The other data from that table (like the total
// "cost", a.k.a. Sum-Squared-Distance) can be stored in
// a much smaller table, like the TwoDimTable here,
// because only the current and the most recent layer of
// the dynamic programming table are accessed during
// the the calculation of the alignments.

template< class T >
class TwoDimTable {
  TwoDimTable(const TwoDimTable &); // dissable, no implicit copying allowed
  TwoDimTable &operator =(const TwoDimTable &); //no implicit copying

  T *aT;
  long size_x, size_y;

public:

  inline TwoDimTable() {
    size_x = size_y = 0;
  }
  inline TwoDimTable(short set_size_x, short set_size_y) {
    Alloc(set_size_x, set_size_y);
  }
  inline ~TwoDimTable() { Dealloc(); }

  //memory managment:
  inline void Alloc(short set_size_x, short set_size_y)
  {
    size_x = set_size_x;
    size_y = set_size_y;
    aT = new T [ size_x * size_y ];
    DEBUG_MSG( DBG_ALLOC, 
               "alocating " << aT
               << " = 2LayerTable[]");
    CHECK_ALLOC(aT);
  }

  inline void Dealloc()
  {
    size_x = size_y = 0;
    if (aT!=NULL)
    {
      DEBUG_MSG( DBG_ALLOC, 
                 "deleting " << aT
                 << " = 2LayerTable[]");
      delete [] aT;
      aT = NULL;
    }
  }

  inline T Get(short i, short j) const
  {
    ASSERT((1<=i) && (i<=size_x) && (1<=j) && (j <= size_y));
    return aT[
              (i-1) * size_y +
              (j-1)
              ];
  }
  
  inline void Set(short i, short j, T set_value)
  {
    ASSERT((1<=i) && (i<=size_x) && (1<=j) && (j <= size_y));
    aT[
       (i-1) * size_y  +
       (j-1)
       ] = set_value;
  }

  //As an alternative, for faster access, you can find the address of the
  //first element of the row by calling the function below, and in the
  //inner-loop do strictly additive pointer arithmetic.
  //Here, you can read or modify the (private) af data, so be careful.
  //    GetAddressOfNIthRow()
  //returns the address of the first element of the ith row of the
  //(n&1 even/odd)th layer of the cost-table.
  //(Note1:The indexing for n and i begins at 1, not 0, although space
  //       for the zero'th element is allocated...
  // Note2:This address is located at an offset of 1 from the
  //       actual beginning of space allocated for that row,
  //       no thanks to the confusing indexing method I'm using.)
  inline T *GetAddressOfIJ(short i, short j) const
  {
    ASSERT((1<=i) && (i<=size_x) && (1<=j) && (j <= size_y));
    return &( aT[
                 (i-1) * size_y +
                 (j-1)
                 ]
              );
  }
}; //class TwoDimTable

} //namespace minrms

#endif //#ifndef _TWO_DIM_TABLE_H
