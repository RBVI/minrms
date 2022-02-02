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


#include <cctype> //required for "isalpha()"
using namespace std;

#include "pdb_id.h"

namespace minrms
{

bool
PDBresID::operator < (const PDBresID& B) const
{
  return
    ( ((this->chainId) < (B.chainId)) ?
      true                                   
      :                                      
      ( ((this->chainId) > (B.chainId)) ?
      false                                
      :                                      
      ( ((this->seqNum) < (B.seqNum)) ?
        true                                   
        :                                      
        ( ((this->seqNum) > (B.seqNum)) ?
            false                                
            :
            ( ((this->insertCode) < (B.insertCode)) ?
              true : false                         
            )
          )                                       
        )
      )
    );
}

bool
PDBresID::operator <= (const PDBresID& B) const
{
  return (((*this) < B) || (*this == B));
}


bool
PDBresID::operator > (const PDBresID& B) const
{
  return ! ((*this) < B);
}

bool
PDBresID::operator >= (const PDBresID& B) const
{
  return (((*this) > B) || (*this == B));
}


bool
PDBresID::operator == (const PDBresID& B) const
{
  return ( (this->chainId    == B.chainId) &&
           (this->seqNum     == B.seqNum) &&
           (this->insertCode == B.insertCode) );
}


bool
PDBresID::operator != (const PDBresID& B) const
{
  return !(*this == B);
}


ostream& operator << (ostream&  out_file, PDBresID  res_id)
{
  out_file << res_id.seqNum;

  if (isalpha(res_id.insertCode))
    out_file << res_id.insertCode;
  else if (res_id.insertCode == '\0')
    out_file << "\\0";
  else if (res_id.insertCode != ' ')
    out_file << "ASC(" << (int)res_id.insertCode << ")";

  if (isalpha(res_id.chainId))
    out_file << "." << res_id.chainId;
  else if (res_id.chainId == '\0')
    out_file << ".\\0";
  else if (res_id.chainId != ' ')
    out_file << ".ASC(" << (int)res_id.chainId << ")";

#if 0
  out_file << "\'";
  if (isalpha(res_id.chainId) || (res_id.chainId == ' '))
    out_file << res_id.chainId;
  else if (res_id.chainId == '\0')
    out_file << "\\0";
  else
    out_file << "ASC(" << (int)res_id.chainId << ")";
    
  out_file << "\'," << res_id.seqNum << ",\'";

  if (isalpha(res_id.insertCode) || (res_id.insertCode == ' '))
    out_file << res_id.insertCode;
  else if (res_id.insertCode == '\0')
    out_file << "\\0";
  else
    out_file << "ASC(" << (int)res_id.insertCode << ")";
  out_file << "\'";
#endif //#if 0

  return out_file;
}

} //namespace minrms
