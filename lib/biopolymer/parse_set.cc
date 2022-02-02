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


#include <cassert> //defines "assert()"
using namespace std;

#include <parse_utils.h>
#include "parse_set.h"

using namespace ParseUtils;

namespace minrms
{


void
ParseChimeraSet(string::const_iterator& start,
                string::const_iterator  stop,
                Biopolymer const& m,
                vector<ClosedInterval>& dest,
                string terminators)
{
  string token;
  string::const_iterator pc = start;
  ParseUtils::NextToken(pc, stop, token, terminators);
  string::const_iterator tbegin = token.begin();
  string::const_iterator tend   = token.end();
  ParseChimeraSet(tbegin, tend, m, dest);
  if (tbegin != tend)
    ERR("Error: Syntax error in residue specifiers:\n"
        "\"" << token << "\"");
  else
    start = pc;
}



static void GetResidues(string::const_iterator& start,
                        string::const_iterator  stop,
                        int&  first_SeqNum,
                        char& first_InsertCode,
                        bool& first_InsertCode_is_begin,
                        bool& first_res_is_begin,
                        int&  last_SeqNum,
                        char& last_InsertCode,
                        bool& last_InsertCode_is_end,
                        bool& last_res_is_end)
{
  if (start == stop)
    ERR("Error: Invalid residue specifiers: missing text.");

  first_InsertCode_is_begin = last_InsertCode_is_end = true;

  string::const_iterator pc = start;
  //First, read in the residue at the beginning of this interval:
  if ((isdigit(*pc)) || ((*pc) == '-') || (*pc == '*'))
  {
    if (*pc == '*')
    {
      first_res_is_begin = true;
      ++pc;
    }
    else
    {
      first_res_is_begin = false;
      if (! isdigit(*pc))
        ERR("Error in selection syntax.");

      NextInt(pc, stop, first_SeqNum);
      //If the character right after the integer is not - . or , 
      //then it's a residue insert-code.
      if (! ((pc == stop) || isspace(*pc) || (*pc == '.') || (*pc == ',') ||
             (*pc == '-')))
      {
        first_InsertCode = *pc;
        first_InsertCode_is_begin = false;
        ++pc;
      }
    }
    if ((pc != stop) && (*pc == '-'))
    {
      ++pc;
      if (*pc == '*')
      {
        last_res_is_end = true;
        ++pc;
      }
      else if ((isdigit(*pc)) || ((*pc) == '-'))
      {
        last_res_is_end = false;
        NextInt(pc, stop, last_SeqNum);
        if (! ((pc == stop) || isspace(*pc) || (*pc == '.') || (*pc == ',')))
        {
          last_InsertCode = *pc;
          last_InsertCode_is_end = false;
          ++pc;
        }
      }
    }
    else
    {
      //If no - is present, the first residue is the last residue.
      last_res_is_end = first_res_is_begin;
      last_SeqNum = first_SeqNum;
      last_InsertCode_is_end = first_InsertCode_is_begin;
      last_InsertCode = first_InsertCode;
    }
  } // if ((isdigit(*pc)) || ((*pc) == '-') || (*pc == '*'))

#if 0
  //Special case: A set beginning with '.' includes all the residues
  //              in each chain specified after the dot.  So use all residues.
  //              This should not be necessary.
  else if (*pc == '.')
  {
    first_res_is_begin = true;
    last_res_is_end = true;
    return;
  }
#endif //#if 0

  start = pc;
} //GetResidues()



static void GetChains(string::const_iterator& start,
                      string::const_iterator  stop,
                      char& first_chain, bool& first_chain_is_begin,
                      char& last_chain,  bool& last_chain_is_end)
{
  string::const_iterator pc = start;
  if (*pc == '-')
    ERR("Syntax error. Chain-specifiers cannot begin with '-'.");

  //First, read in the residue at the beginning of this interval:
  if ((pc == stop) || isspace(*pc) || (*pc == ','))
  {
    last_chain_is_end = first_chain_is_begin = true;
  }
  else
  {
    last_chain_is_end = first_chain_is_begin = false;
    if (*pc == '*')
    {
      first_chain_is_begin = true;
      ++pc;
    }
    else
    {
      first_chain = *pc;
      ++pc;
    }
    if ((pc != stop) && (*pc == '-'))
    {
      ++pc;
      if ((pc == stop) || isspace(*pc) || (*pc == ','))
      {
        string set_str(start, stop);
        ERR("Error: Syntax error in residue specifiers:\n"
            "\"" << set_str << "\"");
      }
      else if (*pc == '*')
      {
        last_chain_is_end = true;
        ++pc;
      }
      else
      {
        last_chain = *pc;
        ++pc;
      }
    }
    else
    {
      //If no - is present, the first chain is the last chain.
      last_chain_is_end = first_chain_is_begin;
      last_chain = first_chain;
    }
  } //else clause for "if ((pc == stop) || isspace(*pc) || (*pc == '*'))"

  start = pc;

#if 0
  if ((pc == stop) || isspace(*pc) || (*pc == ','))
  {
    return; //normal termination
  }
  else
  {
    string set_str(start, stop);
    ERR("Error: Syntax error in residue specifiers:\n"
        "\"" << set_str << "\"");
  }
#endif //#if 0

} //GetChains()



//Examples, and what they mean:
// "100"         all residues whose seqNums are 100
// "100-150.A"   all residues in chain A whose seqNums are in [100,150]
// "*-150.A"     all the residues in chain A whose seqNums are up to 150
// "150-*.A"     all the residues in chain A whose seqNums are at least 150
// "*.A"         all the residues in chain A
// ".A"           "   "     "     "    "   "
// "100-150.A-C" residues 100-150 in chains A through C
// "100-150.*"   residues 100-150 in all chains
// "100-150."       "      "   "    "    "   "
// "*.*"         the entire molecule
  

void ParseChimeraSet(string::const_iterator& pc,
                     string::const_iterator  stop,
                     Biopolymer const& m,
                     vector<ClosedInterval>& dest)
{
  if (m.size() == 0)
    ERR("ERROR:  Molecule is allready empty.\n"\
        "        Cannot find subsets of residues from this molecule.\n");

  dest.clear(); //make sure dest is empty

  //The following variables will indicate which residues and chains
  //belong to the set.
  const char UNSPECIFIED_INSERT_CODE = ' ';
  const char UNSPECIFIED_CHAIN_ID = ' ';

  char first_chain = UNSPECIFIED_CHAIN_ID; //Selects which chains contain
                                           //residues in the set.
  bool first_chain_is_begin = true; //true means: start at the first chain
                                    //in the molecule (overrides "first_chain")

  char last_chain = UNSPECIFIED_CHAIN_ID;
  bool last_chain_is_end = true; //true means: stop at the last chain
                                 //in the molecule (overrides "last_chain")

  int  first_SeqNum;//Specifies the range of residues from
                    //each chain belonging to the set
  char first_InsertCode = UNSPECIFIED_INSERT_CODE;
  bool first_InsertCode_is_begin = true;
  bool first_res_is_begin = true; //start at the first residue in each chain
                                  //(overrides first_SeqNum and first_InsertCode)
  int  last_SeqNum;  
  char last_InsertCode = UNSPECIFIED_INSERT_CODE;
  bool last_InsertCode_is_end = true;
  bool last_res_is_end = true; //stop at the last residue in each chain
                               //(overrides last_SeqNum and last_InsertCode)


  SkipWhiteSpace(pc, stop);
  if (*pc == '.')
  {
    ++pc;
    GetChains(pc, stop,
              first_chain, first_chain_is_begin,
              last_chain, last_chain_is_end);
  }
  else
  {
    GetResidues(pc, stop,
                first_SeqNum, first_InsertCode,
                first_InsertCode_is_begin, first_res_is_begin,
                last_SeqNum, last_InsertCode,
                last_InsertCode_is_end, last_res_is_end);
    if (*pc == '.')
    {
      ++pc;
      GetChains(pc, stop,
                first_chain, first_chain_is_begin,
                last_chain, last_chain_is_end);
    }
  }


  //Now that we know which residues from each chain are in the set,
  //as well as which chains, 
  //figure out explicit endpoints for all the intervals in the set.

  if (first_chain_is_begin)
    first_chain = m[0].id.chainId;

  if (last_chain_is_end)
    last_chain = m[m.size()-1].id.chainId;

  Biopolymer::const_iterator pS = m.begin();
  while ((pS != m.end()) && ((*pS).id.chainId != first_chain))
    ++pS;

  if (pS == m.end())
    ERR("Chain '" << first_chain << "' not found in molecule.");

  //loop over all chains in the molecule
  bool passed_last_chain = false;
  while ((! passed_last_chain) && (pS != m.end()))
  {
    char current_chain = (*pS).id.chainId;

    //---- Find the first residue in the interval. ---

    if (! first_res_is_begin)
    {
      //skip forward in the chain until we get to the "first" residue
      while(
            (pS != m.end())
            &&
            ((*pS).id.chainId == current_chain)
            &&
            (((*pS).id.seqNum != first_SeqNum) ||
             (((*pS).id.insertCode != first_InsertCode) &&
              (! first_InsertCode_is_begin)))
            )
        ++pS;

      if ((pS == m.end())
          ||
          ((*pS).id.chainId != current_chain))
      {
        ERR("residue "  << first_SeqNum << first_InsertCode <<
            " not found in chain '" << current_chain << "' of molecule.\n"
            "(Or residues are specified out of order.)");
      }
    } //if (! first_res_is_begin)


    ClosedInterval interval;
    //Save the chain, InsertCode, and SeqNum of the first residue.
    interval.first = (*pS).id;

    //---- Now, find the last residue in the interval. ---
    //---- (This is harder, -Either that or i can't program worth crap.) ----

    if (! last_res_is_end)
    {
      //Find the last residue for this interval
      //First, skip to the residues up until this sequenceNum.
      while(
          (pS != m.end())
          &&
          ((*pS).id.chainId == current_chain)
          &&
          ((*pS).id.seqNum != last_SeqNum)
          )
        ++pS;

      //Check to see if we went too far.
      if ((pS == m.end())
          ||
          ((*pS).id.chainId != current_chain))
        ERR("residue "  << last_SeqNum << last_InsertCode <<
            " not found in chain '" << current_chain << "' of molecule.\n"
            "(Or residues are specified out of order.)");

      //If not, skip to the last residue with this seqNum
      //and an insertCode no later than last_InsertCode
      while ((pS != m.end())
             &&
             ((*pS).id.chainId == current_chain)
             &&
             ((*pS).id.seqNum == last_SeqNum)
             &&
             (last_InsertCode_is_end ||
              ((*pS).id.insertCode != last_InsertCode)))
        ++pS;

      //If we skipped just past the end of this seqNum, back up one.
      if (last_InsertCode_is_end)
      {
        assert((pS == m.end()) ||
               ((*pS).id.chainId != current_chain) ||
               ((*pS).id.seqNum != last_SeqNum));
        --pS;
      }
      else
      {
        //Else, check to see if we actually found the residue.
        if ((pS == m.end()) ||
            ((*pS).id.chainId != current_chain) ||
            ((*pS).id.seqNum != last_SeqNum))
          ERR("residue "  << last_SeqNum << last_InsertCode <<
              " not found in chain '" << current_chain << "' of molecule.\n"
              "(Or residues are specified out of order.)");
      } //else for "if (last_InsertCode_is_end)"

    } //if (! last_res_is_end)
    else
    {
      assert((*pS).id.chainId == current_chain);
      //Find the last residue for this chain
      while((pS != m.end())
            &&
            ((*pS).id.chainId == current_chain))
        ++pS;
      --pS;
    }

    //save it.
    interval.last = (*pS).id;
    
    //store this interval in the list.
    dest.push_back(interval);

    //if there are multiple chains, skip to the end of this one
    //so the next time through the loop, we start with a new chain.
    while((pS != m.end())
          &&
          ((*pS).id.chainId == current_chain))
      ++pS;

    DEBUG_MSG(DBG_INTERVAL_LISTS,
              "Interval parsed with chain=" << current_chain << " ["
              << interval.first << ".." << interval.last << "]");

    //Figure out when we can stop looping over chains
    passed_last_chain = (current_chain == last_chain);
  } //loop over all chains
} //ParseChimeraSet()



} //namespace minrms





