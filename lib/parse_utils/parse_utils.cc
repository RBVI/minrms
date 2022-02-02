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
#include <cstdlib> //necessary for strtol() and strtod()
#include <cstring> //necessary for strcpy
#include <fstream>
#include <algorithm>
using namespace std;

#include <global_utils.h>
#include <stl_utils.h>
#include "parse_utils.h"


#ifndef DBG_ARG_FILE_PARSER
//Necessary for some of the DEBUG_MSG() macros.  Just ignore this.
//(See "../global_utils/debug.h" for details.)
#define DBG_ARG_FILE_PARSER 0
#endif



namespace ParseUtils
{

void
NextInt(string::const_iterator& start,
        string::const_iterator  stop,
        long long& dest)
{

  if ((start != stop) && isspace(*start))
  {
    ERR("Error:  Text containing integer was preceded by unnexpected\n"
        "        whitespace.  (This could be a bug in the program.)");
  }

  //First, I have to copy this text into a c_string to get around strtol()'s
  //weird syntax: strings only grant "const char *" access to their contents,
  //but strtol() requires "char *" arguments.
  string temp_str(start, stop);
  char *ac_temp_str = new char [temp_str.size() + 1];
  strcpy(ac_temp_str, temp_str.c_str());
  char *pstart = ac_temp_str;
  char *pstop;
  long long result = strtol(pstart, &pstop, 10);

  //If there was at least some valid text was read, then go ahead
  //and use the value returned by strtol(), by setting dest = result.
  //(otherwise discard the value).

  if ((pstop == pstart+1) && (*start == '-')) //first, we safeguard against "-"
    pstop = pstart; //because strtol() incorrectly skips past the '-' in "-"

  if (pstop != pstart) //If there was some valid text, use the value.
  {
    dest = result;
    //...and skip past the text that was read.
    start = start + (pstop-pstart);
  }
  delete [] ac_temp_str;
} //NextInt()



void
NextInt(string::const_iterator& start,
        string::const_iterator  stop,
        long long& dest,
        string terminators)
{
  string token;
  NextToken(start, stop, token, terminators);

  if (token != "")
    //dest = atol(token.c_str()); This works, but it doesn't check the format.
  {
    string::const_iterator got_this_far = token.begin();
    NextInt(got_this_far,
            token.end(),
            dest);

    if (got_this_far != token.end())
      ERR("Error: \"" << token << "\" is not a valid integer.");
  }
} //NextInt()




//StrToReal() is a C++-string compatible version of the standard C-library's
//strtod().  Instead of pointers, StrToReal() takes two arguments
//which are string::const_iterators.
//These two iterators point to the beginning and ending characters demarking
//an interval within a larger string where a number is stored.
//This number contained in the string will be returned to the user.
//Afterwards "stop" will point to the character after the last
//character used in the conversion.

long double
StrToReal(string::const_iterator  start,
          string::const_iterator& stop)
{

  //First, I have to copy this text into a c_string to get around strtod()'s
  //weird syntax: strings only grant "const char *" access to their contents,
  //but strtod() requires "char *" arguments.
  string temp_str(start, stop);
  char *ac_temp_str = new char [temp_str.size() + 1];
  strcpy(ac_temp_str, temp_str.c_str());
  char *pstart = ac_temp_str;
  char *pstop;

  #ifdef STRTOLD_UNSUPPORTED
  long double result = strtod(pstart, &pstop);
  #else
  long double result = strtold(pstart, &pstop);//Useful but not standard ANSI C
  #endif
  //now pstop points past the last valid char

  stop = start + (pstop - pstart);
  delete [] ac_temp_str;

  return result;
} //StrToReal()




void
NextReal(string::const_iterator& start,
         string::const_iterator  stop,
         long double& dest)
{
  if ((start != stop) && isspace(*start))
  {
    ERR("Error:  Text containing floating point number was preceded by\n"
        "        unnexpected whitespace.  (This could be a bug in the program.)");
  }
  string::const_iterator p = stop;
  long double result = StrToReal(start, p);

  //If there was at least some valid text was read, then go ahead
  //and use the value returned by StrToReal(), by setting dest = result.
  //(otherwise discard the value).

  if ((p == start+1) && (*start == '-')) //first, we safeguard against "-"
    p = start; //because StrToReal() incorrectly skips past the '-' in "-"
  else if ((p == start+1) && (*start == '-')) //safeguard against "."
    p = start; //because StrToReal() incorrectly skips past the '.' in "."
  else if ((p == start+2) && (*start == '-') && (*(start+1) == '.'))
    //safeguards against "-."
    p = start; //because StrToReal() incorrectly skips past "-."

  if (p != start) //If there was some valid text, use the value.
  {
    dest = result;
    //...and skip past the text that was read.
    start = p;
  }

} //NextReal()




void
NextReal(string::const_iterator& first,
         string::const_iterator  last,
         long double& dest,
         string terminators)
{
  string token;
  NextToken(first, last, token, terminators);

  if (token != "")
    //dest = atof(token.c_str()); This works, but it doesn't check the format.
  {
    string::const_iterator got_this_far = token.begin();
    NextReal(got_this_far,
             token.end(),
             dest);

    if (got_this_far != token.end())
      ERR("Error: \"" << token << "\" is not a valid floating point number.");
  }
} //NextReal()




void ConvertArgvToVectorOfStrings(int argc,
                                  char **argv,
                                  vector<string>& dest)
{
  dest.resize(argc);
  assert(argv);
  for (int i=0; i < argc; ++i)
  {
    assert(argv[i]);
    dest[i].assign(argv[i]);
  }
}


void
HandleArgFiles(vector<string> const& source_args,
               vector<string> & dest_args,
               string config_file_flag,
               string comment_indicator,
               string merge_lines_indicator,
               string default_args_filename,
               bool whitespace_delimited)
{
  assert(source_args.size() >= 1); //<-should at least contain the name of the
                                   //  executable as the first argument.
                                   //  (as with argv[0] in C)
  #ifdef DEBUG
  #if DBG_ARG_FILE_PARSER
  DEBUG_MSG(DBG_ARG_FILE_PARSER, "Initial argument list (before any processing):");
  DisplayVectorOfStrings(source_args);
  #endif // #if DBG_ARG_FILE_PARSER
  #endif // #ifdef DEBUG

  dest_args = source_args;

  // ***  Read in the arguments from the 
  // ***  default configuration file (if there is one).
  vector<string> default_args;
  if (default_args_filename != "")
  {
    ifstream default_args_file(default_args_filename.c_str(), ios::in);
    if (! default_args_file)
      cerr <<
    "Warning: Unable to open the default arguments file:\""
       << default_args_filename << "\"\n"
    "         No default arguments loaded." << endl;
    else
      ArgFileToArgs(default_args_file,
            default_args,
            whitespace_delimited,
            comment_indicator,
            merge_lines_indicator);
  }

  //...and insert the arguments from the default configuration file right after
  //right after first argument in the list (the name of the executable).
  //(The first argument in the array of arguments is special.  By convention,
  // the first argument should contain the name of the program,
  // as with argv[0] in C.  We want to preserve this convention.)
  dest_args.insert(dest_args.begin() + 1,
                   default_args.begin(),
                   default_args.end());

  //If there was no config_file_flag specified, we are done.
  if (config_file_flag == "")
    return;

  //Otherwise, scan through the list of arguments.
  //If any of them match the text in the config_file_flag variable,
  //then the following argument will be the name of a
  //configuration file.  This file stores additional arguments.
  //Read in the arguments from that file,
  //and insert them into the argument-list at this point.
  //    Recursion-capable:
  //If one of the arguments in the file is also a config_file_flag,
  //then the process is repeated.
  //Warning: there is no safeguarding against infinite loops.

  DEBUG_MSG(DBG_ARG_FILE_PARSER, "At the beginning of loop inside HandleArgFiles()");

  for (vector<string>::iterator p = dest_args.begin();
       p != dest_args.end();
       ++p)
  {
    //If any arguments do match config_file_flag,
    if ((*p) == config_file_flag)
    {
      //insert arguments from the file indicated.
      int insert_location = p - dest_args.begin();
      ++p;
      if (p == dest_args.end())
        ERR("Error: the \"" << config_file_flag <<
            "\" flag must be followed by the name of\n"
            "      a configuration file.");

      ifstream arg_file((*p).c_str());
      if (! arg_file)
        ERR("Error: Cannot open file \"" << (*p) << "\" for reading.");

      //Read in the arguments from config_file into config_file_args
      vector<string> args_in_file;
      ArgFileToArgs(arg_file,
            args_in_file,
            whitespace_delimited,
            comment_indicator,
            merge_lines_indicator);

      //erase the arguments storing the config_file_flag and the file's name.
      dest_args.erase(p - 1, p + 1);

      //and insert them at this point.
      dest_args.insert(dest_args.begin() + insert_location,
                       args_in_file.begin(),
                       args_in_file.end());

      //Then set p to point to the first argument we inserted.
      p = dest_args.begin() + insert_location;
    }
  } //for (p == dest_args.begin(); p != dest_args.end(); ++p)

  #ifdef DEBUG
  #if DBG_ARG_FILE_PARSER
  DEBUG_MSG(DBG_ARG_FILE_PARSER, "Final argument list (before parsing):");
  DisplayVectorOfStrings(dest_args);
  #endif // #if DBG_ARG_FILE_PARSER
  #endif // #ifdef DEBUG

} //HandleArgFiles()




//SplitUpUsingWhitespace() splits up a string containing whitespace into
//multiple strings, as long as the whitespace in question
//is not enclosed in quotes.
//(Note: Text located in nested qutes like "Sean says 'God Bless America.'"
//       is not split up into different arguments.))
//The resulting strings are saved in the "dest" argument.
//This function also deletes leading and trailing whitespace from the
//strings it generates.
static void
SplitUpUsingWhitespace(string const& source, vector<string>& dest)
{
  dest.clear();
  dest.push_back(string());
  char quote_char = '\0';

  DEBUG_MSG(DBG_ARG_FILE_PARSER, "source=\"" << source << "\"");
  
  string::const_iterator p = source.begin();
  SkipWhiteSpace(p, source.end()); //skip leading spaces
  while (p != source.end())
  {
    //#ifdef DEBUG
    //string debug_msg(p, source.end());
    //DEBUG_MSG(DBG_ARG_FILE_PARSER,
    //          "[p,source.end())=\""
    //          << debug_msg << "\"\n"
    //          "                   dest.back()=\""
    //          << dest.back() << "\"");
    //#endif //#ifdef DEBUG

    //If the characters is whitespace, 
    //and is not enclosed in either kind of quote
    //split this argument up into two arguments,
    if (isspace(*p) && (quote_char == '\0'))
    {
      dest.push_back(string());
      SkipWhiteSpace(p, source.end());
    }
    else
    {
      if ((quote_char != '\0') && ((*p) == quote_char))
      {
        quote_char = '\0';
      }
      else if ((quote_char == '\0') && (((*p) == '\'') || ((*p) == '"')))
        quote_char = (*p);
      else
        dest.back().push_back(*p);

      ++p;
    }
  }

  //If there was any whitespace left after the last argument,
  //a new argument will have been created, but it will not
  //contain anything.  See the clause following
  //"if (isspace(*p) && (quote_char == '\0'))"
  //above.  So we delete this text now:
  //(We only have to worry about this happening once at the end.)
  if ((dest.size() > 0) && (dest.back() == ""))
    dest.erase(dest.end()-1);

  if (quote_char != '\0')
    ERR("Error: unmatched quotes in line: \""
        << source << "\"");

  //The following check takes care of strings containing only whitespace.
  if ((dest.size() == 1) && (dest[0] == ""))
    dest.clear();
} //SplitUpUsingWhitespace()




void
ArgFileToArgs(istream& in,
              vector<string>& dest,
              bool whitespace_delimited,
              string comment_indicator,
              string merge_lines_indicator)
{
  dest.clear();

  //1) Load in the lines of the file into a vector of strings,
  string line;
  vector<string> vLines;
  while (getline(in, line))
  {
    dest.push_back(line);
  }
  #ifdef DEBUG
  #if DBG_ARG_FILE_PARSER
  DEBUG_MSG(DBG_ARG_FILE_PARSER, "1) Read an argument file.  Lines:");
  DisplayVectorOfStrings(dest);
  #endif // #if DBG_ARG_FILE_PARSER
  #endif // #ifdef DEBUG


  //2) Now, we have a vector of strings containing text from
  //   separate lines in the file.
  //   Delete all text after a comment indicator in each line.
  for(vector<string>::iterator pLine = dest.begin();
      pLine != dest.end();
      ++pLine)
  {
    if (comment_indicator != "")
    {
      string::iterator pc = search((*pLine).begin(),
                                   (*pLine).end(),
                                   comment_indicator.begin(),
                                   comment_indicator.end());
      (*pLine).erase(pc, (*pLine).end());
    }
  } //handle comment_indicators
  #ifdef DEBUG
  #if DBG_ARG_FILE_PARSER
  DEBUG_MSG(DBG_ARG_FILE_PARSER, "2) Deleted comments.  Lines:");
  DisplayVectorOfStrings(dest);
  #endif // #if DBG_ARG_FILE_PARSER
  #endif // #ifdef DEBUG


  //3) After deleting commented out text,
  //   merge lines ending in a merge_lines_indicator.
  if (merge_lines_indicator != "")
  {

    vector<string>::iterator  pPrevLine = dest.end();

    for(vector<string>::iterator pLine = dest.begin();
        pLine != dest.end();
        ++pLine)
    {
      if (pPrevLine != dest.end())
      {
        string::iterator pc = search((*pPrevLine).begin(),
                                     (*pPrevLine).end(),
                                     merge_lines_indicator.begin(),
                                     merge_lines_indicator.end());
        if (pc != (*pPrevLine).end())
        {
          assert(pc + merge_lines_indicator.size() <= (*pPrevLine).end());

          //delete the merge_lines_indicator from the previous line.
          (*pPrevLine).erase(pc, (*pPrevLine).end());
          //append the contents of the current line to the previous line.
          (*pPrevLine).insert((*pPrevLine).end(),
                              (*pLine).begin(),
                              (*pLine).end());
          //afterwards, delete the redundant current line
          dest.erase(pLine); //this invalidates pLine
          pLine = pPrevLine; //this re-validates pLine,
        }
      }
      pPrevLine = pLine;
    } //loop over all strings in dest
    #ifdef DEBUG
    #if DBG_ARG_FILE_PARSER
    DEBUG_MSG(DBG_ARG_FILE_PARSER, "3) Merged lines.  Lines:");
    DisplayVectorOfStrings(dest);
    #endif // #if DBG_ARG_FILE_PARSER
    #endif // #ifdef DEBUG
  } //handle merge_lines_indicators


  //4) Now, split up any lines into separate
  //   arguments if they contain whitespace (if necessary).
  if (whitespace_delimited)
  {
    vector<string>::iterator pLine     = dest.begin();
    //vector<string>::iterator pPrevLine = dest.end();

    while (pLine != dest.end())
    {
      vector<string> args_in_line;
      SplitUpUsingWhitespace((*pLine), args_in_line);
      #ifdef DEBUG
      #if DBG_ARG_FILE_PARSER
      DEBUG_MSG(DBG_ARG_FILE_PARSER, "3.5) Args in line #"
                << pLine - dest.begin() + 1);
      DisplayVectorOfStrings(args_in_line);
      #endif // #if DBG_ARG_FILE_PARSER
      #endif // #ifdef DEBUG
      //delete pLine from dest, and insert args_in_line
      if (args_in_line.size() > 0)
      {
        long insert_location = pLine - dest.begin();
        dest.insert(pLine,
                    args_in_line.begin(),
                    args_in_line.end());
        //Delete the argument that used to be stored in "pLine"
        dest.erase(dest.begin() + insert_location + args_in_line.size());
        pLine = dest.begin() + insert_location + args_in_line.size();
      }
      else
      {
        int i = pLine - dest.begin();
        dest.erase(pLine);
        pLine = dest.begin() + i;
      }
      #ifdef DEBUG
      #if DBG_ARG_FILE_PARSER
      DEBUG_MSG(DBG_ARG_FILE_PARSER, "3.6) Right after insertion of args in line:#"
                << pLine - dest.begin() + 1);
      DisplayVectorOfStrings(dest);
      #endif // #if DBG_ARG_FILE_PARSER
      #endif // #ifdef DEBUG

      //pPrevLine = pLine;
    } //loop over all strings in dest
    #ifdef DEBUG
    #if DBG_ARG_FILE_PARSER
    DEBUG_MSG(DBG_ARG_FILE_PARSER, "4) Split up using whitespace.  tokens:");
    DisplayVectorOfStrings(dest);
    #endif // #if DBG_ARG_FILE_PARSER
    #endif // #ifdef DEBUG
  } //handle whitespace inside arguments

} //ArgFileToArgs();









void
GetUntil(istream &in,
         string& dest,
         string terminator_chars)
{
  DeleteContents(dest);
  if (! in)
    ERR("(GetUntil)\nError reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  {
    if (::BelongsTo(c, terminator_chars.begin(), terminator_chars.end()))
    {
      in.putback(c);
      break;
    }
    else
      dest.push_back(c);
  }
}



void
GetlineUntil(istream &in,
             string& dest,
             string terminator_chars)
{
  DeleteContents(dest);
  if (! in)
    ERR("(GetlineUntil)\nError reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  {
    if (::BelongsTo(c, terminator_chars.begin(), terminator_chars.end()))
    {
      //in.putback(c);  No. getline() removes the terminator from the stream
      break;
    }
    else
      dest.push_back(c);
  }
}



void
Skip(istream &in, string chars_to_skip)
{
  if (! in)
    ERR("(Skip)\nError reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  {
    if (! ::BelongsTo(c, chars_to_skip.begin(), chars_to_skip.end()))
    {
      in.putback(c);
      break;
    }
  }
}

void
SkipWhiteSpace(istream &in)
{
  if (! in)
    ERR("(SkipWhiteSpaces)\nError reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  {
    if ((! isspace(c)) && in)
    {
      in.putback(c);
      break;
    }
  }
}

void
SkipWhiteSpaceExceptNewlines(istream &in)
{
  if (! in)
    ERR("(SkipWhiteSpacesExceptNewlines)\n"
        "Error reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  { // skip past white space before next line
    if ((! isspace(c)) || (c == '\n'))
    {
      in.putback(c);
      break;
    }
  }
}


//same as above, except it skips everything _but_ whitespace
void
SkipNonWhiteSpace(istream &in)
{
  if (! in)
    ERR("(SkipNonWhiteSpace)\nError reading file: File ends prematurely.");
  char c;
  //  bool all_numeric = true;
  while(in.get(c))
  {
    if (isspace(c))
    {
      in.putback(c);
      break;
    }
    //    else if (! IS_NUMERIC(c))
    //      all_numeric = false;
  }
  //  return all_numeric;
}





// Skips all whitespace and comments in the "in" stream.
void
SkipWhiteSpaceAndComments(istream &in,
                          char comment_character)
{
  if (! in) {
    ERR("(SkipWhiteSpaceAndComments)\n"
        "Error reading file: File ends prematurely.");
  }
  char c;
  while(in.get(c))
  {
    if (c == comment_character)
    {
      string ignore;
      getline(in, ignore);
    }
    else if ((! isspace(c)) && in)
    {
      in.putback(c);
      break;
    }
  }
}



} //namespace ParseUtils




