#ifndef _STREAM_BUGFIXES_H


//The following was necessary because the pgCC compiler as of 6/2002
//seems to have a bug with getline(istream&, string&)

#ifdef GETLINE_UNSUPPORTED
bool SubstituteForGetline(istream &in, string& line);
#define getline(input_stream, line)  SubstituteForGetline(input_stream, line)
#endif


#endif //#ifndef _STREAM_BUGFIXES_H
