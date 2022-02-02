#ifndef _ERR_H
#define _ERR_H


#include <cstdlib>
#include <iostream>
using namespace std;





//The following variables must be defined only once.
//(This usually happens in the "main.cc" file.)
extern const char g_program_name[];
extern const char g_version_string[];
extern const char g_date_string[];


#define ERR( A ) \
    cerr << A << endl, exit(-1)

#define ERR_INTERNAL( A ) \
    cerr << A << "\n(v" << g_version_string << ", file \"" << __FILE__ << "\", line " << __LINE__ << ")" << endl, exit(-1)

//CHECK_ALLOC() is called right after allocating memory with
//new or malloc() to verify that the pointer variable, p, is valid.
#define CHECK_ALLOC(p) \
    if (! (p)) ERR_INTERNAL("Error: Failure to allocate memory\n")

//I don't know if this works yet so I'm commenting it out.
//#define FPRINTF_CHECK(ARGS) \
//if (fprintf(ARGS) < 0) \
// printf("file write failure at line %s in file %s",__LINE__, __FILE__), \
// exit(-1);

#endif //#ifndef _ERR_H
