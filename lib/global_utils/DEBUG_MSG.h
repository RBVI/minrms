#ifndef _DBG_MESSAGE_H
#define _DBG_MESSAGE_H

#include <iostream>
using namespace std;


 /********************************** debug.h **********************************
  * Utilities for catching bugs and printing out messages that are            *
  * that are conditionally compiled, and do not influence the the normal      *
  * working behavior of the program.                                          *
  *   DEBUG_MSG(), DEBUG_EXEC()                                               *
  * ------------------ Using debug-flags in code: -------------------------
  * -If a debug-flag is set to 1, this triggers
  *  diagnostic messages to be printed at any places in the code
  *  where an instance of DEBUG_MSG() was invoked using that flag as 
  *  an argument.
  * -Additionally, if any of the flags in "debug_flags.h" are defined to be 1,
  *  this forces the 'NDEBUG' symbol to also be undefined.
  * -Lastly, if all debug flags are turned off (set to 0), all of this
  *  diagnostic code is actually deleted by the preprocessor elliminating
  *  the conditional branching which could modestly slow down the program.
  ****************************************************************************/

// --- DEBUG_MSG(type, msg) ---
// -Prints out the message contained in 'msg' if 'type' is defined.
// -The argument 'type' should be one of the debug-flags listed above.
//  (Note: To print a message only if the 'type' flag is _not_ defined call
//  "DEBUG_MSG( (! type), msg )", but this is not the intended use,
//   and it does not work if NDEBUG is defined.)
// -If none of the debug-flags are defined then calls to DEBUG_MSG()
//   are deleted by the preprocessor  (optimize for speed).


#ifndef NDEBUG
/* If I'm not mistaken, the following line is a statement, not an expression
   so it cannot be substituted wherever a function call can.  If
   you are having trouble compiling a line with a DEBUG_MSG() macro
   pull it out of the expression it is defined in.  This shouldn't be
   a problem in most cases.
if (type!=0) cerr << #type << " msg: " << msg << endl
*/
/*
I am more comfortable with the following version because it is a C expression
-not a statement, and can be placed anywhere a function-call can.
It does not work on all compilers:
((type!=0) ? (cerr << #type << " msg " << msg << endl) : 0 )
*/
     /* wait, this should work: 
        I believe this is a valid expression with a return value (0).*/


#define DEBUG_MSG( type, msg ) \
if (type) cerr << #type << " msg: " << msg << endl, 0
//cerr <<" MSG: " << msg << endl, 0
#else
#define DEBUG_MSG( type, msg )  //else, do nothing
#endif




// --- DEBUG_EXEC(type, statement) ---
// -Executes the statment or block indicated by 'statement', so long as 'type'
//   is defined.
// -By convention, additional code contained in the 'statement' should not
//  really modify the working behavior of the program, only help display
//  diagnostic information.
// -The argument 'type' should be one of the debug-flags listed above.
//  (Note: To execute a statement only if the 'type' flag is _not_ defined call
//  "DEBUG_EXEC( (! type), statement )", but this is not the intended use,
//   and it does not work if DEBUG is not defined.)
// -If none of the debug-flags are defined then calls to DEBUG_EXEC()
//   are deleted by the preprocessor  (optimized for speed).

#ifndef NDEBUG
/*#define DEBUG_EXEC( type, statement )   if (type!=0) (statement)*/
#define DEBUG_EXEC( type, statement )   if (type!=0) (statement), 0
#else
#define DEBUG_EXEC( type, statement ) //else, do nothing
#endif

#endif /* #ifndef _DBG_MESSAGE_H */
