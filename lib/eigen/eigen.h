#ifndef _EIGEN_H
#define _EIGEN_H

#define NRANSI
#define select NRselect
/* The eigenvalue-eigenvector finding code we use was taken from
   "Numrical Recipees in C".  The following included file ("nr.h") contains
   declarations for all of the functions in the entire Numerical Recipees
   library, even though we only really use one of them (the one in "jacobi.c").
*/
extern "C" {
#include"nr.h"
}
#undef select
#undef NRANSI

#endif //#ifndef _EIGEN_H
