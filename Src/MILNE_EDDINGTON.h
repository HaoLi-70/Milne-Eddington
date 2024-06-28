
#ifndef MILNE_EDDINGTON_h
#define MILNE_EDDINGTON_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include "ALLOCATION.h"
#include "FADDEEVA.h"
#include "PARAMETER.h"
#include "READ_INPUT.h"

/*--------------------------------------------------------------------------------*/

extern int Milne_Eddington(STRUCT_MELINE *Lines, int nline, \
        STRUCT_STK *Stk, STRUCT_FADDEEVA *Fadd, STRUCT_INPUT *Input, 
        bool Deriv);
  
extern int Init_Guess(STRUCT_STK *Stk, STRUCT_PAR *Par, \
        STRUCT_INPUT *Input, STRUCT_LM *LM);

extern int bounds_check(double *Parnew, STRUCT_PAR *Par);

/*--------------------------------------------------------------------------------*/

#endif /* MILNE_EDDINGTON_h */
