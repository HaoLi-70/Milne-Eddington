
#ifndef LM_FIT_h
#define LM_FIT_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "MILNE_EDDINGTON.h"
#include "MPI_CTRL.h"
#include "PARAMETER.h"
#include "RANDOM_NUMBER.h"
#include "READ_INPUT.h"
#include "SVD.h"

/*--------------------------------------------------------------------------------*/

extern int LM_FIT(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_PAR *Par, \
        STRUCT_LM *LM);

extern int INVERSION(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
        STRUCT_PAR *Par, STRUCT_LM *LM, STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/

#endif /* LM_FIT_h */
