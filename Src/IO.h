
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>
#include "ALLOCATION.h"
#include "ME_SOLVER.h"
#include "MPI_INIT.h"
#include "RINPUT.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

/*--------------------------------------------------------------------------------*/

extern int rWavelength(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
    STRUCT_MPI *Mpi);

extern int rProfile(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
    STRUCT_SUBSET *Subset);

extern int rProfilesll(STRUCT_INPUT *Input, STRUCT_STK *Stk);

extern int PixelMV(STRUCT_INPUT *Input, STRUCT_SUBSET *Subset);

extern int CACHE_INIT(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_STK *Stk);

extern int WRITE_RESULT(STRUCT_INPUT *Input, STRUCT_SUBSET *Subset, \
    STRUCT_STK *Stk);

extern int CLOSE_FILES(void);

/*--------------------------------------------------------------------------------*/
