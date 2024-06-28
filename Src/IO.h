
#ifndef IO_H
#define IO_H

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>
#include "ALLOCATION.h"
#include "MPI_CTRL.h"
#include "READ_INPUT.h"

/*--------------------------------------------------------------------------------*/

extern int rWavelength(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
        STRUCT_MPI *Mpi);

extern int rprofile(STRUCT_INPUT *Input, int coord[], STRUCT_STK *Stk);

extern int **cache_init(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
        STRUCT_STK *Stk, int *status);

extern int cache_write(STRUCT_INPUT *Input, int *coord);

/*--------------------------------------------------------------------------------*/

#endif /* IO_H */
