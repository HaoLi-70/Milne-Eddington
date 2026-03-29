
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "LOG_ERROR.h"
#include "ALLOCATION.h"
#include "LM_FIT.h"
#include "ME_SOLVER.h"
#include "RINPUT.h"

/*--------------------------------------------------------------------------------*/

void FREERAM(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_LM *LM, \
    STRUCT_PARA *Para, STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/


