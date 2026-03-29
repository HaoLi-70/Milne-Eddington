
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include "ALLOCATION.h"
#include "LM_FIT.h"
#include "ME_SOLVER.h"
#include "RINPUT.h"

/*--------------------------------------------------------------------------------*/

extern double Gfactor(double J,double L,double S);

extern double Geffect(double Gu, double Gl, double Ju, double Jl);

extern int INIT_INV(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_LM *LM, \
    STRUCT_PARA *Para, STRUCT_MPI *Mpi, STRUCT_SUBSET *Subset);

/*--------------------------------------------------------------------------------*/


