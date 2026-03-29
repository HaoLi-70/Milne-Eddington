
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "ME_SOLVER.h"
#include "MPI_INIT.h"
#include "RANDOM_NUMBER.h"
#include "RINPUT.h"
#include "SVD.h"

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Levenberg_Marquardt{

    // the indx array for LU decomposition.
    int *indx, niter, verboselv, nruns;

    // the Hessian matrix, the b vacror and the solutions
    STRUCT_MATRIX Hessian, Hessian_new, Regul_H, V;
    
    double *Jacfvec, *Sol, *Regul_J, *W;

    // the damping factor, and values to update the facrot.
    double Damp, Lam_reject, Lam_accept, Penalty, Criteria, Threshold;

    // the limits of the damping factor.
    double Lam_Lim[2];

    double ratio, chisq;
  
    bool Regul_flag, HMI_REF;

    STRUCT_RNGState RNG;

}STRUCT_LM;

/*--------------------------------------------------------------------------------*/

extern int INVERSION_MULTI(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
    STRUCT_PARA *Para, STRUCT_LM *LM, STRUCT_SUBSET *Subset);

/*--------------------------------------------------------------------------------*/

