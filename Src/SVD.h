
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include "ALLOCATION.h"
#include "LOG_ERROR.h"

#ifdef USE_LAPACKE
#include <lapacke.h>
#endif

/*--------------------------------------------------------------------------------*/

extern int svdcmp_float(STRUCT_MATRIX *am, float *w, STRUCT_MATRIX *vm);

extern int svdcmp_double(STRUCT_MATRIX *am, double *w, STRUCT_MATRIX *vm);

extern int svbksb_float(STRUCT_MATRIX *am, STRUCT_MATRIX *vm, \
    const float *w, const float *b, float *x);

extern int svbksb_double(STRUCT_MATRIX *am, STRUCT_MATRIX *vm, \
    const double *w, const double *b, double *x);

extern int svd_dbl(STRUCT_MATRIX *am, double *w, STRUCT_MATRIX *vm);


//extern int SVD_solve(double **A, double *B, double *X, int N, 
//        double Threshold);

/*--------------------------------------------------------------------------------*/

