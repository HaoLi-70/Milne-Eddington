
#ifndef SVD_h
#define SVD_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include "ALLOCATION.h"
#include "PARAMETER.h"

/*--------------------------------------------------------------------------------*/

extern int svdcmp(double **a, int m, int n, double w[], double **v);

extern int svbksb(double **u, double w[], double **v, int m, int n, \
        double b[], double x[]);

extern int SVD_solve(double **A, double *B, double *X, int N, \
        double Threshold);

/*--------------------------------------------------------------------------------*/

#endif /* SVD_h */
