
#pragma once

/*--------------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Faddeeva{

    // number of figures and the corresponding index, max number of iteration.
    int nfigures, indxfig, nmax;

    // the smallest value, log(Rmin0), sqrt(log(Rmin0)), a, and a*a;
    double Rmin0, logRmin0, sqrt_logRmin0, a, AA;
    
    // precomputed exp(-a^2n^2)
    double *Expa2n2;

}STRUCT_FADDEEVA;

/*--------------------------------------------------------------------------------*/

int Faddeeva916(double Nu, double y, double *H, double *L, \
    STRUCT_FADDEEVA *Fadd);

extern int Faddeeva_init(STRUCT_FADDEEVA *Fadd);

extern int Faddeeva(double Nu, double y, double *H, double *L, \
    STRUCT_FADDEEVA *Fadd);

/*--------------------------------------------------------------------------------*/
