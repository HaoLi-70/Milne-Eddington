
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_OPEMP
#include <omp.h>
#endif

#ifdef USE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#include "MPI_INIT.h"

/*--------------------------------------------------------------------------------*/

#define NTAB 32

typedef struct Struct_RNG_State{
#ifdef USE_GSL
    gsl_rng *gsl;
#else
    int32_t idum;
    int32_t iy;
    int32_t iv[NTAB];

    bool iset;
    double gset; 
#endif
}STRUCT_RNGState;

/*--------------------------------------------------------------------------------*/

extern void RNG_INIT(STRUCT_RNGState *State, int rank, int thread_id);

extern double RNG_UNIFORM(STRUCT_RNGState *State);

extern double RNG_GAUSS(STRUCT_RNGState *State);

/*--------------------------------------------------------------------------------*/
