
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_OPEMP
#include <omp.h>
#endif

/*--------------------------------------------------------------------------------*/

typedef struct Struct_MPI{

    int rank, size, thread_id, nthreads;
    long *idum;

#ifdef USE_MPI
    MPI_Comm comm;
#endif

}STRUCT_MPI;

/*--------------------------------------------------------------------------------*/

extern int cpu_check(int *cpu_busy, int ncpu);

extern void MPI_SETUP(STRUCT_MPI *Mpi);

extern void OPENMP_SETUP(STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/

