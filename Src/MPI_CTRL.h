
#ifndef MPI_CTRL_h
#define MPI_CTRL_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <mpi.h>

/*--------------------------------------------------------------------------------*/

extern int Mrank;
extern int Mnprocs;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_MPI{
    int rank, nprocs;
}STRUCT_MPI;

/*--------------------------------------------------------------------------------*/

extern int cpu_check(int *cpu_busy, int ncpu);

extern void CONTROL(void);

extern void ABORTED(void);

/*--------------------------------------------------------------------------------*/

#endif /* MPI_CTRL_h */
