
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

/*--------------------------------------------------------------------------------*/

#define MAX_MESSAGE_LENGTH 2000
#define MAX_BUFRER_SIZE 2500

extern char MeSS[MAX_MESSAGE_LENGTH];

typedef enum ERR_LVL {ERR_LVL_WARNING, ERR_LVL_ERROR} ERR_LVL;

/*--------------------------------------------------------------------------------*/

extern void ABORTED(void);

extern void LOG_INIT(const char *filename);

extern void LOG_FINALIZE(void);

extern void LOG_WRITE(const char *msg, bool to_screen, bool verbose_flag);

extern void LOG_ERROR(ERR_LVL lvl, const char *routine, const char *msg);

extern void LOG_MODEL(const double *Par, bool verbose_flag);

extern bool FILE_EXIST(const char *Filename);

/*--------------------------------------------------------------------------------*/


