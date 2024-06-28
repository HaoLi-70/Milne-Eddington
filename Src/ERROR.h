
#ifndef ERROR_h
#define ERROR_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "MPI_CTRL.h"
#include "PARAMETER.h"

/*--------------------------------------------------------------------------------*/

enum error_level {enum_warning, enum_error};

#define MAX_MESSAGE_LENGTH 2000

extern char MeSS[MAX_MESSAGE_LENGTH];

/*--------------------------------------------------------------------------------*/

extern int Error(enum error_level error_lv, const char *rname,\
        const char *Emes, char *Filename);

extern int VerboseM(const char *Emes, char *Filename, bool screen);

extern int Verbose(const char *Emes, char *Filename, bool cut);

extern int Verbose_model(double *Par, char *Filename, bool cut);

extern bool FILE_EXIST(char *Filename);

/*--------------------------------------------------------------------------------*/

#endif /* ERROR_h */
