
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <string.h>
#include <stdbool.h>

#include "LOG_ERROR.h"

/*--------------------------------------------------------------------------------*/

typedef enum data_type{
    enum_int, enum_flt, enum_dbl, enum_cplx, enum_char
}DATA_TYPE;

/*--------------------------------------------------------------------------------*/

typedef struct Structure_Vector{
    long il, ih, ni;
    DATA_TYPE type;
    void *data;
    void *ptr;
}STRUCT_VECTOR;

typedef struct Structure_Matrix{
    long il, ih, ni;
    long jl, jh, nj;
    DATA_TYPE type;
    void *data;
    void **ptr;
}STRUCT_MATRIX;

typedef struct Structure_Tensor{
    long il, ih, ni;
    long jl, jh, nj;
    long kl, kh, nk;
    DATA_TYPE type;
    void *data;
    void ***ptr;
}STRUCT_TENSOR;

/*--------------------------------------------------------------------------------*/

#define VEC_INT(v)  ((int *)((v).prt))
#define VEC_FLT(v)  ((float *)((v).ptr))
#define VEC_DBL(v)  ((double *)((v).ptr))
#define VEC_CPLX(v)  ((complex double *)((v).ptr))
#define VEC_CHAR(v)  ((char *)((v).ptr))

#define MAT_INT(m)  ((int **)((m).ptr))
#define MAT_FLT(m)  ((float **)((m).ptr))
#define MAT_DBL(m)  ((double **)((m).ptr))
#define MAT_CPLX(m)  ((complex double **)((m).ptr))
#define MAT_CHAR(m)  ((char **)((m).ptr))

#define TENSOR_INT(t)  ((int ***)((t).ptr))
#define TENSOR_FLT(t)  ((float ***)((t).ptr))
#define TENSOR_DBL(t)  ((double ***)((t).ptr))
#define TENSOR_CPLX(t)  ((complex double ***)((t).ptr))

/*--------------------------------------------------------------------------------*/

#define FREE_VECTOR(v)                                    \
    do{                                                   \
      if((v).data) {free((v).data); (v).data = NULL;}     \
      if((v).ptr) {(v).ptr = NULL;}                       \
    }while(0)
   

#define FREE_MATRIX(m)                                    \
    do{                                                   \
      if((m).data) {free((m).data); (m).data = NULL;}     \
      if((m).ptr) {                                       \
        free((m).ptr+(m).il); (m).ptr = NULL;             \
      }                                                   \
    }while(0)


#define FREE_TENSOR(t)                                    \
    do {                                                  \
      if((t).data) {free((t).data);(t).data = NULL;}      \
      if((t).ptr) {                                       \
        free((t).ptr[(t).il]+(t).jl);                     \
        free((t).ptr+(t).il);                             \
        (t).ptr = NULL;                                   \
      }                                                   \
    }while(0)


/*--------------------------------------------------------------------------------*/

extern STRUCT_VECTOR VECTOR(long il, long ih, DATA_TYPE type, bool Init);

extern STRUCT_MATRIX MATRIX(long il, long ih, long jl, long jh, \
    DATA_TYPE type, bool Init);

extern STRUCT_MATRIX MATRIX_TRI(long ih, DATA_TYPE type, bool Init);

extern STRUCT_MATRIX MATRIX_RHO(long ih, DATA_TYPE type, bool Init);

extern STRUCT_TENSOR TENSOR(long il, long ih, long jl, long jh, long kl, \
    long kh, DATA_TYPE type, bool Init);

extern STRUCT_TENSOR TENSOR_TRI(long ih, long jh, DATA_TYPE type, bool Init);

extern STRUCT_TENSOR TENSOR_RHO(long ih,long jh, DATA_TYPE type,bool Init);

/*--------------------------------------------------------------------------------*/

