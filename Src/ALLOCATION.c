
#include "ALLOCATION.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        06 Mar. 2026  (Hao Li)
          --- Updates:  
              Redesign the Vector/Matrix/Tensor Allocation Wrapper. 
              Using macros for element access and memory deallocation.
  

        26 Nov. 2024.
          --- Update: new subroutines TENSOR_FLT and FREE_TENSOR_FLT for
              floats (Hao Li)

        30 Oct. 2022.
          --- Initial commit (Hao Li)

    ######################################################################*/

/*--------------------------------------------------------------------------------*/


static inline size_t DTYPE_SIZE(DATA_TYPE t){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        return the data type size.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        t, data type.
      Return:
        return the sizes.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    switch(t){
      case enum_int:  
        return sizeof(int);
      case enum_flt:  
        return sizeof(float);
      case enum_dbl:  
        return sizeof(double);
      case enum_cplx: 
        return sizeof(complex double);
      case enum_char: 
        return sizeof(char);
      default:
        LOG_ERROR(ERR_LVL_ERROR,"DTYPE_SIZE","unknown type");
        return 0;
    }
}

/*--------------------------------------------------------------------------------*/


STRUCT_VECTOR VECTOR(long il, long ih, DATA_TYPE type, bool Init){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a vector v[i] with subscript range nl<=i<=ih.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        il, ih, the subscript range.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_VECTOR v;
    v.il = il; v.ih = ih;
    v.type = type;

    long n = ih-il+1;
    v.ni = n;
    size_t es = DTYPE_SIZE(type);

    if(Init){
      v.data = calloc(n, es);
    }else{
      v.data = malloc(n*es);
    }

    if(!v.data){
      LOG_ERROR(ERR_LVL_ERROR,"VECTOR","allocation failure");
    }

    v.ptr = (char *)(v.data)-il*es;

    return v;
}


/*--------------------------------------------------------------------------------*/

STRUCT_MATRIX MATRIX(long il, long ih, long jl, long jh, \
    DATA_TYPE type, bool Init){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a matrix m[i][j] with subscript range il<=i<=ih, \
            jl<=j<=jh.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        il, ih, jl, jh, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_MATRIX m;
    m.il = il; m.ih = ih;
    m.jl = jl; m.jh = jh;
    m.type = type;

    long nrows = ih-il+1;
    long ncols = jh-jl+1;
    m.ni = nrows;
    m.nj = ncols;

    size_t es = DTYPE_SIZE(type);

    m.ptr = malloc(nrows*sizeof(void *));
    if(!m.ptr){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX", "Error: allocation failure 1\n");
    }
    m.ptr -= il;

    if(Init){
      m.ptr[il] = calloc(nrows*ncols,es);
    }else{
      m.ptr[il] = malloc(nrows*ncols*es);
    }

    if(!m.ptr[il]){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX", "Error: allocation failure 2\n");
    }

    m.data = m.ptr[il];
    m.ptr[il] = (char *)m.data-jl*es;

    for(long ii=il+1; ii<=ih; ii++){
      m.ptr[ii] = (char *)m.ptr[ii-1]+ncols*es;
    }

    return m;
}


/*--------------------------------------------------------------------------------*/

STRUCT_MATRIX MATRIX_TRI(long ih, DATA_TYPE type, bool Init){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a matrix m[i][j] with subscript range 0<=i<=ih, \
            0<=j<=i.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        ih, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_MATRIX m;
    m.il = 0; m.ih = ih;
    m.type = type;

    long nrow = ih+1;
    long total = nrow*(nrow+1)/2; 

    size_t es = DTYPE_SIZE(type);

    m.ptr = malloc(nrow*sizeof(void *));
    if(!m.ptr){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_TRI", "Error: allocation failure 1\n");
    }

    if(Init){
      m.ptr[0] = calloc(total,es);
    }else{
      m.ptr[0] = malloc(total*es);
    }

    if(!m.ptr[0]){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_TRI", "Error: allocation failure 2\n");
    }
  
    m.data = m.ptr[0];

    for(long ii = 1; ii<=ih; ii++){
      m.ptr[ii] = (char *)m.ptr[ii-1]+ii*es;
    }

    return m;
}

/*--------------------------------------------------------------------------------*/

STRUCT_MATRIX MATRIX_RHO(long ih, DATA_TYPE type, bool Init){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a matrix m[i][j] with subscript range 0<=i<=ih, \
            -i<=j<=i.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        ih, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_MATRIX m;
    m.il = 0; m.ih = ih;
    m.type = type;

    long nrow = ih+1;
    long total = nrow*nrow;
    size_t es = DTYPE_SIZE(type);

    m.ptr = malloc(nrow*sizeof(void *));
    if(!m.ptr){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_RHO", "Error: allocation failure 1\n");
    }

    if(Init){
      m.ptr[0] = calloc(total, es);
    }else{
      m.ptr[0] = malloc(total*es);
    }

    if(!m.ptr[0]){ 
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_RHO", "Error: allocation failure 2\n");
    }

    m.data = m.ptr[0];

    for(long ii=1; ii<=ih; ii++){
      m.ptr[ii] = (char *)m.ptr[ii-1]+2*ii*es;
    }

    return m;
}


/*--------------------------------------------------------------------------------*/

STRUCT_TENSOR TENSOR(long il, long ih, long jl, long jh, long kl, \
    long kh, DATA_TYPE type, bool Init){
   
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a tensor t[i][j][k] with subscript range
            il<=i<=ih, jl<=j<=jh, kl<=k<=kh.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        il, ih, jl, jh, kl, kh, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_TENSOR t;
    t.il = il; t.ih = ih;
    t.jl = jl; t.jh = jh;
    t.kl = kl; t.kh = kh;
    t.type = type;

    long in = ih-il+1;
    long jn = jh-jl+1;
    long kn = kh-kl+1;
    
    t.ni = in;
    t.nj = jn;
    t.nk = kn;

    size_t es = DTYPE_SIZE(type);

    t.ptr = malloc(in*sizeof(void **));
    if(!t.ptr){ 
      LOG_ERROR(ERR_LVL_ERROR, "TENSOR", "Error: allocation failure 1\n");
    }
    t.ptr -= il;

    t.ptr[il] = malloc(in*jn*sizeof(void *));
    if(!t.ptr[il]){ 
      LOG_ERROR(ERR_LVL_ERROR, "TENSOR", "Error: allocation failure 2\n");
    }
    t.ptr[il] -= jl;
    for(long ii=il+1; ii<=ih; ii++){
      t.ptr[ii] = t.ptr[ii-1]+jn;
    }

    if(Init){
      t.ptr[il][jl] = calloc(in*jn*kn, es);
    }else{
      t.ptr[il][jl] = malloc(in*jn*kn*es);
    }

    if(!t.ptr[il][jl]){ 
      LOG_ERROR(ERR_LVL_ERROR, "TENSOR", "Error: allocation failure 3\n");
    }

    t.data = t.ptr[il][jl];
    t.ptr[il][jl] = (char*)t.data-kl*es;

    for(long ii=il+1; ii<=ih; ii++){
      t.ptr[ii][jl] = (char*)t.ptr[ii-1][jl]+jn*kn*es;
    }

    for(long ii=il; ii<=ih; ii++){
      for(long jj =jl+1; jj<=jh; jj++){
        t.ptr[ii][jj] = (char*)t.ptr[ii][jj-1]+kn*es;
      }
    }

    return t;
}

/*--------------------------------------------------------------------------------*/

STRUCT_TENSOR TENSOR_TRI(long ih, long jh, DATA_TYPE type, bool Init){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a tensor t[i][j][k] with subscript range
            0<=i<=ih, 0<=j<=jh, 0<=k<=j.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        il, ih, jl, jh, kl, kh, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_TENSOR t;
    t.il = t.jl = t.kl = 0;                         
    t.ih = ih;
    t.jh = t.kh = jh;
    t.type = type;

    long nx = ih+1;
    long ny = jh+1;
    t.ni = nx;
    t.nj = ny;

    size_t es = DTYPE_SIZE(type);

    t.ptr = malloc(nx*sizeof(void **));
    if(!t.ptr){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_TRI", \
          "Error: allocation failure 1");
    }

    t.ptr[0] = malloc(nx*ny*sizeof(void *));
    if(!t.ptr[0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_TRI", \
          "Error: allocation failure 2");
    }

    for(long ii=1; ii<nx; ii++){
      t.ptr[ii] = t.ptr[ii-1]+ny;
    }

    long tri = ny*(ny+1)/2;
    long total = nx*tri;

    if(Init){
      t.ptr[0][0] = calloc(total, es);
    }else{
      t.ptr[0][0] = malloc(total*es);
    }

    if(!t.ptr[0][0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_TRI", \
          "Error: allocation failure 3");
    }
    
    t.data = t.ptr[0][0];

    for(long ii=1; ii<nx; ii++){
      t.ptr[ii][0] = (char *)t.ptr[ii-1][0]+tri*es;
    }

    for(long ii=0; ii<nx; ii++){
      for(long jj=1; jj<ny; jj++){
        t.ptr[ii][jj] = (char*)t.ptr[ii][jj-1]+jj*es;
      }
    }

    return t;
}


/*--------------------------------------------------------------------------------*/

STRUCT_TENSOR TENSOR_RHO(long ih,long jh, DATA_TYPE type,bool Init){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        allocate a tensor t[i][j][k] with subscript range
            0<=i<=ih, 0<=j<=jh, -j<=k<=j.
      Record of revisions:
        6 Mar. 2026
      Input parameters:
        il, ih, jl, jh, kl, kh, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the structure.
      Reference:
        Numerical Recipes in C 2nd edtion.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
    STRUCT_TENSOR t;

    t.il = t.jl = t.kl = 0;                         
    t.ih = ih;
    t.jh = t.kh = jh;
    t.type = type;

    long nx = ih+1;
    long ny = jh+1;

    size_t es = DTYPE_SIZE(type);

    t.ptr = malloc(nx*sizeof(void **));
    if(!t.ptr){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_RHO", \
          "Error: allocation failure 1");
    }

    t.ptr[0] = malloc(nx*ny*sizeof(void *));
    if(!t.ptr[0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_RHO", \
          "Error: allocation failure 2");
    }

    for(long ii=1; ii<nx; ii++) t.ptr[ii] = t.ptr[ii-1]+ny;

    long tri = ny*ny;
    long total = nx*tri;

    if(Init){
      t.ptr[0][0] = calloc(total, es);
    }else{
      t.ptr[0][0] = malloc(total*es);
    }

    if(!t.ptr[0][0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_RHO", \
          "Error: allocation failure 3");
    }

    t.data = t.ptr[0][0];

    for(long ii=1; ii<nx; ii++){
      t.ptr[ii][0] = (char*)t.ptr[ii-1][0]+tri*es;
    }

    for(long ii=0; ii<nx; ii++){
      for(long jj=1; jj<ny; jj++){
        t.ptr[ii][jj] = (char*)t.ptr[ii][jj-1]+(2*jj)*es;
      }
    }

    return t;
}

/*--------------------------------------------------------------------------------*/
