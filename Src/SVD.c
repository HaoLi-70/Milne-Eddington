
#include "SVD.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        06 Mar. 2026  (Hao Li)
          --- Updates:  
              All the functions are moved to SVD_TEMPLATE.h.
              The functions will be generated from the templates.

        28 Jun. 2024
          --- Initial commit (Hao Li)
            
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

#define TYPE double
#define SVD_FUNC svdcmp_double
#define SVB_FUNC svbksb_double
#define PYTHAG_FUNC pythag_double
#include "SVD_TEMPLATE.h"
#undef TYPE
#undef SVD_FUNC
#undef SVB_FUNC
#undef PYTHAG_FUNC

#define TYPE float
#define SVD_FUNC svdcmp_float
#define SVB_FUNC svbksb_float
#define PYTHAG_FUNC pythag_float
#include "SVD_TEMPLATE.h"
#undef TYPE
#undef SVD_FUNC
#undef SVB_FUNC
#undef PYTHAG_FUNC

/*--------------------------------------------------------------------------------*/

int svd_dbl(STRUCT_MATRIX *am, double *w, STRUCT_MATRIX *vm){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        SVD using either LAPACKE or the Numerical Recipes method.
      Record of revisions:
        25 Mar. 2026
      Input parameters:        
        am, a structure storing the matrix A[0...m-1][0...n-1] and the 
          dimensions m, and n.
      Output parameters:        
        am, A matrix is replaced by U matrix.
        w, the diagnonal matrix.
        vm, a structure storing the matrix v[0...n-1][0...n-1].
    ######################################################################*/
/*--------------------------------------------------------------------------------*/


#ifdef USE_LAPACKE

    double *A = (double *)am->data;
    double *V = (double *)vm->data;
    const int ndim = am->ni;
    double *U = (double *)malloc(ndim*ndim*sizeof(double));
    double *VT = (double *)malloc(ndim*ndim*sizeof(double));

    int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'A', ndim, ndim, \
        A, ndim, w, U, ndim, VT, ndim);

    if(info>0) LOG_ERROR(ERR_LVL_ERROR, "svdcmp_flt", \
        "error in LAPACKE_dgesdd\n");
     
    memcpy(A, U, ndim*ndim*sizeof(double));

    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
           V[i*ndim+j] = VT[j*ndim+i]; 

    free(U);
    free(VT);
      
#else
    svdcmp_double(am, w, vm);
#endif

    return 0;
}

/*--------------------------------------------------------------------------------*/
