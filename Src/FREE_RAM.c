
#include "FREE_RAM.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:

        29 Mar. 2026  (Hao Li)
          --- Initial commit
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

void FREERAM(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_LM *LM, \
    STRUCT_PARA *Para, STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/    
    /*######################################################################
      Purpose:
        free ram.
      Record of revisions:
        29 Mar. 2026.
      Input parameters:
        Input, input configuration.
        Stk, a structure storing  the Stokes profiles.
        LM, a structure storing  the hessian matrix.
        Para, a structure storing  the model parameters.
        Mpi, a structure storing mpi info.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    #define FREE_PRT(a)                             \
      do{                                           \
        if(a){                                      \
          free(a);                                  \
          a = NULL;                                 \
        }                                           \
      }while(0)

    if(Mpi->rank == 0) FREE_PRT(Input->cache);

    FREE_PRT(Stk->Fadd.Expa2n2);
    if(Input->output_fit){
      FREE_PRT(Input->fitbuff);
    }
    FREE_PRT(Para->lines);
    FREE_PRT(Input->profbuff);
    FREE_PRT(Input->resbuff);
    if(Mpi->rank==0){
      LOG_FINALIZE();
    }else if(Input->verboselv>0){
      LOG_FINALIZE();
    }
    if(Mpi->size>1 && Mpi->rank==0) return;
    
    FREE_PRT(Stk->Lambda);
    FREE_PRT(Stk->fit_best);
    FREE_PRT(Stk->fit);
    FREE_PRT(Stk->syn);
    FREE_PRT(Para->Par_Best);
    FREE_PRT(Para->Par);
    FREE_PRT(Para->Par_tmp);
    FREE_PRT(LM->Regul_J);
    FREE_PRT(LM->Sol);
    FREE_PRT(LM->indx);
    FREE_PRT(LM->W);
    FREE_PRT(LM->Jacfvec);
    FREE_PRT(Stk->ivnoise);
    FREE_PRT(Stk->Jacobian);
    
    FREE_MATRIX(LM->Hessian);
    FREE_MATRIX(LM->Hessian_new);
    FREE_MATRIX(LM->Regul_H);
    FREE_MATRIX(LM->V);

    return;
}

/*--------------------------------------------------------------------------------*/
