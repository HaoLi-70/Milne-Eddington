
#include "INV_INIT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        18 Apr. 2026  (Hao Li)
          --- Bugfix:  
              A variable for exit-condition checking is added to the struct 
              LM. 

        16 Apr. 2026  (Hao Li)
          --- Updates:  
              rename the keyword fastmode t cache_prof.

        09 Mar. 2026
          --- Initial commit (Hao Li)
              moved some functions here. 
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/


double Gfactor(double J,double L,double S){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Computes the Lande factor.
      Modified:
        25 Apr. 2018
      Input parameters:
        J, The total angular momentum of the electronic cloud.
        L, The total orbital angular momentum of the electronic cloud.
        S, The total spin of the electronic cloud.
      Return:
        Gfactor, The Lande factor.
      Method:
        L-S coupling asumption, and if the number is 0, return 0 for 
        convernience.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
    
    if(J==0) return 0;
        
    return 1.0+(J*(J+1.0)-L*(L+1.0)+S*(S+1.0))/(2.0*J*(J+1.0));
  
}

/*--------------------------------------------------------------------------------*/

double Geffect(double Gu, double Gl, double Ju, double Jl){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Computes the effective Lande factor.
      Modified:
        27 Nov. 2024
      Input parameters:
        Gu, The lande factor of upper level.
        Gl, The lande factor of lower level.
        Ju, The total angular momentum of upper level.
        Jl, The total angular momentum of lower level.
      Return:
        Geffect, The effect Lande factor.
      References:
        LL04 Chapter 3, Equation 3.44.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
         
    return 0.5*(Gu+Gl)+0.25*(Gu-Gl)*(Ju*(Ju+1)-Jl*(Jl+1));
  
}

/*--------------------------------------------------------------------------------*/

int INIT_INV(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_LM *LM, \
    STRUCT_PARA *Para, STRUCT_MPI *Mpi, STRUCT_SUBSET *Subset){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        initialize the inversion.
      Record of revisions:
        6 Mar. 2026. 
      Input parameters:
         Mpi, a structure storing Mpi info.
      Return:
        return the current conts.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    #define INIT_VERBOSE(lvl, fmt, ...)                                  \
      do{                                                                \
        if(Input->verboselv>=lvl){                                       \
          snprintf(MeSS, sizeof(MeSS), fmt, ##__VA_ARGS__);              \
          LOG_WRITE(MeSS, true, true);                                   \
      }                                                                  \
    }while(0)   

    #define L_c 299792458.0

    Para->nline = Input->nline;
    Para->lines = (STRUCT_MELINE *)malloc(Para->nline*sizeof(STRUCT_MELINE));

    /*
    double gu = Gfactor(0.0, 2.0, 2.0);
    double gl = Gfactor(1.0, 1.0, 2.0);
    Lines->Geffect = Geffect(gu, gl, 0.0, 1.0);
    */
    // Lines->Lambda in A, Par->Par[4] in mA, Lines->BShift in mA
    for(int il=0; il<Para->nline; il++){
      Para->lines[il].Lambda0 = Input->Lambda0[il];
      Para->lines[il].Geff = Input->Geff[il];
      Para->lines[il].BShift = 4.6686e-10*Input->Geff[il]\
          *Input->Lambda0[il]*Input->Lambda0[il];
    }
    free(Input->Lambda0);
    Input->Lambda0 = NULL;
    free(Input->Geff);
    Input->Geff = NULL;

    Para->step = Input->step;
    Para->inv = Input->inv;
    Para->regl = Input->regl;
    Para->value_const = Input->value_const;
    Para->Regul_weight = Input->Regul_weight;

    Para->npar = 0;
    for(int ipar=0; ipar<9; ipar++){
      if(Para->inv[ipar]) Para->npar++;
      Para->Limits[ipar][0] = Input->Limits[ipar][0];
      Para->Limits[ipar][1] = Input->Limits[ipar][1];
    }

    Para->VCoeffi = Input->VCoeffi;
    Para->LCoeffi = Input->LCoeffi;

    Para->delta_v = (Stk->Lambda[Stk->nw-1]-Stk->Lambda[0])/Stk->nw \
        /Para->lines->Lambda0*L_c/1e3;

    LM->Damp = Input->Damp;
    LM->Lam_accept = Input->Lam_accept;
    LM->Lam_reject = Input->Lam_reject;
    LM->Lam_Lim[0] = Input->Lam_Lim[0];
    LM->Lam_Lim[1] = Input->Lam_Lim[1];

    LM->Criteria = Input->Chisq_Criteria;
    LM->Convg_Criteria = Input->Convg_Criteria;

    LM->Threshold = Input->Threshold;
    LM->Regul_flag = Input->Regul;
    LM->niter = Input->niter;
    LM->verboselv = Input->verboselv;
    LM->nruns = Input->nruns;
    LM->HMI_REF = Input->HMI_REF;

    Stk->Icriteria = Input->Icriteria;
    Stk->Fadd.nfigures = Input->nfigures;
    Faddeeva_init(&(Stk->Fadd));

    RNG_INIT(&(LM->RNG), Mpi->rank, Mpi->thread_id);

    Input->cache_header.nx = Input->sol_box[0][1] \
        -Input->sol_box[0][0]+1;
    Input->cache_header.ny = Input->sol_box[1][1] \
        -Input->sol_box[1][0]+1;
    Input->counts = Input->cache_header.nx*Input->cache_header.ny;

    size_t nbuf;

    if(Input->nProf<1) Input->nProf=1;
    if(Input->nProf>Input->counts) Input->nProf=Input->counts;
    Subset->nProf = Input->nProf; 
    if(Input->cache_prof){
      Input->nProf = Input->counts;
      if(Mpi->rank == 0){
        nbuf = Input->counts;
        if(Mpi->size == 1) Subset->nProf = Input->counts;
      }else{
        nbuf = Subset->nProf;
      }
    }else{
      nbuf = Subset->nProf;
    }

    Input->cache_header.ncache = Input->cache_header.nx \
        *Input->cache_header.ny/Subset->nProf+1;

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    INIT_VERBOSE(2, " rank %d, tot pixel number: %d, %d, " \
        "prof number in each inv: %d, buffer size: %ld.", Mpi->rank, \
        Input->counts, Input->nProf, Subset->nProf, nbuf);

    Input->profbuff = (double *)malloc(nbuf*Stk->nw*4*sizeof(double));
    Input->resbuff = (double *)malloc(nbuf*10*sizeof(double));

    if(Input->output_fit){
      Input->fitbuff = (double *)malloc(nbuf*Stk->nw*4*sizeof(double));
    }

    Subset->pcounts = 0;
    Subset->coord[0] = Input->sol_box[0][0];
    Subset->coord[1] = Input->sol_box[1][0];
    Input->Subset_flg = (Input->nx == Input->cache_header.nx);

    if(Mpi->size>1 && Mpi->rank==0){ 
      free(Stk->Lambda);
      Stk->Lambda = NULL;
      return 0;
    }

    //Stk->prof = (double *)malloc(4*Stk->nw*sizeof(double));
    Stk->fit_best = (double *)malloc(4*Stk->nw*sizeof(double));
    Stk->fit = (double *)malloc(4*Stk->nw*sizeof(double));
    Stk->syn = (double *)malloc(4*Stk->nw*sizeof(double));
    Para->Par_Best = (double *)malloc(9*sizeof(double));
    Para->Par = (double *)malloc(9*sizeof(double));
    Para->Par_tmp = (double *)malloc(9*sizeof(double));

    LM->Hessian = MATRIX(0, (Para->npar-1), 0, (Para->npar-1),     \
        enum_dbl, false);
    LM->Hessian_new = MATRIX(0, (Para->npar-1), 0, (Para->npar-1), \
        enum_dbl, false);
    LM->Regul_H = MATRIX(0, (Para->npar-1), 0, (Para->npar-1),     \
        enum_dbl, false);
    LM->V = MATRIX(0, (Para->npar-1), 0, (Para->npar-1),           \
        enum_dbl, false);
    LM->Regul_J = (double *)malloc(Para->npar*sizeof(double));
    LM->Sol = (double *)malloc(Para->npar*sizeof(double));
    LM->indx = (int *)malloc(Para->npar*sizeof(int));
    LM->W = (double *)malloc(Para->npar*sizeof(double));
    LM->Jacfvec = (double *)malloc(Para->npar*sizeof(double));

    for(int istk=0; istk<4; istk++){
      Stk->Weights[istk] = Input->Weights[istk];
      Stk->Weights_SQR[istk] = Stk->Weights[istk]*Stk->Weights[istk];
    }

    Stk->ivnoise = (double *)malloc(4*Stk->nw*sizeof(double));
    Stk->Jacobian = (double *)malloc(4*9*Stk->nw*sizeof(double));

    return 0;
}

/*--------------------------------------------------------------------------------*/

