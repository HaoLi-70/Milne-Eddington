
#include "LM_FIT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:

        06 Mar. 2026  (Hao Li)
          --- Updates:  
              slightly improve the performanceß. 

        18 Apr. 2025
          --- Updates: More runs if chisq twice larger than the criteria
                       The random jump of the model parameter is changed 
                       (Hao Li)

        11 Apr. 2025
          --- bugfix:  does not save the best fit (Hao Li)

        27 Nov. 2024
          --- Updates:  the random seeds are moved to the Mpi structure 
                        (Hao Li)

        28 Jun. 2024
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

// Pi Ratio of circumference to diameter
static const double L_Pi = 3.14159265358979323846;

#define M_SWAP(a,b)                                 \
    do{                                             \
      (void) (&a==&b);                              \
      __typeof__(a) tmp=(a);                        \
      (a)=(b);                                      \
      (b)=tmp;                                      \
    }while(0)

#define converg_check(Chisq,Chisq_new)              \
    (((Chisq-Chisq_new)/Chisq)<LM->Criteria)

/*--------------------------------------------------------------------------------*/

int sav_par(STRUCT_PARA *Para, STRUCT_STK *Stk){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Save the best fit parameters.
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        Para, a structure storing parameters.
        Stk, a structure storing the Stokes profiles.
      Output parameter:
        Para, a structure storing parameters.
        Stk, a structure storing the Stokes profiles.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    M_SWAP(Para->Par, Para->Par_Best);
    M_SWAP(Stk->fit, Stk->fit_best);

    Para->Chisq_Best = Para->Chisq;

    return  0;
}

/*--------------------------------------------------------------------------------*/

static inline double Loss_Function(STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Compute the merit function.
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        Input, Input configuration.
        Stk, a structure storing Stokes profiles.
      Return:
        return the chisq.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    double Chisq = 0, diff;
    int nw = Stk->nw;

    const double *syn = Stk->syn;
    const double *prof = Stk->prof;
    const double *ivnoi = Stk->ivnoise;

    for(int iw=0; iw<nw*4; iw++){
      diff = syn[iw]-prof[iw];
      Chisq += diff*diff*ivnoi[iw];
    }
       
    Chisq /= (4.*nw);

    return Chisq;
}

/*--------------------------------------------------------------------------------*/

static inline bool Lambda_propose(STRUCT_LM *LM, bool accepted){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Propose a lambda factor.
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        LM, a structure storing the hessian matrix.
        accepted, accepted or not.
      Output parameters:
        LM, the damping factor.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(accepted){
      LM->Damp /= LM->Lam_accept;
      if(LM->Damp<LM->Lam_Lim[0]) LM->Damp = LM->Lam_Lim[0];
      return true;
    }else{
      
      LM->Damp *= LM->Lam_reject;
      if(LM->Damp>LM->Lam_Lim[1]){
        return false;
      }
    }

    return true;
}

/*--------------------------------------------------------------------------------*/

static int Random_Jump(STRUCT_PARA *Para, STRUCT_LM *LM, int irun){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        jumping the initial azimuth another value.
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        Para, a structure storing the model parameters.
        LM, a structure storing the seeds for RNG.
        irun, index;
      Output parameter:
        Para, a structure storing the model parameters.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    STRUCT_RNGState *State = &(LM->RNG);
    double *modelpara = Para->Par;
    double *parabest = Para->Par_Best;
    double *paraguess = Para->Par_Guess;

    modelpara[1] = paraguess[1]+L_Pi/8.*RNG_GAUSS(State);
    modelpara[2] = paraguess[2]+L_Pi/4.*(irun+0.3*RNG_GAUSS(State));
    modelpara[3] = paraguess[3]+2.*RNG_GAUSS(State);
    modelpara[5] = parabest[5];
    modelpara[7] = parabest[7];
    
    if(Para->Chisq_Best < LM->Criteria*4){
      modelpara[0] = parabest[0]+300*RNG_GAUSS(State);
      modelpara[4] = parabest[4]+10*RNG_GAUSS(State);
      modelpara[6] = parabest[6]+5.*RNG_GAUSS(State);
      modelpara[8] = parabest[8]+0.2*RNG_GAUSS(State);

    }else{ 
      
      modelpara[8] = paraguess[8]+0.2*RNG_GAUSS(State);

      if(parabest[0]>1500 || paraguess[0]>1500){ 
        modelpara[0] = paraguess[0]+800*RNG_GAUSS(State);
        modelpara[4] = paraguess[4]+10*RNG_GAUSS(State);
        modelpara[6] = paraguess[6]+5*RNG_GAUSS(State);
      
      }else if(paraguess[0]>500){
        modelpara[0] = paraguess[0]+100*RNG_GAUSS(State);
        modelpara[4] = paraguess[4]+10*RNG_GAUSS(State);
        modelpara[6] = paraguess[6]+2*RNG_GAUSS(State);
      }else{
        modelpara[0] = paraguess[0]+20*RNG_GAUSS(State);
        modelpara[4] = paraguess[4]+10*RNG_GAUSS(State);
        modelpara[6] = paraguess[6]+2*RNG_GAUSS(State);
      } 
    }
    
    bounds_check(Para->Par, Para);

    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Par_update(double *Par, double *Sol, STRUCT_PARA *Para){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Update the parameters according to the solution.
      Record of revisions:
        09 Mar. 2026
      Input parameters:
        Para, the original parameters.
        Sol, the solutions.
        Para, a structure storing the model parameters.
      Output parameters:
        Para, the updated parameters.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(Para->npar==9){
      for(int ipar=0; ipar<9; ipar++){
        Par[ipar] += Sol[ipar]*Para->step[ipar];
      }
    }else{
      int ipar = 0;
      for(int ii=0; ii<9; ii++){
        if(!Para->inv[ii]) continue;
        Par[ii] += Sol[ipar]*Para->step[ii];
        ipar++;
      }
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

static inline int Hessian_Compute(STRUCT_STK *Stk, STRUCT_PARA *Para, \
    STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Compute the Hessian matrix.
      Record of revisions:
        09 Mar. 2026
      Input parameters:
        Stk, a structure storing the Stokes profiles.
        Para, a structure storing the model parameters.
        LM, a structure storing the hessian matrix.
      Output parameters:
        LM, a structure storing the hessian matrix.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    const int nw = Stk->nw;
    double *ivnoiI = Stk->ivnoise;
    double *ivnoiQ = ivnoiI+nw;
    double *ivnoiU = ivnoiQ+nw;
    double *ivnoiV = ivnoiU+nw;

    double *synI = Stk->syn;
    double *synQ = synI+nw;
    double *synU = synQ+nw;
    double *synV = synU+nw;

    double *profI = Stk->prof;
    double *profQ = profI+nw;
    double *profU = profQ+nw;
    double *profV = profU+nw;

    double *JacobI = Stk->Jacobian;
    double *JacobQ = JacobI+9*nw;
    double *JacobU = JacobQ+9*nw;
    double *JacobV = JacobU+9*nw;
    
    #define JacbI(ipar,iw) JacobI[ipar*nw+iw]
    #define JacbQ(ipar,iw) JacobQ[ipar*nw+iw]
    #define JacbU(ipar,iw) JacobU[ipar*nw+iw]
    #define JacbV(ipar,iw) JacobV[ipar*nw+iw]

    double *Jvec = LM->Jacfvec;
    const int ndim = LM->Hessian.ni;
    double *Hessian = (double *)(LM->Hessian.data);
    #define Hess(i,j) Hessian[i*ndim+j]
    
    double sum;

    if(Para->npar==9){
      for(int ipar=0; ipar<9; ipar++){
        for(int jpar=ipar; jpar<9; jpar++){
          sum = 0;
          for(int iw=0; iw<nw; iw++){
            sum += ivnoiI[iw]*JacbI(ipar,iw)*JacbI(jpar,iw);
            sum += ivnoiQ[iw]*JacbQ(ipar,iw)*JacbQ(jpar,iw);
            sum += ivnoiU[iw]*JacbU(ipar,iw)*JacbU(jpar,iw);
            sum += ivnoiV[iw]*JacbV(ipar,iw)*JacbV(jpar,iw);
          }
          Hess(jpar,ipar) = sum;
          Hess(ipar,jpar) = sum;

        }

        sum = 0;
        for(int iw=0; iw<nw; iw++){
          sum -= ivnoiI[iw]*JacbI(ipar,iw)*(synI[iw]-profI[iw]);
          sum -= ivnoiQ[iw]*JacbQ(ipar,iw)*(synQ[iw]-profQ[iw]);
          sum -= ivnoiU[iw]*JacbU(ipar,iw)*(synU[iw]-profU[iw]);
          sum -= ivnoiV[iw]*JacbV(ipar,iw)*(synV[iw]-profV[iw]);
        }
        Jvec[ipar] = sum;
      }

    }else{
      int ipar = 0;
      for(int ii=0; ii<9; ii++){
        if(!Para->inv[ii]) continue;
        int jpar = ipar;
        for(int jj=ii; jj<9; jj++){
          if(!Para->inv[jj]) continue;
          sum = 0;
          for(int iw=0; iw<nw; iw++){
            sum += ivnoiI[iw]*JacbI(ii,iw)*JacbI(jj,iw);
            sum += ivnoiQ[iw]*JacbQ(ii,iw)*JacbQ(jj,iw);
            sum += ivnoiU[iw]*JacbU(ii,iw)*JacbU(jj,iw);
            sum += ivnoiV[iw]*JacbV(ii,iw)*JacbV(jj,iw);
          }
          Hess(jpar,ipar) = sum;
          Hess(ipar,jpar) = sum;
          jpar++; 
        }

        sum = 0;
    
        for(int iw=0; iw<nw; iw++){
          sum -= ivnoiI[iw]*JacbI(ii,iw)*(synI[iw]-profI[iw]);
          sum -= ivnoiQ[iw]*JacbQ(ii,iw)*(synQ[iw]-profQ[iw]);
          sum -= ivnoiU[iw]*JacbU(ii,iw)*(synU[iw]-profU[iw]);
          sum -= ivnoiV[iw]*JacbV(ii,iw)*(synV[iw]-profV[iw]);
        }
        
        Jvec[ipar] = sum;
        ipar++;
      }   
    }

    return  0;
}

/*--------------------------------------------------------------------------------*/

static inline int Get_Regul(double *Par, STRUCT_PARA *Para, STRUCT_LM *LM, \
    bool matrix_flag){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Computes the regularization.
      Record of revisions:
        09 Mar. 2026
      Input parameters:        
        Par, the model parameters.
        Para, a structure storing regularization coefficients.
        LM, a structure storing the hessian matrix and J vector.
        matrix_flag, return the matrix or not.
      Output parameters: 
        LM, a structure storing the hessian matrix and J vector.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    int ipar = -1;
    double LL = 1., Resi = 0;
    LM->Penalty = 0.0;

    for(int ii=0; ii<9; ii++){
      if(!Para->inv[ii]) continue;
      ipar++;
      if(!Para->regl[ii]) continue;
      Resi = Par[ii]-Para->value_const[ii];

      LM->Penalty += Resi*Resi*Para->Regul_weight[ii];

      if(matrix_flag){
        double **hessian = (double **)(LM->Hessian.ptr);
        double **regu = (double **)(LM->Regul_H.ptr);
        LM->Regul_J[ipar] = LL*Resi*Para->Regul_weight[ii];
        regu[ipar][ipar] = LL*LL*Para->Regul_weight[ii];
        LM->Jacfvec[ipar] -= LM->Regul_J[ipar]*LM->ratio;
        hessian[ipar][ipar] += regu[ipar][ipar]*LM->ratio;
      }
    }
   
    return 0;
}

/*--------------------------------------------------------------------------------*/

static inline int SVD_solve_mem(STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Solving linear equations (A[][]*X[]=B[]) with svd decomposition.
      Record of revisions:
        09 Mar. 2026
      Input parameters:        
        LM, a structure storing the hessian matrix and J vector.
      Output parameters:
        LM, a structure storing the solutions.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/


    const double maxstep = 0.5;
    const int ndim = LM->Hessian.ni;
    double threshold = LM->Threshold;
    double Wthreshold;
    
    svd_dbl(&(LM->Hessian), LM->W, &(LM->V));

    double WMAX = LM->W[0];
#ifndef USE_LAPACKE
    for (int ii=1; ii<ndim; ii++){ 
      WMAX = WMAX > LM->W[ii] ? WMAX : LM->W[ii];
    }
#endif

    for(int ii=0; ii<=ndim; ii++){ 
      if(LM->W[ii]<Wthreshold){
        LM->W[ii] = 0.;
      }
    }

    svbksb_double(&(LM->Hessian), &(LM->V), LM->W, LM->Jacfvec, \
        LM->Sol);

    bool bracked;
    while(1){

      Wthreshold = WMAX*threshold;

      for(int ii=0; ii<=ndim; ii++){ 
        if(LM->W[ii]<Wthreshold){
          LM->W[ii] = 0.;
        }
      }

      svbksb_double(&(LM->Hessian), &(LM->V), LM->W, LM->Jacfvec, \
        LM->Sol);

      bracked = true;
      for(int ii=0; ii<ndim; ii++){ 
        if(LM->Sol[ii]>maxstep || LM->Sol[ii]<-maxstep){ 
          bracked = false;
          break;
        }
      }

      if(bracked) break;

      if(threshold > 2e-3){
        for(int ii=0; ii<ndim; ii++){ 
        if(LM->Sol[ii]>maxstep){ 
          LM->Sol[ii] = maxstep;
        }else if(LM->Sol[ii]<-maxstep){ 
          LM->Sol[ii] = -maxstep;
        }
      }
        break;
      }else{
        threshold *= 5.;
      }
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static inline double Trial_Synthesis(STRUCT_STK *Stk, STRUCT_PARA *Para, \
    STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Trial synthesis
      Record of revisions:
        09 Mar. 2026
      Input parameters:        
        Stk, a structure storing the synthesized Stokes profiles.
        Para, a structure storing the model parameters
        LM, a structure storing the hessian matrix.
      Output parameters:
        Stk, a structure storing the synthesized Stokes profiles.
      return:
        return the chisq.
    ######################################################################*/
/*----------------------------------------------------------------------------*/  

    const int ndim = LM->Hessian.ni;
    double *Hessian = (double *)(LM->Hessian.data);
    double *Hessian_new = (double *)(LM->Hessian_new.data);
    #define Hessnew(i,j) Hessian_new[i*ndim+j]

    for(int ii=0; ii<ndim*ndim; ii++){
      Hessian_new[ii] = Hessian[ii];
    }

    //memcpy((void *)LM->Sol, (void *)LM->Jacfvec, ndim*sizeof(double));

    double damp_factor = 1.0 + LM->Damp;
    for(int ii=0; ii<ndim; ii++){
      Hessnew(ii,ii) *= damp_factor;
    }

    SVD_solve_mem(LM);

    for(int ipar=0; ipar<9; ipar++){
      Para->Par_tmp[ipar] = Para->Par[ipar];
    }

    Par_update(Para->Par_tmp, LM->Sol, Para);
    bounds_check(Para->Par_tmp, Para);

    Milne_Eddington_Single(Para->Par_tmp, Stk, Para, false);

    LM->chisq = Loss_Function(Stk);

    if(LM->Regul_flag){
      Get_Regul(Para->Par_tmp, Para, LM, false);
    }

    return  LM->chisq;
}

/*--------------------------------------------------------------------------------*/

static int LM_FIT(STRUCT_STK *Stk, STRUCT_PARA *Para, STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        the Levenberg Marquardt.
      Record of revisions:
        09 Mar. 2026
      Input parameters:
        Stk, structure with the Stokes profiles
        Par, structure with the model parameters
        LM, structure with the hessian matrix.
      Return:
        return the num of iteration check.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
    

    Milne_Eddington_Single(Para->Par, Stk, Para, true);
    const double limit = 0.1;
    LM->ratio = 1.0;

    double Chisq = Loss_Function(Stk);
    
    if(LM->Regul_flag){
      Get_Regul(Para->Par, Para, LM, false);
      if(LM->Penalty>limit*LM->chisq){
        LM->ratio = limit*Chisq/LM->Penalty;
      }
      sprintf(MeSS, " *** Init chisq = %.6e; Penality = %.6e *** \n", \
          Chisq, LM->Penalty*LM->ratio);

    }else{
      sprintf(MeSS, " *** Init chisq = %.6e; *** \n", Chisq);
    }

    LOG_WRITE(MeSS, true, LM->verboselv>1);
    LOG_MODEL(Para->Par, LM->verboselv>1);

    int it;
    bool accepted = false;
    double Chisq_new;
    double dpfactor = LM->Damp;

    for(it=0; it<LM->niter; it++){
      
      accepted = false;

      sprintf(MeSS, " *** iteration : %d  \n", it);
      LOG_WRITE(MeSS, true, LM->verboselv>2);

      Hessian_Compute(Stk, Para, LM);

      if(LM->Regul_flag){
        Get_Regul(Para->Par, Para, LM, true);
      }

      do{
        Chisq_new = Trial_Synthesis(Stk, Para, LM);

        if (Chisq_new < Chisq){
          if(LM->Regul_flag){
            sprintf(MeSS, "    * trail accepted. chisq = %.6e; " \
                "damping = %.6e; Penalty = %.6e *\n", Chisq_new, \
                LM->Damp, LM->Penalty*LM->ratio);
            if(LM->Penalty>limit*LM->chisq){
              LM->ratio = limit*Chisq/LM->Penalty;
            }

          }else{
            sprintf(MeSS, "    * trail accepted. chisq = %.6e; " \
                "damping = %.6e; *\n", Chisq_new, LM->Damp);
          }
          LOG_WRITE(MeSS, true, LM->verboselv>2);

          accepted = true;

          M_SWAP(Para->Par, Para->Par_tmp);
          M_SWAP(Stk->syn, Stk->fit);

          LOG_MODEL(Para->Par, LM->verboselv>3);
              
          Milne_Eddington_Single(Para->Par, Stk, Para, true);

          Lambda_propose(LM, accepted);

        }else{
          if(LM->Regul_flag){
            sprintf(MeSS, "    * trail chisq = %.6e; damping = %.6e;" \
                " Penalty = %.6e *\n", \
                Chisq_new, LM->Damp, LM->Penalty*LM->ratio);
          }else{
            sprintf(MeSS, "    * trail chisq = %.6e; damping = %.6e" \
                " *\n", Chisq_new, LM->Damp);
          }
          LOG_WRITE(MeSS, true, LM->verboselv>2);

          accepted =false;
          if(!Lambda_propose(LM, accepted)) break;
        }

      }while(!accepted);
      LM->Damp = dpfactor;
      sprintf(MeSS, "    * iteration : %d  chisq : %.6e -> %.6e *\n", it, \
          Chisq, Chisq_new);
      LOG_WRITE(MeSS, true, LM->verboselv>2);

      if(accepted){
        if(converg_check(Chisq, Chisq_new)){
          Chisq = Chisq_new;
          break;
        } 
        Chisq = Chisq_new;
      }else{
        break;
      }
    }

    if(LM->Regul_flag){
      sprintf(MeSS, " --- final chisq = %.6e; Penalty = %.6e --- \n", \
          Chisq, LM->Penalty*LM->ratio);
    }else{
      sprintf(MeSS, " --- final chisq = %.6e --- \n", Chisq);
    }
    LOG_WRITE(MeSS, true, LM->verboselv>1);
    LOG_MODEL(Para->Par, LM->verboselv>1);

    Para->Chisq = Chisq;

    return it;
}

/*--------------------------------------------------------------------------------*/

static int INVERSION(STRUCT_STK *Stk, STRUCT_PARA *Para, STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        ME inversion of Stokes profiles
      Record of revisions:
        09 Mar. 2026
      Input parameters:
        Stk, a structure storing the Stokes profiles
        Para, a structure storing the model parameters
        LM, a structure storing the hessian matrix.
      Output:
        Para, return the best fit parameters.
        Stk, return the best fit profiles if needed.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
  
    if(Init_Guess(Stk, Para)){
      Para->Chisq_Best = -1.;
      Stk->normp = 1.;
      for(int ipar=0; ipar<9; ipar++){
        Para->Par_Best[ipar] = 0;
      }
      return 0;
    }

    Noise_Init(Stk);
    int irun;
    for(irun=0; irun<LM->nruns*4; irun++){
      
      sprintf(MeSS, "\n ### inversion run: %d  \n", irun);
      LOG_WRITE(MeSS, true, LM->verboselv>1);
      if(irun>0){ 
        Random_Jump(Para, LM, irun);
        sprintf(MeSS, "\n  ramdom jumping: %d \n", irun);
        LOG_WRITE(MeSS, true, LM->verboselv>1);
      }

      LM_FIT(Stk, Para, LM);

      if(Para->Chisq<LM->Criteria){
        sav_par(Para, Stk);
        break;
      }else if(irun==0){
        sav_par(Para, Stk);
      }else if(Para->Chisq<Para->Chisq_Best){
        sav_par(Para, Stk);
      }

      if(irun >= LM->nruns*2){
        if(Para->Chisq_Best<LM->Criteria*2) break;
      }
    }

    if(Para->Par_Best[1]>L_Pi){
      Para->Par_Best[1] -= L_Pi;
    }
    if(LM->HMI_REF){
      Para->Par_Best[2] += 0.5*L_Pi;
      if(Para->Par_Best[2]>L_Pi) Para->Par_Best[2] -= L_Pi;
    }

    sprintf(MeSS, "\n  ### best run (%d): ", irun);
    LOG_WRITE(MeSS, true, LM->verboselv>0);
    sprintf(MeSS, " --- Best chisq = %.6e --- \n", Para->Chisq_Best);
    LOG_WRITE(MeSS, true, LM->verboselv>0);
    LOG_MODEL(Para->Par_Best, LM->verboselv>0);

    return irun;
}

/*--------------------------------------------------------------------------------*/

int INVERSION_MULTI(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_PARA *Para, \
    STRUCT_LM *LM, STRUCT_SUBSET *Subset){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        ME inversion of Stokes profiles
      Record of revisions:
        09 Mar. 2026
      Input parameters:
        Input, the input configuration.
        Stk, a structure storing the Stokes profiles
        Para, a structure storing the model parameters
        LM, a structure storing the hessian matrix.
        Subset, a structure storing the pixels to read.
      Output:
        Para, return the best fit parameters.
        Stk, return the best fit profiles if needed.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
  
    const int nw = Stk->nw;
    const int nshift = 4*nw;
    double *ptr = Input->resbuff;
    double *ptrfit = Input->fitbuff;

    Stk->prof = Input->profbuff;
    
    for(int icount=0; icount<Subset->nProf; icount++){
  
      double Imean = 0.0;
      for(int iw=0;iw<nw;iw++){
        Imean += Stk->prof[iw];
      }
      Imean /= nw;

      bool valid = true;
      if(isnan(Imean) || Imean<Input->Icriteria){ 
        valid = false;
      }

      if(valid){
        if(!INVERSION(Stk, Para, LM)) valid = false;
      }
    
      if(valid){
        for(int ipar=0; ipar<9; ipar++){
          ptr[ipar] = Para->Par_Best[ipar];
        }
        ptr[9] = Para->Chisq_Best/Stk->normp;

        if(Input->output_fit){
          memcpy((void *)ptrfit, (void *)Stk->fit_best, nshift*sizeof(double));
        }
      }else{
        memset(ptr, 0, 10*sizeof(double));
        if(Input->output_fit) memset(ptrfit, 0, nshift*sizeof(double));
      }

      ptr += 10;
      Stk->prof += nshift;
      if(Input->output_fit) ptrfit += nshift;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

