
#include "LM_FIT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:

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

static int sav_par(STRUCT_PAR *Par, STRUCT_STK *Stk);

static int Noise_Init(STRUCT_STK *Stk);

static double Loss_Function(STRUCT_STK *Stk);

static int Par_update(double *Par, double *Sol, STRUCT_INPUT *Input);

static int Hessian_Compute(STRUCT_STK *Stk, STRUCT_LM *LM, \
        STRUCT_INPUT *Input);

static double Trial_Synthesis(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
        STRUCT_PAR *Par, STRUCT_LM *LM);

static bool Lambda_propose(STRUCT_LM *LM, bool accepted);

static bool converg_check(double Chisq, double Chisq_new, \
        STRUCT_INPUT *Input);

static int Random_Jump(STRUCT_PAR *Par, int indx, STRUCT_MPI *Mpi, 
        STRUCT_INPUT *Input);

static int Get_Regul(STRUCT_INPUT *Input, double *Par, STRUCT_LM *LM, 
        bool matrix_flag);

static int SVD_solve_mem(double **A, double *B, double *X, \
        STRUCT_INPUT *Input);

/*--------------------------------------------------------------------------------*/

static int sav_par(STRUCT_PAR *Par, STRUCT_STK *Stk){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Save the best fit parameters.
      Record of revisions:
        11 Apr. 2025 (Hao Li)
      Input parameters:
        Par, structure with parameters.
      Output parameter:
        Par, structure with parameters.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double *p, **pp;

    p = Par->Par_Best;
    Par->Par_Best = Par->Par;
    Par->Par = p;
    
    pp = Stk->fit;
    Stk->fit = Stk->syn;
    Stk->syn = pp;

    Par->Chisq_Best = Par->Chisq;

    return  0;

}

/*--------------------------------------------------------------------------------*/

static int Noise_Init(STRUCT_STK *Stk){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Initialize the noise.
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:
        Stk, structure with the Stokes profiles.
      Output parameter:
        Stk, structure with the Stokes profiles.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    int i;
    
    for(i=0;i<Stk->nl;i++){
      Stk->noise[i] = 4.*Stk->norm;
    }
    
    /*
    double sum  = 0;

    for(i=0;i<Stk->nl;i++){
      Stk->noise[i] = sqrt(Stk->prof[0][i]);
      //Stk->noise[i] = Stk->prof[0][i];

      sum += Stk->noise[i];
    }

    sum /= Stk->nl;

    for(i=0;i<Stk->nl;i++){

      Stk->noise[i] *= 10.*Stk->norm/sum;

    }
    */

    return 0;

}

/*--------------------------------------------------------------------------------*/

static double Loss_Function(STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Compute the merit function.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Input, Input configuration.
        Stk, structure with Stokes profiles.
      Return:
        return the chisq.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int i,j;
    double Chisq = 0;
   
    for(i=0; i<4; i++){
      for(j=0; j<Stk->nl; j++){
        Chisq += Stk->Weights_SQR[i]*(Stk->syn[i][j]-Stk->prof[i][j]) \
            *(Stk->syn[i][j]-Stk->prof[i][j]);
      }
    }    

    Chisq /= 4*Stk->nl*Stk->norm;
    //Chisq /= 4*Stk->nl;

    
    return Chisq;

}

/*--------------------------------------------------------------------------------*/

static bool Lambda_propose(STRUCT_LM *LM, bool accepted){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Propose a lambda factor.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        LM, structure with the hessian matrix.
        accepted, accepted or not.
      Output parameters:
        LM, the damping factor.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(accepted){
      LM->Damp /= LM->Lam_accept;
      if (LM->Damp < LM->Lam_Lim[0]) LM->Damp = LM->Lam_Lim[0];
    }else{
      if (LM->Damp > LM->Lam_Lim[1]/LM->Lam_reject){
        return false;
      }else{
        LM->Damp *= LM->Lam_reject;
      }
    }

    return true;

}

/*--------------------------------------------------------------------------------*/

static bool converg_check(double Chisq, double Chisq_new, \
        STRUCT_INPUT *Input){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        check the convergence.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Chisq, the chisq.
        Chisq_new, the new chisq.
        Input, the input configuration.
      Return:
        return the check.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if ((Chisq -Chisq_new)/Chisq < Input->Convg_Criteria){
      return true;
    }else {
      return false;
    }

}

/*--------------------------------------------------------------------------------*/

static int Random_Jump(STRUCT_PAR *Par, int indx, STRUCT_MPI *Mpi, 
    STRUCT_INPUT *Input){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        jumping the initial azimuth another value.
      Record of revisions:
        17 Apr. 2025 (Hao Li)
      Input parameters:
        Par, structure with the model parameters
        indx, index of the run.
        Mpi, structure with tha random seeds.
      Output parameter:
        Par, structure with the model parameters
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    Par->Par[2] = Par->Par_Guess[2]+Par_Pi/6.*GASDEV(Mpi->idum);

    Par->Par[3] = Par->Par_Guess[3]+Par_Pi/2.*GASDEV(Mpi->idum);

    Par->Par[4] = Par->Par_Guess[4]+5.*GASDEV(Mpi->idum);

    Par->Par[6] = Par->Par_Best[6];

    Par->Par[8] = Par->Par_Best[8];
    
    if(Par->Chisq_Best < Input->Chisq_Criteria*1.5){

      Par->Par[1] = Par->Par_Best[1]+300*GASDEV(Mpi->idum);
      Par->Par[5] = Par->Par_Best[5]+0.2*GASDEV(Mpi->idum);
      Par->Par[7] = Par->Par_Best[7]+10.*GASDEV(Mpi->idum);
      Par->Par[9] = Par->Par_Best[9]+0.2*GASDEV(Mpi->idum);

    }else if(Par->Par_Best[1]>1500 || Par->Par_Guess[1]>1500){ 

      Par->Par[1] = Par->Par_Guess[1]+800*GASDEV(Mpi->idum);
      Par->Par[5] = Par->Par_Guess[5]+0.2*GASDEV(Mpi->idum);
      Par->Par[7] = Par->Par_Guess[7]+10*GASDEV(Mpi->idum);
      Par->Par[9] = Par->Par_Guess[9]+0.2*GASDEV(Mpi->idum);
    }
    
    while(Par->Par[3]<0){
      Par->Par[3] += Par_Pi;
    }
    
    while(Par->Par[3]>Par_Pi){
      Par->Par[3] -= Par_Pi;
    }

    bounds_check(Par->Par, Par);

    return 0;

}

/*--------------------------------------------------------------------------------*/

static int Par_update(double *Par, double *Sol, STRUCT_INPUT *Input){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Update the parameters according to the solution.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Par, the original parameters.
        Sol, the solutions.
        Input, Input configuration.
      Output parameters:
        Par, the updated parameters.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int i, ipar;

    if(Input->npar==9){
      for(i=1; i<10; i++){
        Par[i] += Sol[i]*Input->step[i];
      }
    }else{
      ipar = 1;
      for(i=1; i<10; i++){
        if(!Input->inv[i]) continue;
        Par[i] += Sol[ipar]*Input->step[i];
        ipar++;
      }
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

static int Hessian_Compute(STRUCT_STK *Stk, STRUCT_LM *LM, \
        STRUCT_INPUT *Input){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Compute the Hessian matrix.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Stk, structure with Stokes profiles.
        LM, structure with hessian matrix.
        Input, Input configuration.
      Output parameters:
        LM, structure with hessian matrix.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int i, j, k, l, ipar, jpar;

    /*
    long Byte = Input->npar*sizeof(double);
    memset(&LM->Jacfvec[1], 0, Byte);
    Byte *= Input->npar;
    memset(&LM->Hessian[1][1], 0, Byte);
    */

    if(Input->npar==9){
      for(i=1; i<10; i++){
        for(j=i; j<10; j++){
          LM->Hessian[i][j] = 0;
          for(k=0; k<4; k++){
            for(l=0; l<Stk->nl; l++){
              LM->Hessian[i][j] += Stk->Weights_SQR[k] \
                  *Stk->Jacobian[k][i][l]*Stk->Jacobian[k][j][l] \
                  /Stk->noise[l];
            }
          }
          LM->Hessian[j][i] = LM->Hessian[i][j];
        }

        LM->Jacfvec[i] = 0;
        for(k=0; k<4; k++){
          for(l=0; l<Stk->nl; l++){
            LM->Jacfvec[i] -= Stk->Weights_SQR[k] \
                *Stk->Jacobian[k][i][l]*(Stk->syn[k][l]-Stk->prof[k][l]) \
                /Stk->noise[l];
          }
        }
      }

    }else{
      ipar = 1;
      for(i=1; i<10; i++){
        if(!Input->inv[i]) continue;
        jpar = ipar;
        for(j=i; j<10; j++){
          if(!Input->inv[j]) continue;
          LM->Hessian[ipar][jpar] = 0;
          for(k=0; k<4; k++){
            for(l=0; l<Stk->nl; l++){
              LM->Hessian[ipar][jpar] += Stk->Weights_SQR[k] \
                  *Stk->Jacobian[k][i][l]*Stk->Jacobian[k][j][l] \
                  /Stk->noise[l];
            }
          }
          LM->Hessian[jpar][ipar] = LM->Hessian[ipar][jpar];
          jpar++; 
        }

        LM->Jacfvec[ipar] = 0;
        for(k=0; k<4; k++){
          for(l=0; l<Stk->nl; l++){
            LM->Jacfvec[ipar] -= Stk->Weights_SQR[k] \
                *Stk->Jacobian[k][i][l]*(Stk->syn[k][l]-Stk->prof[k][l]) \
                /Stk->noise[l];
          }
        }
        ipar++;
      }
      
    }
  
    return  0;

}

/*--------------------------------------------------------------------------------*/

static int Get_Regul(STRUCT_INPUT *Input, double *Par, STRUCT_LM *LM, 
        bool matrix_flag){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Computes the regularization.
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:        
        Input, structure with the input information.
        Par, structure with model parameters.
        matrix_flag, a flag.
      Output parameters: 
        LM, structure with the Levenberg-Marquardt configuration.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int i, ipar = 0;
    double LL = 1., Resi = 0;
    LM->Penalty = 0.0;

    for(i=1; i<10; i++){
      if(!Input->inv[i]) continue;
      ipar++;
      if(!Input->regl[i]) continue;
      Resi = Par[i] - Input->value_const[i];

      LM->Penalty += Resi*Resi*Input->Regul_weight[i];

      if(matrix_flag){
        LM->Regul_J[ipar] = LL*Resi*Input->Regul_weight[i];
        LM->Regul_H[ipar][ipar] = LL*LL*Input->Regul_weight[i];
        LM->Jacfvec[ipar] -= LM->Regul_J[ipar];
        LM->Hessian[ipar][ipar] += LM->Regul_H[ipar][ipar];
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

static int SVD_solve_mem(double **A, double *B, double *X, \
        STRUCT_INPUT *Input){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Solving linear equations (A[][]*X[]=B[]) with svd decomposition.
      Record of revisions:
        10 May. 2024 (Hao Li)
      Input parameters:        
        A[][], the A matrix.
        B[], the B vector.
        Input, the input configuration.
      Output parameters:
        X[], The solution.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double Threshold = Input->Threshold, WMAX, Wthreshold;
    int i;

    svdcmp(A, Input->npar, Input->npar, Input->W, Input->V);
    
    WMAX = Input->W[1];

    for (i=2; i<=Input->npar; i++){ 
      if(WMAX<Input->W[i]) WMAX = Input->W[i];
    }


    Wthreshold = WMAX*Threshold;
    for(i=1; i<=Input->npar; i++){ 
      if(Input->W[i]<Wthreshold){
        Input->W[i] = 0.;
      }
    }

    svbksb(A, Input->W, Input->V, Input->npar, Input->npar, B, X);

    double maxstep = 0.5;
    for (i=1; i<=Input->npar; i++){ 
      if(X[i]>maxstep){ 
        X[i] = maxstep;
      }else if(X[i]<-maxstep){ 
        X[i] = -maxstep;
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

static double Trial_Synthesis(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
        STRUCT_PAR *Par, STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Trial synthesis
      Record of revisions:
        17 Apr. 2024 (Hao Li)
      Input parameters:        
        Input, the input configuration.
        Stk, structure with the Stokes profiles
        Par, structure with the model parameters
        LM, structure with the hessian matrix.
      Output parameters:
        Stk, structure with the synthesized Stokes profiles.
      return:
        return the chisq.
    ######################################################################*/

/*----------------------------------------------------------------------------*/  

    int i, j;
    double chisq;

    for(i=1;i<=Input->npar;i++){
      for(j=1;j<=Input->npar;j++){
        LM->Hessian_new[i][j] = LM->Hessian[i][j];
      }
    }
    for(i=1;i<=Input->npar;i++){
      LM->Hessian_new[i][i] *= (1.+LM->Damp);
      LM->Sol[i] = LM->Jacfvec[i];
    }

    SVD_solve_mem(LM->Hessian_new, LM->Jacfvec, LM->Sol, Input);

    for(i=1;i<10;i++){
      Par->Par_tmp[i] = Par->Par[i];
    }

    Par_update(Par->Par_tmp, LM->Sol, Input);

    bounds_check(Par->Par_tmp, Par);

    for(i=0;i<Input->nline;i++){
      Input->Lines[i].Par = Par->Par_tmp;
    }
    Milne_Eddington(Input->Lines, Input->nline, Stk, Input->Fadd, \
        Input, false);

    for(i=0;i<Input->nline;i++){
      Input->Lines[i].Par = Par->Par;
    }

    chisq = Loss_Function(Stk);

    if(Input->Regul){
      Get_Regul(Input, Par->Par_tmp, LM, false);
      chisq += LM->Penalty;
    }

    return  chisq;


}

/*--------------------------------------------------------------------------------*/

extern int LM_FIT(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_PAR *Par, \
        STRUCT_LM *LM){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        the Levenberg Marquardt.
      Record of revisions:
        18 Jun. 2024 (Hao Li)
      Input parameters:
        Input, the input configuration.
        Stk, structure with the Stokes profiles
        Par, structure with the model parameters
        LM, structure with the hessian matrix.
      Return:
        return the num of iteration check.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    
    double Chisq, Chisq_new;
    bool accepted = false;
    int i = 0, j;

    for(i=0;i<Input->nline;i++){
      Input->Lines[i].Par = Par->Par; 
    }
    
    Milne_Eddington(Input->Lines, Input->nline, Stk, Input->Fadd, Input, true);

    Chisq = Loss_Function(Stk);

    if(Input->Regul){
      Get_Regul(Input, Par->Par, LM, false);
      Chisq += LM->Penalty;
      sprintf(MeSS, " *** Init chisq = %.6e; Penality = %.6e *** \n", \
          Chisq, LM->Penalty);

    }else{
      sprintf(MeSS, " *** Init chisq = %.6e; *** \n", Chisq);

    }

    Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=1);
    Verbose_model(Par->Par, Input->Verbose_Path, Input->verboselv<=1);

    for(i=0; i<Input->niter; i++){

      accepted = false;

      sprintf(MeSS, " *** iteration : %d  \n", i);
      Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=2);

      Hessian_Compute(Stk, LM, Input);

      if(Input->Regul){
        Get_Regul(Input, Par->Par, LM, true);
      }

      do{

        Chisq_new = Trial_Synthesis(Input, Stk, Par, LM);

        if (Chisq_new < Chisq){

          if(Input->Regul){
            sprintf(MeSS, "    * trail accepted. chisq = %.6e; " \
                "damping = %.6e; Penalty = %.6e *\n", Chisq_new, \
                LM->Damp, LM->Penalty);
          }else{
            sprintf(MeSS, "    * trail accepted. chisq = %.6e; " \
                "damping = %.6e; *\n", Chisq_new, LM->Damp);
          }
          Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=2);

          accepted = true;

          for(j=1;j<10;j++){
            Par->Par[j] = Par->Par_tmp[j];
          }
              
          Milne_Eddington(Input->Lines, Input->nline, Stk, Input->Fadd, \
              Input, true);

          Lambda_propose(LM, accepted);

        }else{
          if(Input->Regul){
            sprintf(MeSS, "    * trail chisq = %.6e; damping = %.6e;" \
                " Penalty = %.6e *\n", \
                Chisq_new, LM->Damp, LM->Penalty);
          }else{
            sprintf(MeSS, "    * trail chisq = %.6e; damping = %.6e" \
                " *\n", Chisq_new, LM->Damp);
          }
          Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=2);

          accepted =false;
          if(!Lambda_propose(LM, accepted)) break;
        }

      }while(!accepted);

      sprintf(MeSS, "    * iteration : %d  chisq : %.6e -> %.6e *\n", i, \
          Chisq, Chisq_new);
      Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=2);

      if(accepted){

        if(converg_check(Chisq, Chisq_new, Input)){
          Chisq = Chisq_new;
          break;
        } 
        Chisq = Chisq_new;
      }else{
        break;
      }
    }

    if(Input->Regul){
      sprintf(MeSS, " --- final chisq = %.6e; Penalty = %.6e --- \n", \
          Chisq, LM->Penalty);
    }else{
      sprintf(MeSS, " --- final chisq = %.6e --- \n", Chisq);
    }
    Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=1);
    Verbose_model(Par->Par, Input->Verbose_Path, Input->verboselv<=1);

    Par->Chisq = Chisq;

    return i;
}

/*--------------------------------------------------------------------------------*/

extern int INVERSION(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
        STRUCT_PAR *Par, STRUCT_LM *LM, STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        ME inversion of Stokes profiles
      Record of revisions:
        18 Apr. 2025 (Hao Li)
      Input parameters:
        Input, the input configuration.
        Stk, structure with the Stokes profiles
        Par, structure with the model parameters
        LM, structure with the hessian matrix.
        Mpi, structure with MPI information.
      Output:
        Par, return the best fit parameters.
        Stk, return the best fit profiles if needed.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
  
    int i, nrun;

    Init_Guess(Stk, Par, Input, LM);

    Noise_Init(Stk);

    if(Par->Par_Guess[1]>400){ 
      nrun = Input->nrun*2;
    }else{
      nrun = Input->nrun;
    }

    for(i=0; i<nrun*2; i++){
      
      sprintf(MeSS, "\n ### inversion run: %d  \n", i);
      Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=1);
      if(i > 0){ 
        Random_Jump(Par, i, Mpi, Input);
        sprintf(MeSS, "\n  ramdom jumping: %d \n", i);
        Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=1);
      }

      LM_FIT(Input, Stk, Par, LM);

      if(Par->Chisq < Input->Chisq_Criteria){
        sav_par(Par, Stk);
        break;
      }else if (i == 0){
        sav_par(Par, Stk);
      }else if (Par->Chisq < Par->Chisq_Best){
        sav_par(Par, Stk);
      }

      if(i>=nrun){
        if(Par->Chisq_Best<Input->Chisq_Criteria*2) break;
      }
    }

    if(Par->Par_Best[3]>Par_Pi){
      Par->Par_Best[3] -= Par_Pi;
    }
    if(Input->Azimuth_Rotate){
      Par->Par_Best[3] += 0.5*Par_Pi;
      if(Par->Par_Best[3]>Par_Pi) Par->Par_Best[3] -= Par_Pi;
    }

    sprintf(MeSS, "\n  ### best run : ");
    Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=1);
    sprintf(MeSS, " --- Best chisq = %.6e --- \n", Par->Chisq_Best);
    Verbose(MeSS, Input->Verbose_Path, Input->verboselv<=1);
    Verbose_model(Par->Par_Best, Input->Verbose_Path, Input->verboselv<=1);

    return i;
}

/*--------------------------------------------------------------------------------*/

