
#include "MILNE_EDDINGTON.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        28 Jun. 2024
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

static double Parabolic(double *X, double *Y, int Indx);

/*--------------------------------------------------------------------------------*/

static double Parabolic(double *X, double *Y, int Indx){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        3 points Parabolic interpolation to get the minimum position.
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:        
        X, array with x values.
        Y, array with y values.
        Indx, the index of minimum Y.
      Return:
        The X position corresponding to the Y minimum .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double X_tmp1 = X[Indx] - X[Indx-1];
    double X_tmp2 = X[Indx] - X[Indx+1];
    double Y_tmp1 = Y[Indx] - Y[Indx-1];
    double Y_tmp2 = Y[Indx] - Y[Indx+1];

    return X[Indx]-0.5*(X_tmp1*X_tmp1*Y_tmp2-X_tmp2*X_tmp2*Y_tmp1) \
        /(X_tmp1*Y_tmp2-X_tmp2*Y_tmp1);
    
}

/*--------------------------------------------------------------------------------*/
    
extern int Milne_Eddington(STRUCT_MELINE *Lines, int nline, \
        STRUCT_STK *Stk, STRUCT_FADDEEVA *Fadd, STRUCT_INPUT *Input, 
        bool Deriv){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Calculate the Stokes profiles under M-E atmosphere model (normal 
            Zeeman effect)
      Record of revisions:
        17 Jun. 2024
      Input parameters:
        Lines, structure with the line information.
        nline, number of lines.
        Fadd, structure with some precomputed coefficients.
        Deriv, output the derivatives of not.
      Output:
        Stk, structure with Stokes profiles
      Reference:
        LL04 5.47, LL04 9.4 9.8, del Toro Iniesta 2003
      Note:
        Par[1] Bmod; Par[2] ThetaB; Par[3] PhiB; Par[4] Vlos; Par[5] Dopp;
        Par[6] Damp; Par[7] Eta; Par[8] S0+S1; Par[9] S0/(S0+S1);
        Par[1], Par[2], Par[3], Par[4], Par[8], Par[9] must be same for 
        all lines. 
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
  
    // precompute some sin and cos functions
    double SIN_THB = sin(Lines[0].Par[2]);
    double COS_THB = cos(Lines[0].Par[2]);
    double SIN_THB_SQ = SIN_THB*SIN_THB;
    double COS_THB_SQ = COS_THB*COS_THB;
    double COS2PHI = cos(2.*Lines[0].Par[3]);
    double SIN2PHI = sin(2.*Lines[0].Par[3]);

    // memory allocation
    double *phi = (double *)VECTOR(-1, 1, enum_dbl, false);
    double *psi = (double *)VECTOR(-1, 1, enum_dbl, false);
    double *eta = (double *)VECTOR(0, 3, enum_dbl, true);
    double *rho = (double *)VECTOR(0, 3, enum_dbl, true);
    double *ShiftB = (double *)malloc(sizeof(double)*nline);
    double *ShiftV = (double *)malloc(sizeof(double)*nline);

    // Atom->Lambda in A, Lines[0].Par[5] in mA, Atom->BShift in mA
    // Atom->BShift = 4.6686e-10*Atom->Geffect*Atom->Lambda*Atom->Lambda;
    int i, j, q;

    for(i=0;i<nline;i++){
      ShiftB[i] = Lines[i].Par[1]*Lines[i].BShift/Lines[i].Par[5];
      ShiftV[i] = 1e3*(-1e3*Lines[0].Par[4]/Par_C*Lines[i].Lambda0 \
        -Lines[i].Lambda0)/Lines[i].Par[5];
    }

    // get the source function
    double B0 = Lines[0].Par[8]*Lines[0].Par[9];
    double B1 = Lines[0].Par[8]-B0;
    
    double Delta = 0, temp, temp1, temp2, temp3;
    double Lam_Shifts, ShiftW, H, L, Hp_lam, Hp_damp, Lp_lam, Lp_damp;

    double *Deltap = NULL, *tempp = NULL, *temp1p = NULL, *temp2p = NULL;
    double *temp3p = NULL, **etap = NULL, **rhop = NULL, **phip = NULL;
    double **psip = NULL;

    if(Deriv){
      phip = (double **)MATRIX(-1, 1, 1, 6, enum_dbl, false);
      psip = (double **)MATRIX(-1, 1, 1, 6, enum_dbl, false);
      Deltap = (double *)VECTOR(1, 7, enum_dbl, false);
      tempp = (double *)VECTOR(1, 7, enum_dbl, false);
      temp1p = (double *)VECTOR(1, 7, enum_dbl, false);
      temp2p = (double *)VECTOR(1, 7, enum_dbl, false);
      temp3p = (double *)VECTOR(1, 7, enum_dbl, false);
      etap = (double **)MATRIX(0, 3, 1, 7, enum_dbl, false);
      rhop = (double **)MATRIX(1, 3, 1, 7, enum_dbl, false);
    }

    for (j = 0; j<Stk->nl; j++){
      for(i=0;i<nline;i++){

        ShiftW = 1e3*Stk->Lambda[j]/Lines[i].Par[5]+ShiftV[i];

        // compute the phi and psi
        for (q=-1; q<=1; q++){
          // q=Mu-Ml,(1 sigma bule, -1 sigma red, 0 pi)
          Lam_Shifts = ShiftW+q*ShiftB[i];
          Faddeeva(Lam_Shifts, Lines[i].Par[6], &H, &L, Fadd);

          if(i==0){
            phi[q] = H/Par_SqrtPi;
            psi[q] = L/Par_SqrtPi;
          }else{
            phi[q] += H/Par_SqrtPi;
            psi[q] += L/Par_SqrtPi;
          }

          if(Deriv){
            Hp_lam = 2*(-Lam_Shifts*H+Lines[i].Par[6]*L);
            Lp_damp = Hp_lam;
            Hp_damp = 2*(-1/Par_SqrtPi+Lines[i].Par[6]*H+Lam_Shifts*L);
            Lp_lam = -Hp_damp;
            if(i==0){
              phip[q][1] = q*Lines[i].BShift/Lines[i].Par[5] \
                  *Hp_lam/Par_SqrtPi;
              psip[q][1] = q*Lines[i].BShift/Lines[i].Par[5] \
                  *Lp_lam/Par_SqrtPi;
              phip[q][4] = -1e6/Par_C*Lines[i].Lambda0/Lines[i].Par[5] \
                  *Hp_lam/Par_SqrtPi;
              psip[q][4] = -1e6/Par_C*Lines[i].Lambda0/Lines[i].Par[5] \
                  *Lp_lam/Par_SqrtPi;
              phip[q][5] = -Lam_Shifts/Lines[i].Par[5]*Hp_lam/Par_SqrtPi;
              psip[q][5] = -Lam_Shifts/Lines[i].Par[5]*Lp_lam/Par_SqrtPi;
              phip[q][6] = Hp_damp/Par_SqrtPi;
              psip[q][6] = Lp_damp/Par_SqrtPi;
            }else{
              phip[q][1] += q*Lines[i].BShift/Lines[i].Par[5] \
                  *Hp_lam/Par_SqrtPi;
              psip[q][1] += q*Lines[i].BShift/Lines[i].Par[5] \
                  *Lp_lam/Par_SqrtPi;
              phip[q][4] += -1e6/Par_C*Lines[i].Lambda0/Lines[i].Par[5] \
                  *Hp_lam/Par_SqrtPi;
              psip[q][4] += -1e6/Par_C*Lines[i].Lambda0/Lines[i].Par[5] \
                  *Lp_lam/Par_SqrtPi;
              phip[q][5] += -Lam_Shifts/Lines[i].Par[5]*Hp_lam/Par_SqrtPi;
              psip[q][5] += -Lam_Shifts/Lines[i].Par[5]*Lp_lam/Par_SqrtPi;
              phip[q][6] += Hp_damp/Par_SqrtPi;
              psip[q][6] += Lp_damp/Par_SqrtPi;
            }       
          }
        }
      }

      // Stokes profiles LL04 Page 414 Eq 9.110
      eta[0] = 0.5*Lines[0].Par[7]*(phi[0]*SIN_THB_SQ \
          +0.5*(phi[-1]+phi[1])*(1+COS_THB_SQ));
      eta[1] = 0.5*Lines[0].Par[7]*(phi[0]-0.5*(phi[-1]+phi[1])) \
          *SIN_THB_SQ*COS2PHI;
      eta[2] = 0.5*Lines[0].Par[7]*(phi[0]-0.5*(phi[-1]+phi[1])) \
          *SIN_THB_SQ*SIN2PHI;
      eta[3] = 0.5*Lines[0].Par[7]*(phi[-1]-phi[1])*COS_THB;
      rho[1] = 0.5*Lines[0].Par[7]*(psi[0]-0.5*(psi[-1]+psi[1])) \
          *SIN_THB_SQ*COS2PHI;
      rho[2] = 0.5*Lines[0].Par[7]*(psi[0]-0.5*(psi[-1]+psi[1])) \
          *SIN_THB_SQ*SIN2PHI;
      rho[3] = 0.5*Lines[0].Par[7]*(psi[-1]-psi[1])*COS_THB;

      // some coefficients
      temp = 1+eta[0];
      temp1 = eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
      temp2 = rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
      temp3 = eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
      Delta = (temp*temp*(temp*temp-temp1+temp2)-temp3*temp3);

      // conpute the derivatives of eta and rho with respect 
      // to some parameters
      if(Deriv){
        //der  bmod
        etap[0][1] = 0.5*Lines[0].Par[7]*(phip[0][1]*SIN_THB_SQ \
            +0.5*(phip[-1][1]+phip[1][1])*(1+COS_THB_SQ));
        etap[1][1] = 0.5*Lines[0].Par[7]*(phip[0][1] \
            -0.5*(phip[-1][1]+phip[1][1]))*SIN_THB_SQ*COS2PHI;
        etap[2][1] = 0.5*Lines[0].Par[7]*(phip[0][1] \
            -0.5*(phip[-1][1]+phip[1][1]))*SIN_THB_SQ*SIN2PHI;
        etap[3][1] = 0.5*Lines[0].Par[7]*(phip[-1][1]-phip[1][1])*COS_THB;
        rhop[1][1] = 0.5*Lines[0].Par[7]*(psip[0][1] \
            -0.5*(psip[-1][1]+psip[1][1]))*SIN_THB_SQ*COS2PHI;
        rhop[2][1] = 0.5*Lines[0].Par[7]*(psip[0][1] \
            -0.5*(psip[-1][1]+psip[1][1]))*SIN_THB_SQ*SIN2PHI;
        rhop[3][1] = 0.5*Lines[0].Par[7]*(psip[-1][1]-psip[1][1])*COS_THB;

        //der thetab
        etap[0][2] = 0.5*Lines[0].Par[7]*(phi[0]*(2.0*SIN_THB*COS_THB) \
            +0.5*(phi[-1]+phi[1])*(-2.0*COS_THB*SIN_THB));
        etap[1][2] = 0.5*Lines[0].Par[7]*(phi[0]-0.5*(phi[-1]+phi[1])) \
            *(2.0*SIN_THB*COS_THB)*COS2PHI;
        etap[2][2] = 0.5*Lines[0].Par[7]*(phi[0]-0.5*(phi[-1]+phi[1])) \
            *(2.0*SIN_THB*COS_THB)*SIN2PHI;
        etap[3][2] = 0.5*Lines[0].Par[7]*(phi[-1]-phi[1])*(-SIN_THB);
        rhop[1][2] = 0.5*Lines[0].Par[7]*(psi[0]-0.5*(psi[-1]+psi[1])) \
            *(2.0*SIN_THB*COS_THB)*COS2PHI;
        rhop[2][2] = 0.5*Lines[0].Par[7]*(psi[0]-0.5*(psi[-1]+psi[1])) \
            *(2.0*SIN_THB*COS_THB)*SIN2PHI;
        rhop[3][2] = 0.5*Lines[0].Par[7]*(psi[-1]-psi[1])*(-SIN_THB);

        //der phib
        etap[0][3] = 0;
        etap[1][3] = 0.5*Lines[0].Par[7]*(phi[0]-0.5*(phi[-1]+phi[1])) \
            *SIN_THB_SQ*(-2.0*SIN2PHI);
        etap[2][3] = 0.5*Lines[0].Par[7]*(phi[0]-0.5*(phi[-1]+phi[1])) \
            *SIN_THB_SQ*(2.0*COS2PHI);
        etap[3][3] = 0; 
        rhop[1][3] = 0.5*Lines[0].Par[7]*(psi[0]-0.5*(psi[-1]+psi[1])) \
            *SIN_THB_SQ*(-2.0*SIN2PHI);
        rhop[2][3] = 0.5*Lines[0].Par[7]*(psi[0]-0.5*(psi[-1]+psi[1])) \
            *SIN_THB_SQ*(2.0*COS2PHI);
        rhop[3][3] = 0;

        //der to velocit  doppler width  dampping parameter
        for(i = 4; i < 7; i++){
          etap[0][i] = 0.5*Lines[0].Par[7]*(phip[0][i]*SIN_THB_SQ \
              +0.5*(phip[-1][i]+phip[1][i])*(1+COS_THB_SQ));
          etap[1][i] = 0.5*Lines[0].Par[7]*(phip[0][i] \
              -0.5*(phip[-1][i]+phip[1][i]))*SIN_THB_SQ*COS2PHI;
          etap[2][i] = 0.5*Lines[0].Par[7]*(phip[0][i] \
              -0.5*(phip[-1][i]+phip[1][i]))*SIN_THB_SQ*SIN2PHI;
          etap[3][i] = 0.5*Lines[0].Par[7]*(phip[-1][i]-phip[1][i])*COS_THB;
          rhop[1][i] = 0.5*Lines[0].Par[7]*(psip[0][i] \
              -0.5*(psip[-1][i]+psip[1][i]))*SIN_THB_SQ*COS2PHI;
          rhop[2][i] = 0.5*Lines[0].Par[7]*(psip[0][i] \
              -0.5*(psip[-1][i]+psip[1][i]))*SIN_THB_SQ*SIN2PHI;
          rhop[3][i] = 0.5*Lines[0].Par[7]*(psip[-1][i]-psip[1][i])*COS_THB;
        }

        etap[0][7] = eta[0]/Lines[0].Par[7];
        for(i = 1; i < 4; i++){
          etap[i][7] = eta[i]/Lines[0].Par[7];
          rhop[i][7] = rho[i]/Lines[0].Par[7];
        }

        for(i = 1; i < 8; i++){
          tempp[i] = etap[0][i];
          temp1p[i] = 2.0*eta[1]*etap[1][i]+2.0*eta[2]*etap[2][i] \
              +2.0*eta[3]*etap[3][i];
          temp2p[i] = 2.0*rho[1]*rhop[1][i]+2.0*rho[2]*rhop[2][i] \
              +2.0*rho[3]*rhop[3][i];
          temp3p[i] = etap[1][i]*rho[1]+eta[1]*rhop[1][i] \
              +etap[2][i]*rho[2]+eta[2]*rhop[2][i]+etap[3][i]*rho[3] \
              +eta[3]*rhop[3][i];
          Deltap[i] = 2*temp*tempp[i]*(temp*temp-temp1+temp2) \
              +temp*temp*(2*temp*tempp[i]-temp1p[i]+temp2p[i]) \
              -2.0*temp3*temp3p[i];
        }
      }

      Stk->syn[0][j] = B0+B1/Delta*(temp*(temp*temp+temp2));
      Stk->syn[1][j] = -B1/Delta*(temp*temp*eta[1] \
          +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3);
      Stk->syn[2][j] = -B1/Delta*(temp*temp*eta[2] \
          +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3);
      Stk->syn[3][j] = -B1/Delta*(temp*temp*eta[3] \
          +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3);

      // get the jacobian
      if(Deriv){
        for(i = 1; i < 8; i++){
          Stk->Jacobian[0][i][j] = (B1/Delta/Delta*(Delta \
              *(tempp[i]*(temp*temp+temp2)+temp*(2*temp*tempp[i]+temp2p[i])) \
              -Deltap[i]*(temp*(temp*temp+temp2))))*Input->step[i];
          Stk->Jacobian[1][i][j] = (-B1/Delta/Delta*(Delta \
              *(2*temp*tempp[i]*eta[1]+temp*temp*etap[1][i] \
              +tempp[i]*(eta[3]*rho[2]-eta[2]*rho[3]) \
              +temp*(etap[3][i]*rho[2]+eta[3]*rhop[2][i] \
              -etap[2][i]*rho[3]-eta[2]*rhop[3][i]) \
              +rhop[1][i]*temp3+rho[1]*temp3p[i]) \
              -Deltap[i]*(temp*temp*eta[1]+temp*(eta[3]*rho[2] \
              -eta[2]*rho[3])+rho[1]*temp3)))*Input->step[i];
          Stk->Jacobian[2][i][j] = (-B1/Delta/Delta*(Delta \
              *(2*temp*tempp[i]*eta[2]+temp*temp*etap[2][i] \
              +tempp[i]*(eta[1]*rho[3]-eta[3]*rho[1]) \
              +temp*(etap[1][i]*rho[3]+eta[1]*rhop[3][i] \
              -etap[3][i]*rho[1]-eta[3]*rhop[1][i]) \
              +rhop[2][i]*temp3+rho[2]*temp3p[i]) \
              -Deltap[i]*(temp*temp*eta[2]+temp*(eta[1]*rho[3] \
              -eta[3]*rho[1])+rho[2]*temp3)))*Input->step[i];
          Stk->Jacobian[3][i][j] = (-B1/Delta/Delta*(Delta \
              *(2*temp*tempp[i]*eta[3]+temp*temp*etap[3][i] \
              +tempp[i]*(eta[2]*rho[1]-eta[1]*rho[2]) \
              +temp*(etap[2][i]*rho[1]+eta[2]*rhop[1][i] \
              -etap[1][i]*rho[2]-eta[1]*rhop[2][i]) \
              +rhop[3][i]*temp3+rho[3]*temp3p[i]) \
              -Deltap[i]*(temp*temp*eta[3]+temp*(eta[2]*rho[1] \
              -eta[1]*rho[2])+rho[3]*temp3)))*Input->step[i];
        }

        Stk->Jacobian[0][8][j] = (Lines[0].Par[9]+(1-Lines[0].Par[9])/Delta \
            *(temp*(temp*temp+temp2)))*Input->step[8];
        Stk->Jacobian[1][8][j] = (-(1-Lines[0].Par[9])/Delta*(temp*temp*eta[1] \
            +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3))*Input->step[8];
        Stk->Jacobian[2][8][j] = (-(1-Lines[0].Par[9])/Delta*(temp*temp*eta[2] \
            +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3))*Input->step[8];
        Stk->Jacobian[3][8][j] = (-(1-Lines[0].Par[9])/Delta*(temp*temp*eta[3] \
            +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3))*Input->step[8];

        Stk->Jacobian[0][9][j] = ((Lines[0].Par[8]-Lines[0].Par[8]/Delta \
            *(temp*(temp*temp+temp2))))*Input->step[9];
        Stk->Jacobian[1][9][j] = (Lines[0].Par[8]/Delta*(temp*temp*eta[1] \
            +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3))*Input->step[9];
        Stk->Jacobian[2][9][j] = (Lines[0].Par[8]/Delta*(temp*temp*eta[2] \
            +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3))*Input->step[9];
        Stk->Jacobian[3][9][j] = (Lines[0].Par[8]/Delta*(temp*temp*eta[3] \
            +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3))*Input->step[9];
       
      }
    }
    
    // free the memory
    if(Deriv){
      FREE_MATRIX((void *)phip, -1, 1, enum_dbl);
      FREE_MATRIX((void *)psip, -1, 1, enum_dbl);
      FREE_MATRIX((void *)etap, 0, 1, enum_dbl);
      FREE_MATRIX((void *)rhop, 1, 1, enum_dbl);
      FREE_VECTOR((void *)Deltap, 1, enum_dbl);
      FREE_VECTOR((void *)tempp, 1, enum_dbl);
      FREE_VECTOR((void *)temp1p, 1, enum_dbl);
      FREE_VECTOR((void *)temp2p, 1, enum_dbl);
      FREE_VECTOR((void *)temp3p, 1, enum_dbl);
    }

    FREE_VECTOR(phi, -1, enum_dbl);
    FREE_VECTOR(psi, -1, enum_dbl);
    FREE_VECTOR(eta, 0, enum_dbl);
    FREE_VECTOR(rho, 0, enum_dbl);

    free(ShiftB);
    free(ShiftV);

    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Init_Guess(STRUCT_STK *Stk, STRUCT_PAR *Par, \
        STRUCT_INPUT *Input, STRUCT_LM *LM){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Initialize the model patameters
      Record of revisions:
        18 Jun. 2024 (Hao Li)
      Input parameters:
        Stk, structure with Stokes profiles.
        Par, structure with parameters.
        Input, the input configuration.
        LM, structure with the hessian matrix.
      Output:
        Par, structure with parameters
    ######################################################################*/

/*--------------------------------------------------------------------------------*/  

    int i, iV1 = 0, iV2 = 0, iL = 0; 
    double IMax = Stk->prof[0][0], VMax = Stk->prof[3][0]/Stk->prof[0][0];
    double IMin = IMax, VMin = VMax, delta_v, sumQ = 0., sumU = 0.;
    double dtmp, LP, LV, sumI, sumL = 0.0; 
    double weight = 0.0, Blos = 1.0, Bpos =1.0;

    LP = sqrt(Stk->prof[1][0]*Stk->prof[1][0]+Stk->prof[2][0] \
        *Stk->prof[2][0])/Stk->prof[0][0]; 
    LV = fabs(Stk->prof[3][0]/Stk->prof[0][0]);

    double LPmean = LP, LVmean = LV;

    sumI = Stk->prof[0][0];
    for(i=1; i<Stk->nl; i++){
      sumI += Stk->prof[0][i];
      sumQ += Stk->prof[1][i];
      sumU += Stk->prof[2][i];

      if(IMax < Stk->prof[0][i]){ 
        IMax = Stk->prof[0][i];
      }else if(IMin > Stk->prof[0][i]){
        IMin = Stk->prof[0][i];
      } 

      dtmp = Stk->prof[3][i]/Stk->prof[0][i];
      if(VMin > dtmp){ 
        VMin = dtmp;
        iV1 = i;
      }else if(VMax < dtmp){ 
        VMax = dtmp;
        iV2 = i;
      }

      dtmp = sqrt(Stk->prof[1][i]*Stk->prof[1][i]+Stk->prof[2][i] \
        *Stk->prof[2][i])/Stk->prof[0][i]; 

      LPmean += dtmp;
      LVmean += fabs(Stk->prof[3][i]/Stk->prof[0][i]);
      if(LP < dtmp){
        LP = dtmp;
        iL = i;
      }
    }

    LPmean /= Stk->nl;

    // normalization factor
    Input->step[8] = sumI/Stk->nl*0.2;
    Stk->norm = sumI/100./Stk->nl;
    Stk->norm *= Stk->norm;
    //Stk->norm = sumI/Stk->nl;

    // Src
    Par->Par_Guess[8] = IMax;
    Par->Limits[8][0] = Par->Par_Guess[8]*0.6;
    Par->Limits[8][1] = Par->Par_Guess[8]*1.4;


    int nminimum = 0, iminimum = 0, ntmp;
    sumL = 0.0, weight = 0.0;

    for(i=1; i<Stk->nl-1; i++){
      if(Stk->prof[0][i]<Stk->prof[0][i-1] \
          && Stk->prof[0][i]<Stk->prof[0][i+1]){
        nminimum++;
        iminimum = i;
      }
    }

    if(nminimum == 1 && LV < 0.05 && LP < 0.1 && false){
      /*
      ntmp = iminimum<Stk->nl-1-iminimum?iminimum:Stk->nl-1-iminimum;
      for(i=iminimum-ntmp; i<=iminimum-ntmp; i++){
        sumL += Stk->Lambda[i]*(Par->Par_Guess[8]-Stk->prof[0][i]);
        weight += Par->Par_Guess[8]-Stk->prof[0][i];
      }
      dtmp = sumL/weight;
      */
      dtmp = Parabolic(Stk->Lambda, Stk->prof[0], iminimum);

    }else{

      for(i=0; i<Stk->nl; i++){
        sumL += Stk->Lambda[i]*(Par->Par_Guess[8]-Stk->prof[0][i]);
        weight += Par->Par_Guess[8]-Stk->prof[0][i];
      }
      dtmp = sumL/weight;
    }

    Par->Par_Guess[4] = (dtmp-Input->Lines[0].Lambda0) \
        /Input->Lines[0].Lambda0*Par_C/1e3;


    Par->Limits[4][0] = Par->Par_Guess[4]-Input->delta_v*5;
    Par->Limits[4][1] = Par->Par_Guess[4]+Input->delta_v*5;

    

    if(LV>0.001){
      if(iV1<iV2){
          Blos = -Input->VCoeffi*LV;
          if(Blos<-1300) Blos = -1300.;
        }else{
          Blos = Input->VCoeffi*LV;
          if(Blos>1300) Blos = 1300.;
      }
    }else{
      Blos = 0;
    }

    Bpos = 20.+Input->LCoeffi*LP;
    if(Bpos>1500){ 
      Bpos = 1500.;
    }else if(Bpos<150){
      Bpos = 150;
    }

    // Bstrength
    Par->Par_Guess[1] = sqrt(Blos*Blos+Bpos*Bpos);

    // inclination
    Par->Par_Guess[2] = acos(Blos/Par->Par_Guess[1]);
    //Par->Par_Guess[2] = Par_Pi/2.;

    // tan(2*alpha) = U/Q;
    Par->Par_Guess[3] = 0.5*atan2(Stk->prof[2][iL], Stk->prof[1][iL]);
    //Par->Par_Guess[3] = 0.5*atan2(sumU, sumQ);

    if(Par->Par_Guess[3]<0) Par->Par_Guess[3] += Par_Pi;

    // Damp
    Par->Par_Guess[6] = 0.5;
    
    if(Par->Par_Guess[1]>1100){
      // Dopp
      Par->Par_Guess[5] = 20;
      // Eta
      Par->Par_Guess[7] = 20;
      // Beta
      Par->Par_Guess[9] = 0.3;
    }else{
      // Dopp
      Par->Par_Guess[5] = 30;
      // Eta
      Par->Par_Guess[7] = 5;
      // Beta
      Par->Par_Guess[9] = 0.15;
    }

    Input->value_const[7] = Par->Par_Guess[7];
    
    bounds_check(Par->Par_Guess, Par);

    for(i=1; i<10; i++){
      Par->Par[i] = Par->Par_Guess[i];
    }

    for(i=0;i<Input->nline;i++){
      Input->Lines[i].Par = Par->Par;
    }

    return  0;
}

/*--------------------------------------------------------------------------------*/

extern int bounds_check(double *Parnew, STRUCT_PAR *Par){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Check the bounds.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Par, structure with parameters
      Output:
        Par, structure with parameters
    ######################################################################*/

/*--------------------------------------------------------------------------------*/  
    
    int i;

    
    if(Parnew[2]<-0.8*Par_Pi) Parnew[2] = -0.8*Par_Pi;
    if(Parnew[2]>1.8*Par_Pi) Parnew[2] = 1.8*Par_Pi;

    if(Parnew[2]<0){
      Parnew[2] = -Parnew[2];
      //Parnew[3] += Par_Pi;
    }else if(Parnew[2]>Par_Pi){
      Parnew[2] = Par_Pi*2.-Parnew[2];
      //Parnew[3] += Par_Pi;
    }

    for(i=1;i<10;i++){
      if(Parnew[i] < Par->Limits[i][0]){
        Parnew[i] = Par->Limits[i][0];
      }else if(Parnew[i] > Par->Limits[i][1]){
        Parnew[i] = Par->Limits[i][1];
      }
    }

    while(Parnew[3]<0){
      Parnew[3] += Par_Pi;
    }

    while(Parnew[3]>Par_Pi){
      Parnew[3] -= Par_Pi;
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/