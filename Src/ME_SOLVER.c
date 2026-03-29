
#include "ME_SOLVER.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        09 Mar. 2026  (Hao Li)
          --- Updates:  
              removed the support of solver for multi lines to improve. 

        18 Apr. 2025
          --- Updates: Update the initial guess (Hao Li)

        11 Apr. 2025
          --- Updates:  Update the initial guess for Fe I 15648 (Hao Li)

        28 Jun. 2024
          --- Initial commit (Hao Li)
     
        to do list:
            multi lines. filling factor.
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

// Pi Ratio of circumference to diameter
static const double L_Pi = 3.14159265358979323846;
// Light speed
static const double L_C = 299792458.0;

/*--------------------------------------------------------------------------------*/

static inline double clamp(double x, double *bounds) {

  /*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        restricts a variable to lie within a specified range.
      Record of revisions:
        09 Mar. 2026
      Input parameters:        
        x, the variable.
        bounds, the limits.
      Return:
        the claped value.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    return x<bounds[0] ? bounds[0] : (x>bounds[1] ? bounds[1] : x);
}

/*--------------------------------------------------------------------------------*/

static inline double Parabolic(const double *X, const double *Y, int Indx){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        3 points Parabolic interpolation to get the minimum position.
      Record of revisions:
        28 Jun. 2024
      Input parameters:        
        X, array with x values.
        Y, array with y values.
        Indx, the index of minimum Y.
      Return:
        The X position corresponding to the Y minimum .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    const double X_tmp1 = X[Indx] - X[Indx-1];
    const double X_tmp2 = X[Indx] - X[Indx+1];
    const double Y_tmp1 = Y[Indx] - Y[Indx-1];
    const double Y_tmp2 = Y[Indx] - Y[Indx+1];

    return X[Indx]-0.5*(X_tmp1*X_tmp1*Y_tmp2-X_tmp2*X_tmp2*Y_tmp1) \
        /(X_tmp1*Y_tmp2-X_tmp2*Y_tmp1);
}

/*--------------------------------------------------------------------------------*/

int Milne_Eddington_Single(double *modelpar, STRUCT_STK *Stk, \
    STRUCT_PARA *Para, bool Deriv){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Calculate the Stokes profiles under M-E atmosphere model (normal 
            Zeeman effect)
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        modelpar, the model parameters.
        Stk, a structure storing the Stokes profiles.
        Para, a structure with the model parameters.
        Deriv, output the derivatives of not.
      Output:
        Stk, structure with Stokes profiles
      Reference:
        LL04 5.47, LL04 9.4 9.8, del Toro Iniesta 2003
      Note:
        Par[0] Bmod; Par[1] ThetaB; Par[2] PhiB; Par[3] Vlos; Par[4] Dopp;
        Par[5] Damp; Par[6] Eta; Par[7] S0+S1; Par[8] S0/(S0+S1);
        Par[0], Par[1], Par[2], Par[3], Par[7], Par[8] must be same for 
        all lines. 
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
  
    static const double L_SqrtPi = 1.77245385090551588191;

    STRUCT_FADDEEVA *Fadd = &(Stk->Fadd);
    STRUCT_MELINE *Lines = Para->lines;

    // Atom->Lambda in A, modelpar[4] in mA, Atom->BShift in mA
    // Atom->BShift = 4.6686e-10*Atom->Geffect*Atom->Lambda*Atom->Lambda;
    const double Para4 = modelpar[4];
    const double Para5 = modelpar[5];
    const double Para6 = modelpar[6];
    const double Para6_2 = 0.5*Para6;
    const double Para7 = modelpar[7];
    const double Para8 = modelpar[8];
    // get the source function
    const double B0 = Para7*Para8;
    const double B1 = Para7-B0;
    const double BshiftNorm = Lines->BShift/Para4;
    const double LambdaNorm = Lines->Lambda0/Para4;

    // precompute some sin and cos functions
    const double SIN_THB = sin(modelpar[1]);
    const double COS_THB = cos(modelpar[1]);
    const double SIN_THB_SQ = SIN_THB*SIN_THB;
    const double COS_THB_SQ_1 = COS_THB*COS_THB+1;
    const double COS2PHI = cos(2.*modelpar[2]);
    const double SIN2PHI = sin(2.*modelpar[2]);
    const double ScSqCs2 = Para6_2*SIN_THB_SQ*COS2PHI;
    const double ScSqSc2 = Para6_2*SIN_THB_SQ*SIN2PHI;
    const double ScCsB2 = Para6_2*2.0*SIN_THB*COS_THB;
    const double CsB2 = Para6_2*COS_THB;

    const double ShiftB = modelpar[0]*BshiftNorm;
    const double ShiftV = 1e3*(-1e3*modelpar[3]/L_C*Lines->Lambda0 \
        -Lines->Lambda0)/Para4;

    const int nw = Stk->nw;
    double *synI = Stk->syn;
    double *synQ = synI+nw;
    double *synU = synQ+nw;
    double *synV = synU+nw;

    double *JacobI = Stk->Jacobian;
    double *JacobQ = JacobI+9*nw;
    double *JacobU = JacobQ+9*nw;
    double *JacobV = JacobU+9*nw;
    
    #define JacbI(ipar,iw) JacobI[ipar*nw+iw]
    #define JacbQ(ipar,iw) JacobQ[ipar*nw+iw]
    #define JacbU(ipar,iw) JacobU[ipar*nw+iw]
    #define JacbV(ipar,iw) JacobV[ipar*nw+iw]

    const double *step = Para->step;
    const double step7 = Para->step[7], step8 = Para->step[8];
    const double *Lambda = Stk->Lambda;

    //for derivatives
    double Deltap[7], tempp[7], temp1p[7], temp2p[7], temp3p[7];
    double etap[4][7], rhop[4][7], phip[3][6], psip[3][6];
    double phi[3], psi[3], eta[4], rho[4];
    double H=0, L=0;

    for(int iw=0; iw<nw; iw++){
      double ShiftW = 1e3*Lambda[iw]/Para4+ShiftV;
      // compute the phi and psi
      for(int q=-1; q<=1; q++){
        int iq = q+1;
        // q=Mu-Ml,(1 sigma bule, -1 sigma red, 0 pi)
        double Lam_Shifts = ShiftW+q*ShiftB;
        Faddeeva(Lam_Shifts, Para5, &H, &L, Fadd);

        phi[iq] = H/L_SqrtPi;
        psi[iq] = L/L_SqrtPi;
       
        if(Deriv){
          // the derivatives are devided L_SqrtPi
          double Hp_lam = 2.*(-Lam_Shifts*H+Para5*L)/L_SqrtPi;
          double Lp_lam = -2.*(-1./L_SqrtPi+Para5*H+Lam_Shifts*L)/L_SqrtPi;
          double *ptr1 = phip[iq];
          double *ptr2 = psip[iq];
          ptr1[0] = q*BshiftNorm*Hp_lam;
          ptr2[0] = q*BshiftNorm*Lp_lam;
          ptr1[3] = -1e6/L_C*LambdaNorm*Hp_lam;
          ptr2[3] = -1e6/L_C*LambdaNorm*Lp_lam;
          ptr1[4] = -Lam_Shifts/Para4*Hp_lam;
          ptr2[4] = -Lam_Shifts/Para4*Lp_lam;
          ptr1[5] = -Lp_lam; //  Hp_damp = -Lp_lam;
          ptr2[5] = Hp_lam;  //  Hp_lam = Lp_damp;                
        }
      }

      double phisum = 0.5*(phi[0]+phi[2]);
      double phitot = phi[1]-phisum;
      double psitot = psi[1]-0.5*(psi[0]+psi[2]);
      // Stokes profiles LL04 Page 414 Eq 9.110
      eta[0] = Para6_2*(phi[1]*SIN_THB_SQ+phisum*COS_THB_SQ_1);
      eta[1] = phitot*ScSqCs2;
      eta[2] = phitot*ScSqSc2;
      eta[3] = (phi[0]-phi[2])*CsB2;
      rho[1] = psitot*ScSqCs2;
      rho[2] = psitot*ScSqSc2;
      rho[3] = (psi[0]-psi[2])*CsB2;

      // some coefficients
      double temp = 1+eta[0];
      double tempsq = temp*temp;
      double temp1 = eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
      double temp2 = rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
      double temp3 = eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
      double Delta = (tempsq*(tempsq-temp1+temp2)-temp3*temp3);

      double B1D = B1/Delta;
      double C1 = (eta[3]*rho[2]-eta[2]*rho[3]);
      double C2 = (eta[1]*rho[3]-eta[3]*rho[1]);
      double C3 = (eta[2]*rho[1]-eta[1]*rho[2]); 
      synI[iw] = B0+B1D*(temp*(tempsq+temp2));
      synQ[iw] = -B1D*(tempsq*eta[1]+temp*C1+rho[1]*temp3);
      synU[iw] = -B1D*(tempsq*eta[2]+temp*C2+rho[2]*temp3);
      synV[iw] = -B1D*(tempsq*eta[3]+temp*C3+rho[3]*temp3);

      // conpute the derivatives of eta and rho with respect 
      // to some parameters
      if(Deriv){
        //der  bmod
        double *phip0 = phip[0];
        double *phip1 = phip[1];
        double *phip2 = phip[2];
        double *psip0 = psip[0];
        double *psip1 = psip[1];
        double *psip2 = psip[2];
        double phipsum = 0.5*(phip0[0]+phip2[0]);
        etap[0][0] = Para6_2*(phip1[0]*SIN_THB_SQ \
            +phipsum*COS_THB_SQ_1);
        etap[1][0] = (phip1[0]-phipsum);
        etap[2][0] = etap[1][0]*ScSqSc2;
        etap[1][0] *= ScSqCs2;
        etap[3][0] = (phip0[0]-phip2[0])*CsB2;
        rhop[1][0] = (psip1[0]-0.5*(psip0[0]+psip2[0]));
        rhop[2][0] = rhop[1][0]*ScSqSc2;
        rhop[1][0] *= ScSqCs2;
        rhop[3][0] = (psip0[0]-psip2[0])*CsB2;

        //der thetab
        etap[0][1] = phitot*ScCsB2;
        etap[1][1] = etap[0][1]*COS2PHI;
        etap[2][1] = etap[0][1]*SIN2PHI;
        etap[3][1] = Para6_2*(phi[0]-phi[2])*(-SIN_THB);
        rhop[1][1] = psitot*ScCsB2*COS2PHI;
        rhop[2][1] = psitot*ScCsB2*SIN2PHI;
        rhop[3][1] = Para6_2*(psi[0]-psi[2])*(-SIN_THB);

        //der phib
        etap[0][2] = 0;
        etap[1][2] = -2.*eta[2];
        etap[2][2] = 2.*eta[1];
        etap[3][2] = 0; 
        rhop[1][2] = -2.*rho[2];
        rhop[2][2] = 2.*rho[1];
        rhop[3][2] = 0;

        //der to velocit  doppler width  dampping parameter
        for(int ipar=3; ipar<6; ipar++){
          phipsum = 0.5*(phip0[ipar]+phip2[ipar]);
          etap[0][ipar] = Para6_2*(phip1[ipar]*SIN_THB_SQ \
              +phipsum*COS_THB_SQ_1);
          etap[1][ipar] = phip1[ipar]-phipsum;
          etap[2][ipar] = etap[1][ipar]*ScSqSc2;
          etap[1][ipar] *= ScSqCs2;
          etap[3][ipar] = (phip0[ipar]-phip2[ipar])*CsB2;
          rhop[1][ipar] = psip1[ipar]-0.5*(psip0[ipar]+psip2[ipar]);
          rhop[2][ipar] = rhop[1][ipar]*ScSqSc2;
          rhop[1][ipar] *= ScSqCs2;
          rhop[3][ipar] = (psip0[ipar]-psip2[ipar])*CsB2;
        }

        etap[0][6] = eta[0]/Para6;
        for(int iStk=1; iStk<4; iStk++){
          etap[iStk][6] = eta[iStk]/Para6;
          rhop[iStk][6] = rho[iStk]/Para6;
        }

        for(int ipar=0; ipar<7; ipar++){
          tempp[ipar] = etap[0][ipar];
          temp1p[ipar] = 2.0*eta[1]*etap[1][ipar]+2.0*eta[2]*etap[2][ipar] \
              +2.0*eta[3]*etap[3][ipar];
          temp2p[ipar] = 2.0*rho[1]*rhop[1][ipar]+2.0*rho[2]*rhop[2][ipar] \
              +2.0*rho[3]*rhop[3][ipar];
          temp3p[ipar] = etap[1][ipar]*rho[1]+eta[1]*rhop[1][ipar] \
              +etap[2][ipar]*rho[2]+eta[2]*rhop[2][ipar]+etap[3][ipar]*rho[3] \
              +eta[3]*rhop[3][ipar];
          Deltap[ipar] = 2*temp*tempp[ipar]*(tempsq-temp1+temp2) \
              +tempsq*(2*temp*tempp[ipar]-temp1p[ipar]+temp2p[ipar]) \
              -2.0*temp3*temp3p[ipar];
  
          double B1DD = B1D/Delta;;
          JacbI(ipar,iw) = (B1DD*(Delta*(tempp[ipar]*(tempsq+temp2) \
              +temp*(2*temp*tempp[ipar]+temp2p[ipar])) \
              -Deltap[ipar]*(temp*(tempsq+temp2))))*step[ipar];
          JacbQ(ipar,iw) = (-B1DD*(Delta \
              *(2*temp*tempp[ipar]*eta[1]+tempsq*etap[1][ipar] \
              +tempp[ipar]*C1 \
              +temp*(etap[3][ipar]*rho[2]+eta[3]*rhop[2][ipar] \
              -etap[2][ipar]*rho[3]-eta[2]*rhop[3][ipar]) \
              +rhop[1][ipar]*temp3+rho[1]*temp3p[ipar]) \
              -Deltap[ipar]*(tempsq*eta[1]+temp*(eta[3]*rho[2] \
              -eta[2]*rho[3])+rho[1]*temp3)))*step[ipar];
          JacbU(ipar,iw) = (-B1DD*(Delta \
              *(2*temp*tempp[ipar]*eta[2]+tempsq*etap[2][ipar] \
              +tempp[ipar]*C2 \
              +temp*(etap[1][ipar]*rho[3]+eta[1]*rhop[3][ipar] \
              -etap[3][ipar]*rho[1]-eta[3]*rhop[1][ipar]) \
              +rhop[2][ipar]*temp3+rho[2]*temp3p[ipar]) \
              -Deltap[ipar]*(tempsq*eta[2]+temp*(eta[1]*rho[3] \
              -eta[3]*rho[1])+rho[2]*temp3)))*step[ipar];
          JacbV(ipar,iw) = (-B1DD*(Delta \
              *(2*temp*tempp[ipar]*eta[3]+tempsq*etap[3][ipar] \
              +tempp[ipar]*C3 \
              +temp*(etap[2][ipar]*rho[1]+eta[2]*rhop[1][ipar] \
              -etap[1][ipar]*rho[2]-eta[1]*rhop[2][ipar]) \
              +rhop[3][ipar]*temp3+rho[3]*temp3p[ipar]) \
              -Deltap[ipar]*(tempsq*eta[3]+temp*(eta[2]*rho[1] \
              -eta[1]*rho[2])+rho[3]*temp3)))*step[ipar];
        }

        double T = (1-Para8)/Delta;
        JacbI(7,iw) = (Para8+T*(temp*(tempsq+temp2)))*step7;
        JacbQ(7,iw) = (-T*(tempsq*eta[1]+temp*C1+rho[1]*temp3))*step7;
        JacbU(7,iw) = (-T*(tempsq*eta[2]+temp*C2+rho[2]*temp3))*step7;
        JacbV(7,iw) = (-T*(tempsq*eta[3]+temp*C3+rho[3]*temp3))*step7;

        T = Para7/Delta;
        JacbI(8,iw) = (Para7-T*(temp*(tempsq+temp2)))*step8;
        JacbQ(8,iw) = (T*(tempsq*eta[1]+temp*C1+rho[1]*temp3))*step8;
        JacbU(8,iw) = (T*(tempsq*eta[2]+temp*C2+rho[2]*temp3))*step8;
        JacbV(8,iw) = (T*(tempsq*eta[3]+temp*C3+rho[3]*temp3))*step8;
      }
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int bounds_check(double *Parnew, STRUCT_PARA *Para){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Check the bounds.
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        Parnew, the parameters.
        Para, a structure storing the limits on the model parameters
      Output:
        Parnew, the parameters.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/  
    
    static const double L_2Pi = 2.0*L_Pi;

    Parnew[0] = clamp(Parnew[0], Para->Limits[0]);

    double p1 = Parnew[1];
    p1 = fmod(p1, L_2Pi);
    if(p1<0) p1 += L_2Pi;
    if(p1>L_Pi) p1 = L_2Pi-p1;
    Parnew[1] = p1;

    for(int ii=3; ii<9; ii++){
      Parnew[ii] = clamp(Parnew[ii], Para->Limits[ii]);
    }

    double p2 = Parnew[2];
    p2 = fmod(p2, L_Pi);
    if (p2 < 0.0) p2 += L_Pi;
    Parnew[2] = p2;

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Noise_Init(STRUCT_STK *Stk){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Initialize the noise if noise is nor provided.
      Record of revisions:
        09 Mar. 2026.
      Input parameters:
        Stk, a structure storing the Stokes profiles.
      Output parameter:
        Stk, a structure storing the Stokes profiles.
      To do:
        Only wavelength independent noise is adopted.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    double *ptr = Stk->ivnoise;
    const int nw = Stk->nw;
    double ivnoi; 
    const double noisesq = 0.25/Stk->norm0;
    
    for(int istk=0; istk<4; istk++){
      ivnoi = Stk->Weights_SQR[istk]*noisesq;
      for(int iw=0; iw<nw; iw++){
        ptr[iw] = ivnoi;
      }
      ptr += nw;
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

int Init_Guess(STRUCT_STK *Stk, STRUCT_PARA *Para){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Initialize the model patameters
      Record of revisions:
        9 Mar. 2026.
      Input parameters:
        Stk, a structure storing the Stokes profiles.
        Para, a structure storing the model parameters.
      Output:
        Para, a structure storing the initial guess.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/  

    const int nw = Stk->nw;
    const double *profI = Stk->prof;
    const double *profQ = profI+nw;
    const double *profU = profQ+nw;
    const double *profV = profU+nw;
    const double *Lambda = Stk->Lambda;

    double Isum = 0., Qsum = 0., Usum = 0., Vsum = 0., Lsum = 0., weight;
    double Imax = 0, Lmax = 0., Vmax = 0., Vmin = 0., Imean, Lmean, Vmean;
    int iI = 0, iVmax = 0, iVmin = 0, iL = 0;
    double tmp, tmp1, tmp2, Blos = 0., Bpos = 0., thresholdV = 0., thresholdL = 0.;
    double Lp[nw];

    for(int iw=0; iw<Stk->nw; iw++){

      Isum += profI[iw];
      if(Imax<profI[iw]){
        iI = iw;
        Imax = profI[iw];
      }
   
      Qsum += fabs(profQ[iw]);
      Usum += fabs(profU[iw]);
      Lp[iw] = sqrt(profQ[iw]*profQ[iw]+profU[iw]*profU[iw]);
      if(Lmax<Lp[iw]){
        iL = iw;
        Lmax = Lp[iw];
      }
      Lsum += Lp[iw];

      Vsum += fabs(profV[iw]);
      if(Vmin>profV[iw]){
        iVmin = iw;
        Vmin = profV[iw];
      }else if(Vmax<profV[iw]){
        iVmax = iw;
        Vmax = profV[iw];
      }
    }

    Imean = Isum/nw;
    if(Imean<Stk->Icriteria) return 1;

    // normalization factor
    Stk->norm0 = 1e-2*Isum/Stk->nw;
    Stk->norm0 *= Stk->norm0;
    Stk->normp = 1+Qsum*Qsum/(Isum*Isum)+Usum*Usum/(Isum*Isum) \
        +Vsum*Vsum/(Isum*Isum);

    Lmax = 100*Lmax/Imean;
    Vmin = 100*Vmin/Imean;
    Vmax = 100*Vmax/Imean;
    Lmean = 100*Lsum/Imean;
    Vmean = 100*Vsum/Imean;
   
    if(Para->lines->Lambda0>15640&&Para->lines->Lambda0<15660){
      thresholdV = .8;
      thresholdL = 1.2;
    }else{
      thresholdV = 3.0;
      thresholdL = 3.0;
    }

    double LV = 0.5*(Vmax-Vmin);  
    if(LV >0.1){
      Blos = Para->VCoeffi*LV;
      if(Blos>1300) Blos = 1300.;
      if(iVmin<iVmax) Blos = -Blos;
    }

    Bpos = 20+Para->LCoeffi*Lmax;
    Bpos = Bpos<1500 ? Bpos : 1500;
    if(Bpos<150) Bpos = 150;

    tmp1 = 0;
    Para->Par_Guess[3] = 0;
    if(LV>thresholdV || Lmax>thresholdL){

      tmp1 = 0.5*(Lambda[iVmin]+Lambda[iVmax])-Para->lines->Lambda0;

      Lsum = 0;
      weight = 0.;
      for(int iw =0; iw<nw; iw++){
        if(Lp[iw]>0.75*thresholdL*Imean){
          Lsum += Lambda[iw]*Lp[iw];
          weight += Lp[iw];
        }
      }
      tmp2 = weight>0? Lsum/weight-Para->lines->Lambda0:0;

      if(LV>thresholdV && Lmax>thresholdL){
        if(LV>4*thresholdV){
          tmp = tmp1;
        }else if(LV<2*Lmax && Lmax<2*LV){
          tmp = 0.5*(tmp1+tmp2);
        }else{
          tmp = tmp2;
        }

      }else if (LV>thresholdV){
        tmp = tmp1;
      }else{
        tmp = tmp2;
      }

      Para->Par_Guess[3] = tmp/Para->lines->Lambda0*L_C/1e3;
      if(Para->Par_Guess[3]>Para->delta_v*4){
        Para->Par_Guess[3]=Para->delta_v*4;
      }else if(Para->Par_Guess[3]<-Para->delta_v*4){
        Para->Par_Guess[3]=-Para->delta_v*4;
      }
    }else{ 
      if(iL>=0 && iL<nw-2){
        if(profI[iL-1]<profI[iL-2] \
            && profI[iL]<profI[iL-1] \
            && profI[iL]<profI[iL+1] \
            && profI[iL+1]<profI[iL+2]){
          tmp = Lambda[iL]-Para->lines->Lambda0;
        }
      }
      Para->Par_Guess[3] = tmp/Para->lines->Lambda0*L_C/1e3;
      if(Para->Par_Guess[3]>Para->delta_v){
        Para->Par_Guess[3]=Para->delta_v;
      }else if(Para->Par_Guess[3]<-Para->delta_v){
        Para->Par_Guess[3]=-Para->delta_v;
      }
    }

    // Bstrength
    Para->Par_Guess[0] = sqrt(Blos*Blos+Bpos*Bpos);

    // Inclination
    Para->Par_Guess[1] = acos(Blos/Para->Par_Guess[0]);
    
    // tan(2*alpha) = U/Q;
    weight = 0.;
    double phi = 0.;
    double tmpphi;
    for(int iw =0; iw<nw; iw++){
      if(Lp[iw] > thresholdL){
        tmpphi = 0.5*atan2(profU[iw], profQ[iw]);
        tmpphi = tmpphi>0? tmpphi : L_Pi+tmpphi;
        phi += tmpphi*Lp[iw];
        weight += Lp[iw];
      }
    }
    if(weight > 2.*thresholdL){
      Para->Par_Guess[2] = phi/weight;
    }else{
      Para->Par_Guess[2] = L_Pi/2.;
    }



    if(Para->lines->Lambda0>15640&&Para->lines->Lambda0<15660){
      // Doppler width
      Para->Par_Guess[4] = 60.;
      // Eta
      Para->Par_Guess[6] = 10;
      // Beta
      Para->Par_Guess[8] = 0.7;
    }else{
      if(Para->Par_Guess[0]>1100){
        // Doppler width
        Para->Par_Guess[4] = 20;
        // Eta
        Para->Par_Guess[6] = 20;
        // Beta
        Para->Par_Guess[8] = 0.3;
      }else{
        // Doppler width
        Para->Par_Guess[4] = 30;
        // Eta
        Para->Par_Guess[6] = 5;
        // Beta
        Para->Par_Guess[8] = 0.15;
      }
    }
    //Para->Par_Guess[3] = 0.;

    // Damp
    Para->Par_Guess[5] = 0.5;

    // Src
    Para->Par_Guess[7] = Imax;
    Para->Limits[7][0] = Imax*0.6;
    Para->Limits[7][1] = Imax*1.4;
    Para->step[7] = Imax*0.1;

    //Input->value_const[6] = Para->Par_Guess[6];
    
    bounds_check(Para->Par_Guess, Para);

    for(int iw=0; iw<9; iw++){
      Para->Par[iw] = Para->Par_Guess[iw];
    }

    return  0;
}

/*--------------------------------------------------------------------------------*/
