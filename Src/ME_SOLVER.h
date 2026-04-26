
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include "ALLOCATION.h"
#include "FADDEEVA.h"
#include "RINPUT.h"

/*--------------------------------------------------------------------------------*/

typedef struct Struct_MEline{

    // center of the wavelength, effective lande factor, a precomputed 
    // coeffcient related to the Zeeman splitting.
    double Lambda0, Geff, BShift;

    // Milne-Eddington model parameters Par[0:8]
    //double Par[9];

}STRUCT_MELINE;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Parameter{

    STRUCT_MELINE *lines;
    
    int nline; 
    // 9 model parameters, their best in the inversion, and the guess
    // Bmod, ThetaB, PhiB, V, Dopp, Damp, Eta, Src, Beta;
    double *Par, *Par_Best, *Par_tmp, Par_Guess[9];

    double *step;
    // limits on the parameter
    double Limits[9][2];

    // iteration number;
    int Niter;

    // chisq and best chisq
    double Chisq, Chisq_Best, delta_v;

    bool *inv, *regl;
    double *value_const, *Regul_weight;
    // number of free parameter.
    int npar;

    double VCoeffi, LCoeffi;

}STRUCT_PARA;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Stokes{

    // number of spectral lines and wavelength points;
    int nw, ncut, ic, i0, i1, ni;
    // nomalization factor, the weights for Stokes parameters and 
    // the squre of the weights
    double norm0, normp, Weights[4], Weights_SQR[4];
    // the wavelength, input profiles, synthesized profiles, fitting, 
    // best the fit, noise, and the Jacobian used for the inversion.
    double *Lambda, *prof, *syn, *fit, *fit_best, *ivnoise, *Jacobian;

    double Icriteria;

    STRUCT_FADDEEVA Fadd;

}STRUCT_STK;

/*--------------------------------------------------------------------------------*/

extern int Milne_Eddington_Single(double *modelpar, STRUCT_STK *Stk, \
    STRUCT_PARA *Para, bool Deriv);
  
extern int bounds_check(double *Parnew, STRUCT_PARA *Para);

extern int Noise_Init(STRUCT_STK *Stk);

extern int Init_Guess(STRUCT_STK *Stk, STRUCT_PARA *Para);

/*--------------------------------------------------------------------------------*/

/*
extern int Milne_Eddington(STRUCT_MELINE *Lines, int nline, \
        STRUCT_STK *Stk, STRUCT_FADDEEVA *Fadd, STRUCT_INPUT *Input, 
        bool Deriv){
  

    #define LOCAL_SqrtPi  1.77245385090551588191


    // precompute some sin and cos functions
    const double SIN_THB = sin(Lines[0].Par[1]);
    const double COS_THB = cos(Lines[0].Par[1]);
    const double SIN_THB_SQ = SIN_THB*SIN_THB;
    const double COS_THB_SQ = COS_THB*COS_THB;
    const double COS2PHI = cos(2.*Lines[0].Par[2]);
    const double SIN2PHI = sin(2.*Lines[0].Par[2]);

    // memory allocation
    double *phi = (double *)VECTOR(-1, 1, enum_dbl, false);
    double *psi = (double *)VECTOR(-1, 1, enum_dbl, false);

    double *eta = (double *)alloca(4*sizeof(double));
    double *rho = (double *)alloca(4*sizeof(double));
    double *ShiftB = (double *)alloca(nline*sizeof(double));
    double *ShiftV = (double *)alloca(nline*sizeof(double));

    // Atom->Lambda in A, Lines[0].Par[4] in mA, Atom->BShift in mA
    // Atom->BShift = 4.6686e-10*Atom->Geffect*Atom->Lambda*Atom->Lambda;

    for(int il=0; il<nline; il++){
      ShiftB[il] = Lines[il].Par[0]*Lines[il].BShift/Lines[il].Par[4];
      ShiftV[il] = 1e3*(-1e3*Lines[il].Par[3]/Par_C*Lines[il].Lambda0 \
        -Lines[il].Lambda0)/Lines[il].Par[4];
    }

    const double Para7 = Lines[0].Par[7];
    const double Para8 = Lines[0].Par[8];

    // get the source function
    double B0 = Para7*Para8;
    double B1 = Para7-B0;
    
    double Delta = 0, temp, temp1, temp2, temp3;
    double Lam_Shifts, ShiftW, H, L, Hp_lam, Hp_damp, Lp_lam, Lp_damp;

    double *Deltap = NULL, *tempp = NULL, *temp1p = NULL, *temp2p = NULL;
    double *temp3p = NULL, **etap = NULL, **rhop = NULL, **phip = NULL;
    double **psip = NULL;

    if(Deriv){
      phip = (double **)MATRIX(-1, 1, 0, 5, enum_dbl, false);
      psip = (double **)MATRIX(-1, 1, 0, 5, enum_dbl, false);
      etap = (double **)MATRIX(0, 3, 0, 6, enum_dbl, false);
      rhop = (double **)MATRIX(1, 3, 0, 6, enum_dbl, false);
      Deltap = (double *)alloca(7*sizeof(double))
      tempp = (double *)alloca(7*sizeof(double))
      temp1p = (double *)alloca(7*sizeof(double))
      temp2p = (double *)alloca(7*sizeof(double))
      temp3p = (double *)alloca(7*sizeof(double))
    }

    for(int iw=0; iw<Stk->nl; iw++){

      for(int il=0; il<nline; il++){

        const double Para4 = Lines[il].Par[4];
        const double Para5 = Lines[il].Par[5];
        ShiftW = 1e3*Stk->Lambda[iw]/Para4+ShiftV[il];

        // compute the phi and psi
        for (int iq=-1; iq<=1; iq++){
          // q=Mu-Ml,(1 sigma bule, -1 sigma red, 0 pi)
          Lam_Shifts = ShiftW+q*ShiftB[il];
          Faddeeva(Lam_Shifts, Para5, &H, &L, Fadd);

          if(il==0){
            phi[iq] = H/LOCAL_SqrtPi;
            psi[iq] = L/LOCAL_SqrtPi;
          }else{
            phi[iq] += H/LOCAL_SqrtPi;
            psi[iq] += L/LOCAL_SqrtPi;
          }

          if(Deriv){
            Hp_lam = 2*(-Lam_Shifts*H+Para5*L);
            Lp_damp = Hp_lam;
            Hp_damp = 2*(-1/LOCAL_SqrtPi+Para5*H+Lam_Shifts*L);
            Lp_lam = -Hp_damp;
            if(il==0){
              phip[iq][0] = q*Lines[il].BShift/Para4 \
                  *Hp_lam/LOCAL_SqrtPi;
              psip[iq][0] = q*Lines[il].BShift/Para4 \
                  *Lp_lam/LOCAL_SqrtPi;
              phip[iq][3] = -1e6/Par_C*Lines[il].Lambda0/Para4 \
                  *Hp_lam/LOCAL_SqrtPi;
              psip[iq][3] = -1e6/Par_C*Lines[il].Lambda0/Para4 \
                  *Lp_lam/LOCAL_SqrtPi;
              phip[iq][4] = -Lam_Shifts/Para4*Hp_lam/LOCAL_SqrtPi;
              psip[iq][4] = -Lam_Shifts/Para4*Lp_lam/LOCAL_SqrtPi;
              phip[iq][5] = Hp_damp/LOCAL_SqrtPi;
              psip[iq][5] = Lp_damp/LOCAL_SqrtPi;
            }else{
              phip[iq][0] += q*Lines[il].BShift/Para4 \
                  *Hp_lam/LOCAL_SqrtPi;
              psip[iq][0] += q*Lines[il].BShift/Para4 \
                  *Lp_lam/LOCAL_SqrtPi;
              phip[iq][3] += -1e6/Par_C*Lines[il].Lambda0/Para4 \
                  *Hp_lam/LOCAL_SqrtPi;
              psip[iq][3] += -1e6/Par_C*Lines[il].Lambda0/Para4 \
                  *Lp_lam/LOCAL_SqrtPi;
              phip[iq][4] += -Lam_Shifts/Para4*Hp_lam/LOCAL_SqrtPi;
              psip[iq][4] += -Lam_Shifts/Para4*Lp_lam/LOCAL_SqrtPi;
              phip[iq][5] += Hp_damp/LOCAL_SqrtPi;
              psip[iq][5] += Lp_damp/LOCAL_SqrtPi;
            }       
          }
        }
      }

      const double Para6 = Lines[0].Par[6];
      const double Para6_2 = 0.5*Para6;
      const double phisum = 0.5*(phi[-1]+phi[1]);
      const double psisum = 0.5*(psi[-1]+psi[1]);
      // Stokes profiles LL04 Page 414 Eq 9.110
      eta[0] = Para6_2*(phi[0]*SIN_THB_SQ+phisum*(1+COS_THB_SQ));
      eta[1] = Para6_2*(phi[0]-phisum)*SIN_THB_SQ*COS2PHI;
      eta[2] = Para6_2*(phi[0]-phisum)*SIN_THB_SQ*SIN2PHI;
      eta[3] = Para6_2*(phi[-1]-phi[1])*COS_THB;
      rho[1] = Para6_2*(psi[0]-psisum)*SIN_THB_SQ*COS2PHI;
      rho[2] = Para6_2*(psi[0]-psisum)*SIN_THB_SQ*SIN2PHI;
      rho[3] = Para6_2*(psi[-1]-psi[1])*COS_THB;

      // some coefficients
      temp = 1+eta[0];
      temp1 = eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
      temp2 = rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
      temp3 = eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
      Delta = (temp*temp*(temp*temp-temp1+temp2)-temp3*temp3);

      Stk->syn[0][iw] = B0+B1/Delta*(temp*(temp*temp+temp2));
      Stk->syn[1][iw] = -B1/Delta*(temp*temp*eta[1] \
          +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3);
      Stk->syn[2][iw] = -B1/Delta*(temp*temp*eta[2] \
          +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3);
      Stk->syn[3][iw] = -B1/Delta*(temp*temp*eta[3] \
          +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3);

      // conpute the derivatives of eta and rho with respect 
      // to some parameters
      if(Deriv){
        //der  bmod
        etap[0][0] = Para6_2*(phip[0][0]*SIN_THB_SQ \
            +0.5*(phip[-1][0]+phip[1][0])*(1+COS_THB_SQ));
        etap[1][0] = Para6_2*(phip[0][0] \
            -0.5*(phip[-1][0]+phip[1][0]))*SIN_THB_SQ*COS2PHI;
        etap[2][0] = Para6_2*(phip[0][0] \
            -0.5*(phip[-1][0]+phip[1][0]))*SIN_THB_SQ*SIN2PHI;
        etap[3][0] = Para6_2*(phip[-1][0]-phip[1][0])*COS_THB;
        rhop[1][0] = Para6_2*(psip[0][0] \
            -0.5*(psip[-1][0]+psip[1][0]))*SIN_THB_SQ*COS2PHI;
        rhop[2][0] = Para6_2*(psip[0][0] \
            -0.5*(psip[-1][0]+psip[1][0]))*SIN_THB_SQ*SIN2PHI;
        rhop[3][0] = Para6_2*(psip[-1][0]-psip[1][0])*COS_THB;

        //der thetab
        etap[0][1] = Para6_2*(phi[0]*(2.0*SIN_THB*COS_THB) \
            +phisum*(-2.0*COS_THB*SIN_THB));
        etap[1][1] = Para6_2*(phi[0]-phisum) \
            *(2.0*SIN_THB*COS_THB)*COS2PHI;
        etap[2][1] = Para6_2*(phi[0]-phisum) \
            *(2.0*SIN_THB*COS_THB)*SIN2PHI;
        etap[3][1] = Para6_2*(phi[-1]-phi[1])*(-SIN_THB);
        rhop[1][1] = Para6_2*(psi[0]-psisum) \
            *(2.0*SIN_THB*COS_THB)*COS2PHI;
        rhop[2][1] = Para6_2*(psi[0]-psisum) \
            *(2.0*SIN_THB*COS_THB)*SIN2PHI;
        rhop[3][1] = Para6_2*(psi[-1]-psi[1])*(-SIN_THB);

        //der phib
        etap[0][2] = 0;
        etap[1][2] = Para6_2*(phi[0]-phisum)*SIN_THB_SQ*(-2.0*SIN2PHI);
        etap[2][2] = Para6_2*(phi[0]-phisum)*SIN_THB_SQ*(2.0*COS2PHI);
        etap[3][2] = 0; 
        rhop[1][2] = Para6_2*(psi[0]-psisum)*SIN_THB_SQ*(-2.0*SIN2PHI);
        rhop[2][2] = Para6_2*(psi[0]-psisum)*SIN_THB_SQ*(2.0*COS2PHI);
        rhop[3][2] = 0;

        //der to velocit  doppler width  dampping parameter
        for(int ipar=3; ipar<6; ipar++){
          etap[0][ipar] = Para6_2*(phip[0][ipar]*SIN_THB_SQ \
              +0.5*(phip[-1][ipar]+phip[1][ipar])*(1+COS_THB_SQ));
          etap[1][ipar] = Para6_2*(phip[0][ipar] \
              -0.5*(phip[-1][ipar]+phip[1][ipar]))*SIN_THB_SQ*COS2PHI;
          etap[2][ipar] = Para6_2*(phip[0][ipar] \
              -0.5*(phip[-1][ipar]+phip[1][ipar]))*SIN_THB_SQ*SIN2PHI;
          etap[3][ipar] = Para6_2*(phip[-1][ipar]-phip[1][ipar])*COS_THB;
          rhop[1][ipar] = Para6_2*(psip[0][ipar] \
              -0.5*(psip[-1][ipar]+psip[1][ipar]))*SIN_THB_SQ*COS2PHI;
          rhop[2][ipar] = Para6_2*(psip[0][ipar] \
              -0.5*(psip[-1][ipar]+psip[1][ipar]))*SIN_THB_SQ*SIN2PHI;
          rhop[3][ipar] = Para6_2*(psip[-1][ipar]-psip[1][ipar])*COS_THB;
        }

        etap[0][6] = eta[0]/Para6;
        for(int iStk=1; iStk < 4; iStk++){
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
          Deltap[ipar] = 2*temp*tempp[ipar]*(temp*temp-temp1+temp2) \
              +temp*temp*(2*temp*tempp[ipar]-temp1p[ipar]+temp2p[ipar]) \
              -2.0*temp3*temp3p[ipar];
  
          Stk->Jacobian[0][ipar][iw] = (B1/Delta/Delta*(Delta \
              *(tempp[ipar]*(temp*temp+temp2) \
              +temp*(2*temp*tempp[ipar]+temp2p[ipar])) \
              -Deltap[ipar]*(temp*(temp*temp+temp2))))*Input->step[ipar];
          Stk->Jacobian[1][ipar][iw] = (-B1/Delta/Delta*(Delta \
              *(2*temp*tempp[ipar]*eta[1]+temp*temp*etap[1][ipar] \
              +tempp[ipar]*(eta[3]*rho[2]-eta[2]*rho[3]) \
              +temp*(etap[3][ipar]*rho[2]+eta[3]*rhop[2][ipar] \
              -etap[2][ipar]*rho[3]-eta[2]*rhop[3][ipar]) \
              +rhop[1][ipar]*temp3+rho[1]*temp3p[ipar]) \
              -Deltap[ipar]*(temp*temp*eta[1]+temp*(eta[3]*rho[2] \
              -eta[2]*rho[3])+rho[1]*temp3)))*Input->step[ipar];
          Stk->Jacobian[2][ipar][iw] = (-B1/Delta/Delta*(Delta \
              *(2*temp*tempp[ipar]*eta[2]+temp*temp*etap[2][ipar] \
              +tempp[ipar]*(eta[1]*rho[3]-eta[3]*rho[1]) \
              +temp*(etap[1][ipar]*rho[3]+eta[1]*rhop[3][ipar] \
              -etap[3][ipar]*rho[1]-eta[3]*rhop[1][ipar]) \
              +rhop[2][ipar]*temp3+rho[2]*temp3p[ipar]) \
              -Deltap[ipar]*(temp*temp*eta[2]+temp*(eta[1]*rho[3] \
              -eta[3]*rho[1])+rho[2]*temp3)))*Input->step[ipar];
          Stk->Jacobian[3][ipar][iw] = (-B1/Delta/Delta*(Delta \
              *(2*temp*tempp[ipar]*eta[3]+temp*temp*etap[3][ipar] \
              +tempp[ipar]*(eta[2]*rho[1]-eta[1]*rho[2]) \
              +temp*(etap[2][ipar]*rho[1]+eta[2]*rhop[1][ipar] \
              -etap[1][ipar]*rho[2]-eta[1]*rhop[2][ipar]) \
              +rhop[3][ipar]*temp3+rho[3]*temp3p[ipar]) \
              -Deltap[ipar]*(temp*temp*eta[3]+temp*(eta[2]*rho[1] \
              -eta[1]*rho[2])+rho[3]*temp3)))*Input->step[ipar];
        }

        Stk->Jacobian[0][7][iw] = (Para8+(1-Para8)/Delta \
            *(temp*(temp*temp+temp2)))*Input->step[7];
        Stk->Jacobian[1][7][iw] = (-(1-Para8)/Delta*(temp*temp*eta[1] \
            +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3))*Input->step[7];
        Stk->Jacobian[2][7][iw] = (-(1-Para8)/Delta*(temp*temp*eta[2] \
            +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3))*Input->step[7];
        Stk->Jacobian[3][7][iw] = (-(1-Para8)/Delta*(temp*temp*eta[3] \
            +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3))*Input->step[7];

        Stk->Jacobian[0][8][iw] = ((Para7-Para7/Delta \
            *(temp*(temp*temp+temp2))))*Input->step[8];
        Stk->Jacobian[1][8][iw] = (Para7/Delta*(temp*temp*eta[1] \
            +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3))*Input->step[8];
        Stk->Jacobian[2][8][iw] = (Para7/Delta*(temp*temp*eta[2] \
            +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3))*Input->step[8];
        Stk->Jacobian[3][8][iw] = (Para7/Delta*(temp*temp*eta[3] \
            +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3))*Input->step[8];
       
      }
    }
    
    // free the memory
    if(Deriv){
      FREE_MATRIX((void *)phip, -1, 1, enum_dbl);
      FREE_MATRIX((void *)psip, -1, 1, enum_dbl);
      FREE_MATRIX((void *)etap, 0, 1, enum_dbl);
      FREE_MATRIX((void *)rhop, 1, 1, enum_dbl);

    }

    FREE_VECTOR(phi, -1, enum_dbl);
    FREE_VECTOR(psi, -1, enum_dbl);

    return 0;
}
*/
/*--------------------------------------------------------------------------------*/


