
#include "FADDEEVA.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        9 Jun. 2024
          --- update: add a structure STRUCT_FADDEEVA with some precomputed 
              coefficients.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

static double DAWSON(double x);

static double Fadd_erfcx(double a);

static int Hui_p6(double Nu, double y, double *H, double *L);

static int Hum_W4(double Nu, double y, double *H, double *L);

static int Faddeeva_4digits(double Nu, double y, double *H, double *L);

static int Faddeeva_5digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd);

static int Faddeeva_6digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd);

static int Faddeeva_7digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd);

static int Faddeeva_8digits(double Nu, double y, double *H, double *L, \
       STRUCT_FADDEEVA *Fadd);

/*--------------------------------------------------------------------------------*/

static double DAWSON(double x){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Dawson’s Integral
      Record of revisions:
        7 May 2018 (Hao Li)
      Input parameters:
        x, the real number.
      Output parameters:
        ans,the integral.
      Reference:
        Numerical recipes in C 2ed.
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double A1 = 2.0/3.0, A2 = 0.4, A3 = 2.0/7.0, H = 0.4;
    int NMAX = 6;
    int i, n0;
    double d1, d2, e1, e2, sum, x2, xp, xx, ans, sqrarg, TMP;
    double c[NMAX+1];
    
    xx = fabs(x);
    
    if(xx<0.2){  
      x2 = x*x; 
      ans = x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));    
    }else {
      for(i=1;i<=NMAX;i++){ 
        TMP=(2.0*i-1.0)*H;    
        sqrarg=(sqrarg=(TMP)) == 0.0 ? 0.0 : sqrarg*sqrarg;    
        c[i]=exp(-sqrarg);    
      }
        
      n0 = 2*(int)(0.5*xx/H+0.5);
      xp = xx-n0*H;
      e1 = exp(2.0*xp*H);
      e2 = e1*e1;  
      d1 = n0+1;
      d2 = d1-2.0;  
      sum = 0.0;
        
      for(i=1; i<=NMAX; i++,d1+=2.0,d2-=2.0,e1*=e2) \
          sum += c[i]*(e1/d1+1.0/(d2*e1));
      TMP = (x) >= 0.0 ? fabs(exp(-xp*xp)) : -fabs(exp(-xp*xp));
      ans = 0.5641895835*TMP*sum;
    }
    
    return ans;
}

/*--------------------------------------------------------------------------------*/

static double Fadd_erfcx(double a){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Scaled complementary error function
      Record of revisions::
        8 Jun 2024 (Hao Li)
      Input parameters:
        a, the real number.
      Output parameters:
        erfcx, the result.
      Reference:
        aghloul 2007 MNRAS.
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(a>26.6){

      double asqr = a*a;

      return ((((((162.421875/asqr-29.53125)/asqr+6.5625)/asqr-1.875) \
          /asqr+0.75)/asqr-0.5)/asqr+1)/Par_SqrtPi/a;

    }else{

      return exp(a*a)*erfc(a);

    }
    
}

/*--------------------------------------------------------------------------------*/

static int Hui_p6(double Nu, double y, double *H, double *L){

/*--------------------------------------------------------------------------------*/
   
    /*######################################################################
      Purpose:
        Voigt function.
      Record of revisions:
        07 May 2023 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Hui 1978 Journal of Quantitative Spectroscopy and Radiative Transfer
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double A[7] = {122.607931777104326, 214.382388694706425, \
        181.928533092181549, 93.155580458138441, 30.180142196210589, \
        5.912626209773153, 0.564189583562615};
    double B[7] = {122.60793177387535, 352.730625110963558, \
        457.334478783897737, 348.703917719495792, 170.354001821091472, \
        53.992906912940207, 10.479857114260399};
    
    double complex Z = y-Nu*I;
    double complex sum1=A[6], sum2=B[6]+Z;
    double complex tmp;
    int i;

    for(i=5; i>-1; i--){
      sum1 = sum1*Z+A[i];
      sum2 = sum2*Z+B[i];
    }

    tmp = sum1/sum2;
    *H = creal(tmp);
    *L = cimag(tmp);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Hum_W4(double Nu, double y, double *H, double *L){

/*--------------------------------------------------------------------------------*/
 
    /*######################################################################
      Purpose:
        Voigt function.
      Record of revisions:
        28 Dec, 2018 (Hao Li)
      Input parameters:
        Nu, wavelength.
        y, damping parameter.
      Output parameters:
        erfcx, the result.
      Reference:
        HUMLÍČEK 1982 J. Quant. Spectrosc. Radiat. Transfer.
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    complex double W4, U, T = y-Nu*I;
    double S = fabs(Nu)+y;
    
    if(S >= 15){
      W4 = T*0.56418958355/(0.5+T*T);
    }else if(S >= 5.5){
      U = T*T;
      W4 = T*(1.4104739589+U*0.5641896)/(0.75+U*(3.+U));
    }else if(y >= 1.95*fabs(Nu)-0.176){
      W4 = (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*0.5642236)))) \
        /(16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))));
    }else {
      U = T*T;
      W4 = cexp(U)-T*(36183.30536-U*(3321.990492-U*(1540.786893-U \
          *(219.0312964-U*(35.7668278-U*(1.320521697-U*0.5641900381)))))) \
          /(32066.59372-U*(24322.84021-U*(9022.227659-U*(2186.181081 \
          -U*(364.2190727-U*(61.57036588-U*(1.841438936-U))))))); \
    }
    *H = creal(W4);
    *L = cimag(W4);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/


static int Faddeeva_4digits(double Nu, double y, double *H, double *L){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function with 4-digit accuracy.
      Record of revisions:
        28 Dec, 2018 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H=exp(-Nu*Nu);
      *L=2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y, TMP_y = y*y;
    
    if(TMP>=1.6e4){
      W = I/Z/Par_SqrtPi;
    }else if(TMP>=160.){
      W = I*Z/Par_SqrtPi/(TMP1-0.5);
    }else if(TMP>=107){
      W = (TMP1-1)/(TMP1-1.5)*I/Z/Par_SqrtPi;
    }else if(TMP>=28.5){
      if(TMP_y>=6e-14){
        W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
      }else {
        Hum_W4(Nu, y, H, L);
        return 0;
      }
    }else if((TMP>=3.5)&&(TMP_y<0.026)){
      Hum_W4(Nu, y, H, L);
      return 0;
    }else {
      Hui_p6(Nu, y, H, L);
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_5digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function with 5-digit accuracy.
      Record of revisions:
        9 Jun, 2024 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y, TMP_y = y*y;
    
    if(TMP>=1.5e5){
      W = I/Z/Par_SqrtPi;
    }else if(TMP>=510.){
      W = I*Z/Par_SqrtPi/(TMP1-0.5);
    }else if(TMP>=110){
      W = (TMP1-1)/(TMP1-1.5)*I/Z/Par_SqrtPi;
    }else if(TMP>=109){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP>=39){
      if(TMP_y>=1e-9){
        W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
      }else {
        Hum_W4(Nu, y, H, L);
        return 0;
      }
    }else {
      if(TMP_y>=0.27){
        Hui_p6(Nu, y, H, L);
      }else if(TMP_y>=1e-9){
        Faddeeva916(Nu, y, H, L, Fadd);
      }else{
        Hum_W4(Nu, y, H, L);
      }
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_6digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function with 6-digit accuracy.
      Record of revisions:
        9 Jun, 2024 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I, TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y, TMP_y = y*y;
    
    if(TMP>=1.451e6){
      W = I/(Z*Par_SqrtPi);
    }else if(TMP>=1.6e3){
      W = I*Z/(Par_SqrtPi*(TMP1-0.5));
    }else if(TMP>=180){
      W = (TMP1-1)/(TMP1-1.5)*I/(Z*Par_SqrtPi);
    }else if(TMP>=110){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP_y<1){
      Faddeeva916(Nu, y, H, L, Fadd);
      return 0;
    }else {
      Hui_p6(Nu, y, H, L);
      return 0;
    }
    
    *H=creal(W);
    *L=cimag(W);
    return 0;   
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_7digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function with 7-digit accuracy.
      Record of revisions:
        9 Jun, 2024 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }
    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y;
    
    if(TMP>=1.5e7){
      W = I/(Z*Par_SqrtPi);
    }else if(TMP>=5.01e3){
      W = I*Z/(Par_SqrtPi*(TMP1-0.5));
    }else if(TMP>=380){
      W = (TMP1-1)/(TMP1-1.5)*I/(Z*Par_SqrtPi);
    }else if(TMP>=115){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP>=114){
      W = (TMP1*(TMP1-4.5)+2.0)/(TMP1*(TMP1-5)+3.75)/Par_SqrtPi/Z*I;
    }else {
      Faddeeva916(Nu, y, H, L, Fadd);
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);
    return 0;
    
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_8digits(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function with 8-digit accuracy.
      Record of revisions:
        9 Jun, 2024 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y;
    
    if(TMP>=1.3e8){
      W = I/(Z*Par_SqrtPi);
    }else if(TMP>=1.6e4){
      W = I*Z/(Par_SqrtPi*(TMP1-0.5));
    }else if(TMP>=810){
      W = (TMP1-1)/(TMP1-1.5)*I/(Z*Par_SqrtPi);
    }else if(TMP>=195){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP>=116){
      W = (TMP1*(TMP1-4.5)+2.0)/(TMP1*(TMP1-5)+3.75)/Par_SqrtPi/Z*I;
    }else {
      Faddeeva916(Nu, y, H, L, Fadd);
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);

    return 0;    
}

/*--------------------------------------------------------------------------------*/

extern int Faddeeva(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function.
      Record of revisions:
        9 Jun. 2024 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft.
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

    switch(Fadd->nfigure){
      case 4:
        Faddeeva_4digits(Nu, y, H, L);
        break;
      case 5:
        Faddeeva_5digits(Nu, y, H, L, Fadd);
        break;
      case 6:
        Faddeeva_6digits(Nu, y, H, L, Fadd);
        break;
      case 7:
        Faddeeva_7digits(Nu, y, H, L, Fadd);
        break;
      case 8:
        Faddeeva_8digits(Nu, y, H, L, Fadd);
        break;
      default:
        Faddeeva916(Nu, y, H, L, Fadd);
        break;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Faddeeva_init(STRUCT_FADDEEVA *Fadd){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Precompute some coefficients for Faddeyeva function (algorithm 916).
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        Fadd, a structure with precomputed coefficients .
     ######################################################################*/
    
/*--------------------------------------------------------------------------------*/


    Fadd->Rmin0 = pow(10.,-Fadd->nfigure);
    Fadd->a = Par_Pi/sqrt(-log(0.5*Fadd->Rmin0));
    Fadd->AA = Fadd->a*Fadd->a;
    Fadd->logRmin0 = -log(Fadd->Rmin0);
    Fadd->sqrt_logRmin0 = sqrt(Fadd->logRmin0);
    Fadd->nmax = 50;
    Fadd->Expa2n2 = (double *)malloc(Fadd->nmax*sizeof(double));

    int i;

    Fadd->Expa2n2[0] = 1.;
    for(i=1;i<Fadd->nmax;i++){
      Fadd->Expa2n2[i] = exp(-Fadd->AA*i*i);
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int Faddeeva916(double Nu, double y, double *H, double *L, \
        STRUCT_FADDEEVA *Fadd){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Faddeyeva function (algorithm 916).
      Record of revisions:
        9 Jun 2024 (Hao Li)
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure with precomputed coefficients .
      Output parameters:
        *H, Voigt function.
        *L, associated dispersion profile.
      Reference:
        Zaghloul 2011.
     ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    if(Nu == 0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }else if(y == 0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    const double x = fabs(Nu), ya = fabs(y);

    if(x+ya>1e7){
      if(x<ya){
        double xs = x/ya;
        *H = Par_OneSqrtPi/(xs*xs+1.);
        *L = Par_OneSqrtPi*xs/(xs*xs+1.);
      }else{
        double ys = ya/x;
        *H = Par_OneSqrtPi*ys/(ys*ys+1.);
        *L = Par_OneSqrtPi/(ys*ys+1.);
      }
      return 0;
    }

    const double XX = x*x, YY = y*y, XY = x*y;
    const double XY2 = 2*XY;
    const double SINXY = sin(XY), SIN2XY = sin(XY2), COS2XY = cos(XY2);
    const double EXPXX = exp(-XX), EXP2AX = exp(2*Fadd->a*x);
    const double const1 = EXPXX*Fadd_erfcx(y);

    const int nMax1 = (int)(sqrt(Fadd->logRmin0-x*x)/Fadd->a+0.5);
    const int nMax2 = (int)((Fadd->sqrt_logRmin0-x)/Fadd->a+0.5);
    int nMax = nMax1>nMax2?nMax1:nMax2;
    if(nMax>=50) nMax = 49;

    int Backward = (int)(x/Fadd->a+0.5);
    int Forward = Backward+1;

    const double dx1 = Fadd->a*Forward-x;
    const double dx2 = x-Fadd->a*Backward;
    int nMaxF = (int)((Fadd->sqrt_logRmin0-dx1)/Fadd->a);
    int nMaxB = (Fadd->sqrt_logRmin0-dx2)/Fadd->a;
    if(nMaxF>=50) nMaxF = 49;
    if(nMaxB>=50) nMaxB = 49;

    const int istart = (Backward-nMaxB)>1?(Backward-nMaxB):1;

    int n;
    double Re = const1*COS2XY+2*Fadd->a*x*SINXY*EXPXX*SINXY/XY/Par_Pi;
    double Im = -const1*SIN2XY+2*Fadd->a*x*EXPXX*SIN2XY/XY2/Par_Pi;
    double EXP1, EXP2, EPtmp, delta, sum1=0, sum2=0, sum3=0, sum4=0, sum5=0;

    EXP1 = exp(-dx1*dx1);
    EXP2 = exp(-2*Fadd->a*dx1);
    EPtmp = 1;
    for(n = Forward; n <= Forward+nMaxF; n++){
      delta = EXP1*Fadd->Expa2n2[n-Forward]/(Fadd->AA*n*n+YY);
      sum3 += delta*EPtmp;
      sum5 += delta*EPtmp*Fadd->a*n;
      EPtmp *= EXP2;            
    }

    EXP1 = exp(-dx2*dx2);
    EXP2 = exp(-2*Fadd->a*dx2);
    EPtmp = 1;
    for(n = Backward; n >= istart; n--){
      delta = EXP1*Fadd->Expa2n2[Backward-n]/(Fadd->AA*n*n+YY);
      sum3 += delta*EPtmp;
      sum5 += delta*EPtmp*Fadd->a*n;
      EPtmp *= EXP2;
    }

    Re += Fadd->a*y*sum3/Par_Pi;
    Im += Fadd->a*sum5/Par_Pi;

    if(x < Fadd->sqrt_logRmin0 || x < 10){

      EPtmp = 1;
      for(n = 1; n <= nMax; n++){
        delta = EXPXX*Fadd->Expa2n2[n]/(Fadd->AA*n*n+YY);
        sum1 += delta;
        EPtmp /= EXP2AX;
        sum2 += delta*EPtmp;
        sum4 += delta*EPtmp*Fadd->a*n;
      }

      Re += Fadd->a*y*(-2*COS2XY*sum1+sum2)/Par_Pi;
      Im += Fadd->a*(2*y*SIN2XY*sum1-sum4)/Par_Pi;
    }

    *H = Re;
    if(Nu>0){
      *L = Im;
    }else{
      *L = -Im;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

