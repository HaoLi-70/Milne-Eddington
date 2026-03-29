
#include "FADDEEVA.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        07 Mar. 2026  (Hao Li)
          --- Updates: 
              Use doubles for real/imag parts to avoid complex double. 
              Change some functions to "static inline".
              Redesign Faddeeva function.

        9 Jun. 2024
          --- update: add a structure STRUCT_FADDEEVA with some precomputed 
              coefficients.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

// square root of Pi
static const double L_SqrtPi = 1.77245385090551588191;

/*--------------------------------------------------------------------------------*/

static inline double DAWSON(double x){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Dawson’s Integral.
      Record of revisions:
        6  Mar. 2026.
      Input parameters:
        x, the real number.
      Output parameters:
        ans, the integral.
      Reference:
        Numerical recipes in C, 2nd Edition.
     ######################################################################*/
/*--------------------------------------------------------------------------------*/

    const double A1 = 2.0/3.0, A2 = 0.4, A3 = 2.0/7.0, H = 0.4;
    const int NMAX = 6;
    double c[NMAX+1];
    double xx = fabs(x), ans;
    
    if(xx<0.2){  
      double x2 = x*x; 
      ans = x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));    
    }else {
      int i;
      double TMP;
      for(i=1;i<=NMAX;i++){ 
        TMP = (2.0*i-1.0)*H;    
        c[i] = exp(-TMP*TMP);    
      }
        
      int n0 = 2*(int)round(xx/(2*H));
      
      double xp = xx-n0*H;
      double e1 = exp(2.0*xp*H);
      double e2 = e1*e1;  
      double d1 = n0+1;
      double d2 = d1-2.0;  
      double sum = 0.0;
        
      for(i=1; i<=NMAX; i++,d1+=2.0,d2-=2.0,e1*=e2){
        sum += c[i]*(e1/d1+1.0/(d2*e1));
      }
      TMP = exp(-xp*xp);
      ans = 0.5641895835*copysign(TMP, x)*sum;
    }
    
    return ans;
}

/*--------------------------------------------------------------------------------*/

static inline double ERFCX(double a){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Scaled complementary error function (erfcx(x) = e^{x^2} * erfc(x))
      Record of revisions::
        6  Mar. 2026.
      Input parameters:
        a, the real number.
      Output parameters:
        erfcx, the result.
      Reference:
        Aghloul 2007 MNRAS.
     ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(a>26.6){
      double asqr = a*a;
      return ((((((162.421875/asqr-29.53125)/asqr+6.5625)/asqr-1.875) \
          /asqr+0.75)/asqr-0.5)/asqr+1)/L_SqrtPi/a;

    }else{
      return exp(a*a)*erfc(a);
    }
}

/*--------------------------------------------------------------------------------*/

static inline int Hui_p6(double Nu, double y, double *H, double *L){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Voigt function（Hui 1978）.
      Record of revisions:
        6  Mar. 2026.
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

    static const double A[7] = {122.607931777104326,                   \
        214.382388694706425, 181.928533092181549, 93.155580458138441,  \
        30.180142196210589, 5.912626209773153, 0.564189583562615};
    static const double B[7] = {122.60793177387535,                    \
        352.730625110963558, 457.334478783897737, 348.703917719495792, \
        170.354001821091472, 53.992906912940207, 10.479857114260399};

    double Z_re = y, Z_im = -Nu;
    double sum1_re = A[6], sum1_im = 0.0;
    double sum2_re = B[6]+Z_re, sum2_im = Z_im;
    double tmp_re, tmp_im;

    for(int ii=5; ii>=0; ii--){
      // sum1 = sum1*Z+A[i]
      tmp_re = sum1_re*Z_re-sum1_im*Z_im+A[ii];
      tmp_im = sum1_re*Z_im+sum1_im*Z_re;
      sum1_re = tmp_re; 
      sum1_im = tmp_im;

      // sum2 = sum2*Z+B[i]
      tmp_re = sum2_re*Z_re-sum2_im*Z_im+B[ii];
      tmp_im = sum2_re*Z_im+sum2_im*Z_re;
      sum2_re = tmp_re; 
      sum2_im = tmp_im;
    }

    // sum1/sum2;
    double denom = sum2_re*sum2_re+sum2_im*sum2_im;
    *H = (sum1_re*sum2_re+sum1_im*sum2_im)/denom;
    *L = (sum1_im*sum2_re-sum1_re*sum2_im)/denom;

    return 0;
}

/*--------------------------------------------------------------------------------*/

static inline int Hum_W4(double Nu, double y, double *H, double *L){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Voigt function (HUMLÍČEK 1982).
      Record of revisions:
        6 Mar. 2026.
      Input parameters:
        Nu, wavelength.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        HUMLÍČEK 1982 J. Quant. Spectrosc. Radiat. Transfer.
     ######################################################################*/
/*--------------------------------------------------------------------------------*/

    const double absNu = fabs(Nu);
    const double S = absNu +y;
    const double T_re = y, T_im = -Nu;

    if(S >= 15){
      const double U_re = T_re*T_re-T_im*T_im;
      const double U_im = 2.0*T_re*T_im;
      double D_re = 0.5+U_re, D_im = U_im;
      double denom = D_re*D_re+ D_im*D_im;

      *H = 0.56418958355*(T_re*D_re+T_im*D_im)/denom;
      *L = 0.56418958355*(T_im*D_re-T_re*D_im)/denom;

    }else if(S >= 5.5){
      const double U_re = T_re*T_re - T_im*T_im;
      const double U_im = 2.0*T_re*T_im;
      double N_re = 1.4104739589 + 0.5641896*U_re;
      double N_im = 0.5641896*U_im;
      double D_re = 0.75+3.0*U_re+(U_re*U_re-U_im*U_im);
      double D_im = 3.0*U_im+2.0*U_re*U_im;
      double denom = D_re*D_re+D_im*D_im;
      double tmp_re = (N_re*D_re+N_im*D_im)/denom;
      double tmp_im = (N_im*D_re-N_re*D_im)/denom;

      *H = T_re*tmp_re-T_im*tmp_im;
      *L = T_re*tmp_im+T_im*tmp_re;

    }else if(y >= 1.95*absNu-0.176){
      const double c[5] = {0.5642236, 3.778987, 11.96482, 20.20933,     \
          16.4955};  
      const double d[6] = {1.0, 6.699398, 21.69274, 39.27121, 38.82363, \
          16.4955};

      double N_re = c[0], N_im = 0;
      double D_re = d[0], D_im = 0;
      double tmp_re, tmp_im;

      for(int ii=1; ii<5; ii++){
        tmp_re = N_re*T_re-N_im*T_im+c[ii];
        tmp_im = N_re*T_im+N_im*T_re;
        N_re = tmp_re;
        N_im = tmp_im;
      }

      for(int ii=1; ii<6; ii++){
        tmp_re = D_re*T_re-D_im*T_im + d[ii];
        tmp_im = D_re*T_im+D_im*T_re;
        D_re = tmp_re;
        D_im = tmp_im;
      }
      double denom = D_re*D_re+D_im*D_im;

      *H = (N_re*D_re+N_im*D_im)/denom;
      *L = (N_im*D_re-N_re*D_im)/denom;

    }else{
      const double c[7] = {0.5641900381, 1.320521697, 35.7668278,    \
          219.0312964, 1540.786893, 3321.990492, 36183.30536};
      const double d[8] = {1, 1.841438936, 61.57036588, 364.2190727, \
          2186.181081, 9022.227659, 24322.84021, 32066.59372};

      const double U_re = T_re*T_re-T_im*T_im;
      const double U_im = 2.0*T_re*T_im;
      double expU = exp(U_re);
      double exp_re = expU*cos(U_im);
      double exp_im = expU*sin(U_im);
      double N_re = c[0], N_im = 0;
      double D_re = d[0], D_im = 0;
      double tmp_re, tmp_im;

      for(int ii=1; ii<7; ii++){
        tmp_re = -(N_re*U_re-N_im*U_im)+c[ii];
        tmp_im = -(N_re*U_im+N_im*U_re);
        N_re = tmp_re;
        N_im = tmp_im;
      }

      tmp_re = N_re*T_re-N_im*T_im;
      tmp_im = N_re*T_im+N_im*T_re;
      N_re = tmp_re;
      N_im = tmp_im;

      for(int ii=1; ii<8; ii++){
        tmp_re = -(D_re*U_re-D_im*U_im)+d[ii];
        tmp_im = -(D_re*U_im+D_im*U_re);
        D_re = tmp_re;
        D_im = tmp_im;
      }
      double denom = D_re*D_re+D_im*D_im;

      *H = exp_re-(N_re*D_re+N_im*D_im)/denom;
      *L = exp_im-(N_im*D_re-N_re*D_im)/denom;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/


int Faddeeva916(double Nu, double y, double *H, double *L, \
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
      *H = ERFCX(y);
      *L = 0;
      return 0;
    }else if(y == 0){
      *H = exp(-Nu*Nu);
      *L = 2.0/L_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    const double x = fabs(Nu), ya = fabs(y);

    if(x+ya>1e7){
      if(x<ya){
        double xs = x/ya;
        *H = 1./L_SqrtPi/(xs*xs+1.);
        *L = xs/L_SqrtPi/(xs*xs+1.);
      }else{
        double ys = ya/x;
        *H = ys/L_SqrtPi/(ys*ys+1.);
        *L = 1./L_SqrtPi/(ys*ys+1.);
      }
      return 0;
    }

    const double XX = x*x, YY = y*y, XY = x*y;
    const double XY2 = 2*XY;
    const double SINXY = sin(XY), SIN2XY = sin(XY2), COS2XY = cos(XY2);
    const double EXPXX = exp(-XX), EXP2AX = exp(2*Fadd->a*x);
    const double const1 = EXPXX*ERFCX(y);

    const int nMax1 = (int)(sqrt(Fadd->logRmin0-x*x)/Fadd->a+0.5);
    const int nMax2 = (int)((Fadd->sqrt_logRmin0-x)/Fadd->a+0.5);
    int nMax = nMax1>nMax2?nMax1:nMax2;
    if(nMax>=Fadd->nmax) nMax = Fadd->nmax-1;

    int nBackward = (int)(x/Fadd->a+0.5);
    int nForward = nBackward+1;

    const double dxForward = Fadd->a*nForward-x;
    const double dxBackward = x-Fadd->a*nBackward;

    int nMaxF = (int)((Fadd->sqrt_logRmin0-dxForward)/Fadd->a);
    int nMaxB = (int)((Fadd->sqrt_logRmin0-dxBackward)/Fadd->a);
    if(nMaxF>=Fadd->nmax) nMaxF = Fadd->nmax-1;
    if(nMaxB>=Fadd->nmax) nMaxB = Fadd->nmax-1;

    const int istart = (nBackward-nMaxB)>1?(nBackward-nMaxB):1;
    
    double Re = const1*COS2XY+2*Fadd->a*x*SINXY*EXPXX*SINXY/XY/M_PI;
    double Im = -const1*SIN2XY+2*Fadd->a*x*EXPXX*SIN2XY/XY2/M_PI;
    double EXP1, EXP2, EPtmp, delta; 
    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

    EXP1 = exp(-dxForward*dxForward);
    EXP2 = exp(-2*Fadd->a*dxForward);
    EPtmp = 1.;
    for(int n = nForward; n <= nForward+nMaxF; n++){
      delta = EXP1*Fadd->Expa2n2[n-nForward]/(Fadd->AA*n*n+YY);
      sum3 += delta*EPtmp;
      sum5 += delta*EPtmp*Fadd->a*n;
      EPtmp *= EXP2;            
    }

    EXP1 = exp(-dxBackward*dxBackward);
    EXP2 = exp(-2*Fadd->a*dxBackward);
    EPtmp = 1.;
    for(int n = nBackward; n >= istart; n--){
      delta = EXP1*Fadd->Expa2n2[nBackward-n]/(Fadd->AA*n*n+YY);
      sum3 += delta*EPtmp;
      sum5 += delta*EPtmp*Fadd->a*n;
      EPtmp *= EXP2;
    }

    Re += Fadd->a*y*sum3/M_PI;
    Im += Fadd->a*sum5/M_PI;

    if(x < Fadd->sqrt_logRmin0 || x < 10){
      EPtmp = 1.;
      for(int n = 1; n <= nMax; n++){
        delta = EXPXX*Fadd->Expa2n2[n]/(Fadd->AA*n*n+YY);
        sum1 += delta;
        EPtmp /= EXP2AX;
        sum2 += delta*EPtmp;
        sum4 += delta*EPtmp*Fadd->a*n;
      }

      Re += Fadd->a*y*(-2*COS2XY*sum1+sum2)/M_PI;
      Im += Fadd->a*(2*y*SIN2XY*sum1-sum4)/M_PI;
    }

    *H = Re;
    *L = Nu>0 ? Im : -Im;

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Faddeeva_init(STRUCT_FADDEEVA *Fadd){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Precompute some coefficients for Faddeyeva function (algorithm 916).
      Record of revisions:
        6  Mar. 2026.
      Input parameters:
        Fadd, a structure with precomputed coefficients.
      Output parameters:
        Fadd, a structure with precomputed coefficients.
     ######################################################################*/
/*--------------------------------------------------------------------------------*/

    Fadd->Rmin0 = pow(10.,-Fadd->nfigures);
    Fadd->a = M_PI/sqrt(-log(0.5*Fadd->Rmin0));
    Fadd->AA = Fadd->a*Fadd->a;
    Fadd->logRmin0 = -log(Fadd->Rmin0);
    Fadd->sqrt_logRmin0 = sqrt(Fadd->logRmin0);
    Fadd->nmax = 50;
    Fadd->Expa2n2 = (double *)malloc(Fadd->nmax*sizeof(double));
    Fadd->indxfig = Fadd->nfigures>4?Fadd->nfigures-4:0;

    Fadd->Expa2n2[0] = 1.;
    for(int ii=1; ii<Fadd->nmax; ii++){
      Fadd->Expa2n2[ii] = exp(-Fadd->AA*ii*ii);
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Faddeeva(double Nu, double y, double *H, double *L, \
    STRUCT_FADDEEVA *Fadd){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Faddeyeva function with required accuracy (nfigures). 
      Record of revisions:
        7 Mar. 2026.
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        Fadd, a structure storing precomputed coefficients and nfigures.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft.
     ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(Nu == 0){
      *H = ERFCX(y);
      *L = 0;
      return 0;
    }else if(y == 0){
      *H = exp(-Nu*Nu);
      *L = 2.0/L_SqrtPi*DAWSON(Nu);
      return 0;
    }

    static const double THRESHOLD[5][5] = {         \
        {1.6e4, 160.0, 107.0, 28.5, 3.5 },          \
        {1.5e5, 510.0, 110.0, 39., -0.1 },          \
        {1.451e6, 1.6e3, 180.0, 110.0, -0.1 },      \
        {1.5e7, 5.01e3, 380.0, 115.0, 114.0 },      \
        {1.3e8, 1.6e4, 810.0, 195.0, 116.0 }        \
    };

    const double Z_re = Nu; const double Z_im = y;
    const double Z2_re = Z_re*Z_re-Z_im*Z_im; 
    const double Z2_im = 2.0*Z_re*Z_im;
    const double TMP = Z_re*Z_re+Z_im*Z_im;
    double TMP_y = y*y;

    if(TMP >= THRESHOLD[Fadd->indxfig][0]){
      //I/Z/L_SqrtPi;
      double denom = L_SqrtPi*(Z_re*Z_re+Z_im*Z_im);
      *H = Z_im/denom;
      *L = Z_re/denom;

    }else if(TMP >= THRESHOLD[Fadd->indxfig][1]){
      //I*Z/L_SqrtPi/(Z2-0.5);
      double Z2n_re = Z2_re-0.5;
      double denom = L_SqrtPi*(Z2n_re*Z2n_re+Z2_im*Z2_im);
      *H = (-Z_im*Z2n_re+Z_re*Z2_im)/denom;
      *L = (Z_re*Z2n_re+Z_im*Z2_im)/denom;

    }else if(TMP >= THRESHOLD[Fadd->indxfig][2]){
      //(Z2-1)/(Z2-1.5)*I/Z/L_SqrtPi;
      double Z2n_re = Z2_re-1.5;
      double Z2Z_re = Z2n_re*Z_re-Z2_im*Z_im;
      double Z2Z_im = Z2n_re*Z_im+Z2_im*Z_re;                                 
      double denom = L_SqrtPi*(Z2Z_re*Z2Z_re+Z2Z_im*Z2Z_im);
      Z2n_re = Z2_re-1.0;
      *H = (-Z2_im*Z2Z_re+Z2n_re*Z2Z_im)/denom;
      *L = (Z2n_re*Z2Z_re+Z2_im*Z2Z_im)/denom;

    }else if(TMP >= THRESHOLD[Fadd->indxfig][3]){

      if(Fadd->nfigures<=4 && TMP_y<6e-14){
        Hum_W4(Nu, y, H, L);
        return 0;
      }else if(Fadd->nfigures==5 && TMP<39 && TMP_y<1e-9){
        Hum_W4(Nu, y, H, L);
        return 0;
      }

      //(Z2-2.5)/(Z2*(Z2-3)+0.75)/L_SqrtPi*Z*I;
      double Z2n_re = Z2_re-2.5;
      // A = (Z2-2.5)*Z*I;
      double A_re = -(Z2n_re*Z_im+Z2_im*Z_re);
      double A_im = Z2n_re*Z_re-Z2_im*Z_im;

      Z2n_re = Z2_re-3.;
      // B = Z2*(Z2-3)+0.75
      double B_re = Z2n_re*Z2_re-Z2_im*Z2_im+0.75;
      double B_im = (Z2n_re+Z2_re)*Z2_im;

      double denom = L_SqrtPi*(B_re*B_re+B_im*B_im);

      *H = (A_re*B_re+A_im*B_im)/denom;
      *L = (-A_re*B_im+A_im*B_re)/denom;

    }else if(TMP >= THRESHOLD[Fadd->indxfig][4]){
      if(Fadd->nfigures<=4){
        if(TMP_y<0.026){
          Hum_W4(Nu, y, H, L);
        }else{
          Hui_p6(Nu, y, H, L);
        }
        return 0;
      }else if(Fadd->nfigures==5){
        if(TMP_y>=0.27){
          Hui_p6(Nu, y, H, L);
        }else if(TMP_y>=1e-9){
          Faddeeva916(Nu, y, H, L, Fadd);
        }else{
          Hum_W4(Nu, y, H, L);
        }
        return 0;
      }else if(Fadd->nfigures==6){
        if(TMP_y>=1){
          Hui_p6(Nu, y, H, L);
        }else{
          Faddeeva916(Nu, y, H, L, Fadd);
        }
        return 0;
      }

      //(Z2*(Z2-4.5)+2.0)/(Z2*(Z2-5)+3.75)/L_SqrtPi/Z*I;
      double Z2n_re = Z2_re-4.5;
      // A = (Z2*(Z2-4.5)+2.0)*I
      double A_re = -(Z2n_re+Z2_re)*Z2_im;
      double A_im = Z2n_re*Z2_re-Z2_im*Z2_im+2.;

      Z2n_re = Z2_re-5.;
      // B = Z2*(Z2-5)+3.75
      double B_re = Z2n_re*Z2_re-Z2_im*Z2_im+3.75;
      double B_im = (Z2n_re+Z2_re)*Z2_im;

      // C = B*Z
      double C_re = B_re*Z_re-B_im*Z_im;
      double C_im = B_re*Z_im+B_im*Z_re;

      double denom = L_SqrtPi*(C_re*C_re+C_im*C_im);

      *H = (A_re*C_re+A_im*C_im)/denom;
      *L = (-A_re*C_im+A_im*C_re)/denom;

    }else{
      if(Fadd->nfigures<=4){
        Hui_p6(Nu, y, H, L);
        return 0;
      }else{
        Faddeeva916(Nu, y, H, L, Fadd);
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/
