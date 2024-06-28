
#include "SVD.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        28 Jun. 2024
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

static double pythag(double a, double b){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Computes sqrt(a^2 + b^2) without destructive underflow or overflow.
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:        
        a, parameter a.
        b, parameter b.
      Return:
        sqrt(a^2 + b^2) .
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double absa = fabs(a);
    double absb = fabs(b);

    if (absa > absb){

      return absa*sqrt(1.0+M_SQR(absb/absa));

    }else{

      return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+M_SQR(absa/absb)));

    }
}

/*--------------------------------------------------------------------------------*/

extern int svdcmp(double **a, int m, int n, double w[], double **v){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        singlar value decomposition (SVD).
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:        
        a, matrix a[m][n].
        m, n, the sizes.
      Output parameters:        
        a, the U matrix.
        w, the diagnonal matrix.
        v, the matrix V.
      Reference:
        numerical recipes in C 2ed.
        Given a matrix a[1..m][1..n], this routine computes its singular 
        value decomposition, A = U ·W ·V T . The matrix U replaces a 
        on output. The diagonal matrix of singular values W is output 
        as a vector w[1..n]. The matrix V (not the transpose V T ) 
        is output as v[1..n][1..n].
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int flag, i, its, j, jj, k, l, nm, itmp;
    double anorm, c, f, g, h, s, scale, x, y, z, *rv1, dtmp;
    
    rv1 = (double *)VECTOR(1, n, enum_dbl, true);

    g = scale = anorm = 0.0;
    for (i=1;i<=n;i++) {
      l = i+1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
      if (i <= m) {
        for (k=i;k<=m;k++) scale += fabs(a[k][i]);
        if (scale) {
          for (k=i;k<=m;k++){
            a[k][i] /= scale;
            s += a[k][i]*a[k][i];
          }
          f = a[i][i];
          //g = -M_SIGN(sqrt(s),f);
          g = f >= 0? -sqrt(s):sqrt(s);
          h=f*g-s;
          a[i][i] = f-g;
          for (j=l;j<=n;j++) {
            for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
            f=s/h;
            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
          }
          for (k=i;k<=m;k++) a[k][i] *= scale;
        }
      }
      w[i] = scale*g;
      g = s = scale = 0.0;
      if (i <= m && i != n) {
        for (k=l;k<=n;k++) scale += fabs(a[i][k]);
        if (scale) {
          for (k=l;k<=n;k++) {
            a[i][k] /= scale;
            s += a[i][k]*a[i][k];
          }
          f = a[i][l];
          //g = -M_SIGN(sqrt(s),f);
          g = f >= 0? -sqrt(s):sqrt(s);

          h = f*g-s;
          a[i][l] = f-g;
          for (k=l;k<=n;k++) rv1[k] = a[i][k]/h;
          for (j=l;j<=m;j++) {
            for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
          }
          for (k=l;k<=n;k++) a[i][k] *= scale;
        }
      }
      //anorm = M_MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
      anorm = anorm>=(fabs(w[i])+fabs(rv1[i]))? \
          anorm : (fabs(w[i])+fabs(rv1[i]));

    }
    for (i=n;i>=1;i--) {
      if (i < n) {
        if (g) {
          for (j=l;j<=n;j++) v[j][i] = (a[i][j]/a[i][l])/g; 
          for (j=l;j<=n;j++) {
              for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
              for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
          }
        }
        for (j=l;j<=n;j++) v[i][j] = v[j][i] = 0.0;
      }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }
    itmp = m<=n? m:n;
    for (i=itmp;i>=1;i--) {
      l = i+1;
      g = w[i];
      for (j=l;j<=n;j++) a[i][j] = 0.0;
      if (g) {
        g = 1.0/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
          f = (s/a[i][i])*g;
          for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
        }
        for (j=i;j<=m;j++) a[j][i] *= g;
      }else{
        for (j=i;j<=m;j++) a[j][i]=0.0;
      }
      ++a[i][i];
    }for (k=n;k>=1;k--) {
      for (its=1;its<=30;its++){
        flag = 1;
        for (l=k;l>=1;l--) {
          nm = l-1;
          if ((double)(fabs(rv1[l])+anorm) == anorm) {
            flag=0;
            break;
          }
          if ((double)(fabs(w[nm])+anorm) == anorm) break;
        }
        if (flag) {
          c = 0.0;
          s = 1.0;
          for (i=l;i<=k;i++) {
            f = s*rv1[i];
            rv1[i] = c*rv1[i];
            if ((double)(fabs(f)+anorm) == anorm) break;
            g=w[i];
            h = pythag(f,g);
            w[i] = h;
            h = 1.0/h;
            c = g*h;
            s = -f*h;
            for (j=1;j<=m;j++) {
              y = a[j][nm];
              z = a[j][i];
              a[j][nm] = y*c+z*s;
              a[j][i] = z*c-y*s;
            }
          }
        }
        z = w[k];
        if (l == k) {
          if (z < 0.0) {
            w[k] = -z;
            for (j=1;j<=n;j++) v[j][k] = -v[j][k];
          }
          break;
        }
        if (its == 30){ 
          // nrerror("no convergence in 30 svdcmp iterations");
          fprintf(stderr,"svdcmp fail 1 \n");
        }
        x = w[l];
        nm = k-1;
        y = w[nm];
        g = rv1[nm];
        h = rv1[k]; f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g = pythag(f,1.0);
        //tmp = M_SIGN(g,f);
        dtmp = f >= 0?g:-g;
        f = ((x-z)*(x+z)+h*((y/(f+dtmp))-h))/x;
        c = s = 1.0;
        for (j=l;j<=nm;j++) {
          i = j+1;
          g = rv1[i];
          y = w[i];
          h = s*g;
          g = c*g;
          z = pythag(f,h);
          rv1[j] = z;
          c = f/z;
          s = h/z;
          f = x*c+g*s;
          g = g*c-x*s;
          h = y*s;
          y *= c;
          for (jj=1;jj<=n;jj++) {
            x = v[jj][j];
            z = v[jj][i];
            v[jj][j] = x*c+z*s;
            v[jj][i] = z*c-x*s;
          }
          z=pythag(f,h);
          w[j]=z;
          if (z) {
            z = 1.0/z;
            c = f*z;
            s = h*z;
          }
          f = c*g+s*y;
          x = c*y-s*g;
          for (jj=1;jj<=m;jj++) {
            y = a[jj][j];
            z = a[jj][i];
            a[jj][j] = y*c+z*s;
            a[jj][i] = z*c-y*s;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }
    FREE_VECTOR(rv1, 1, enum_dbl);
      
    return  0;
}

/*--------------------------------------------------------------------------------*/

extern int svbksb(double **u, double w[], double **v, int m, int n, \
        double b[], double x[]){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Solveing A·X = B for a vector X.
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:        
        u, the U matrix.
        w, the diagnonal matrix.
        v, the matrix V.        
        m, n, the sizes.
        b, the B vector.
      Output parameters:        
        x, the solutions.
      Reference:
        numerical recipes in C 2ed.
        Solves A·X = B for a vector X, where A is specified by the arrays 
        u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m 
        and n are the dimensions of a, and will be equal for square 
        matrices. b[1..m] is the input right-hand side. x[1..n] is the 
        output solution vector. No input quantities are destroyed, so 
        the routine may be called sequentially with different b’s.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int jj,j,i;
    double s,  ZERO_TINY = 1e-50;
    double *tmp = (double *)VECTOR(1, n, enum_dbl, true);
  
    for (j=1;j<=n;j++){
      s = 0.0;
      if (w[j]>ZERO_TINY){
        for (i=1;i<=m;i++) s += u[i][j]*b[i];
        s /= w[j];
      }
      tmp[j] = s;
    }
    for (j=1;j<=n;j++) {
      s = 0.0;
      for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
      x[j] = s;
    }
    FREE_VECTOR(tmp, 1, enum_dbl);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int SVD_solve(double **A, double *B, double *X, int N, \
        double Threshold){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Solveing A·X = B for a vector X.
      Record of revisions:
        28 Jun. 2024 (Hao Li)
      Input parameters:        
        A, the A matrix.
        B, the B vector.      
        N, the size.
        Threshold, the threshold for the the diagnonal matrix w.
      Output parameters:        
        X, the solutions.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    double *W = (double *)VECTOR(1, N, enum_dbl, false);
    double **V = (double **)MATRIX(1, N, 1, N, enum_dbl, false);
    
    svdcmp(A, N, N, W, V);
    
    double WMAX = W[1];
    int i;
    
    for (i=2; i<=N; i++){ 
      if(WMAX<W[i]) WMAX = W[i];
    }

    WMAX *= Threshold;

    for (i=1; i<=N; i++){ 
      if(W[i]<WMAX){
        fprintf(stderr,"%d %e \n",i,W[i]);
        W[i] = 0.;
      }
    }

    svbksb(A, W, V, N, N, B, X);

    FREE_VECTOR(W, 1, enum_dbl);
    FREE_MATRIX((void *)V, 1, 1, enum_dbl);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/
