
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "ALLOCATION.h"
#include "LOG_ERROR.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        06 Mar. 2026  (Hao Li)
          --- Initial commit:  (moved from SVD.c) 
              1d array is used in SVD instead of the 2d array. 
              The subscripts starts from 0.

        28 Jun. 2024
          --- Initial commit (Hao Li)

     ######################################################################*/

/*--------------------------------------------------------------------------------*/

#ifndef TYPE
#error "TYPE macro must be defined before including svd_template.h"
#endif

#ifndef SVD_FUNC
#error "SVD_FUNC macro must be defined before including svd_template.h"
#endif

#ifndef SVB_FUNC
#error "SVB_FUNC macro must be defined before including svd_template.h"
#endif

#ifndef PYTHAG_FUNC
#error "PYTHAG_FUNC macro must be defined before including svd_template.h"
#endif


#define SVD_MAX_ITER 50

/*--------------------------------------------------------------------------------*/

static inline TYPE PYTHAG_FUNC(TYPE a, TYPE b){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Computes sqrt(a^2+b^2) avoiding destructive underflow or overflow.
      Record of revisions:
        8 Mar. 2026
      Input parameters:        
        a, parameter a.
        b, parameter b.
      Return:
        sqrt(a^2+b^2) .
      Reference:
        Numerical Recipes in C 2nd Edition.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    TYPE absa = a>0 ? a : -a;
    TYPE absb = b>0 ? b : -b;

    TYPE max = absa > absb ? absa : absb;

    if(max==0) return 0;

    if(max < 1e150) return sqrt(a*a+b*b);

    TYPE min = absa>absb ? absb : absa;

    TYPE r = min/max;
    return max*sqrt(1.0+r*r);
}

/*--------------------------------------------------------------------------------*/

int SVD_FUNC(STRUCT_MATRIX *am, TYPE *w, STRUCT_MATRIX *vm){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        singlar value decomposition (SVD).
      Record of revisions:
        8 Mar. 2026
      Input parameters:        
        am, a structure storing the matrix A[0...m-1][0...n-1] and the 
          dimensions m, and n.
      Output parameters:        
        am, A matrix is replaced by U matrix.
        w, the diagnonal matrix.
        vm, a structure storing the matrix v[0...n-1][0...n-1].
      Reference:
        Numerical Recipes in C 2nd Edition.
        Given a matrix A[0...m-1][0...n-1], this routine computes its 
        singular value decomposition, A = U ·W ·V T . The matrix U replaces 
        a on output. The diagonal matrix of singular values W is output 
        as a vector w[0...n-1]. The matrix V (not the transpose V T ) 
        is output as V[0...n-1][0...n-1].
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    int flag, i, its, j, jj, k, l, nm;
    TYPE anorm = 0., c, f, g = 0., h, s, scale = 0., x, y, z, tmp;

    TYPE *a = (TYPE *)am->data;
    TYPE *v = (TYPE *)vm->data;
    int m = am->ni, n = am->nj;
    TYPE rv1[n]; 

    #define a(i,j) a[i*n+j]
    #define v(i,j) v[i*n+j]

    for(i=0; i<n; i++){
      l = i+1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
      if(i<m){
        for(k=i;k<m;k++){
          scale += fabs(a(k,i));
        }
        if(scale){
          for(k=i; k<m; k++){
            a(k,i) /= scale;
            s += a(k,i)*a(k,i);
          }
          f = a(i,i);
          g = f>=0 ? -sqrt(s) : sqrt(s);
          h = f*g-s;
          a(i,i) = f-g;
          for(j= l; j<n; j++){
            for(s=0.0, k=i; k<m; k++){ 
              s += a(k,i)*a(k,j);
            }
            f = s/h;
            for(k=i; k<m; k++){ 
              a(k,j) += f*a(k,i);
            }
          }
          for(k=i; k<m; k++){ 
            a(k,i) *= scale;
          }
        }
      }
      w[i] = scale*g;
      g = s = scale = 0.0;
      if(i<m && i!=n-1){
        for(k=l; k<n; k++){ 
          scale += fabs(a(i,k));
        }
        if(scale){
          for(k = l; k<n; k++){
            a(i,k) /= scale;
            s += a(i,k)*a(i,k);
          }
          f = a(i,l);
          g = f>=0 ? -sqrt(s) : sqrt(s);
          h = f*g-s;
          a(i,l) = f-g;
          for(k=l; k<n; k++){ 
            rv1[k] = a(i,k)/h;
          }
          for(j=l; j<m; j++){
            for(s = 0.0, k=l; k<n; k++){ 
              s += a(j,k)*a(i,k);
            }
            for(k=l; k<n; k++){ 
              a(j,k) += s*rv1[k];
            }
          }
          for(k=l; k<n; k++){ 
            a(i,k) *= scale;
          }
        }
      }
      tmp = fabs(w[i])+fabs(rv1[i]);
      if(anorm<tmp){
        anorm = tmp;
      }
    }

    for(i=n-1; i>=0; i--){
      if(i<n-1){
        if(g!= 0.0){
          for(j=l; j<n; j++){
            v(j,i) = (a(i,j)/a(i,l))/g;
          }
          for(j=l; j<n; j++){
            for(s=0.0, k=l; k<n; k++){
              s += a(i,k)*v(k,j);
            }
            for(k=l; k<n; k++){
              v(k,j) += s*v(k,i);
            }
          }
        }
        for(j=l; j<n; j++){
          v(i,j) = v(j,i) = 0.0;
        }
      }
      v(i,i) = 1.0;
      g = rv1[i];
      l = i;
    }

    int itmp = m<=n ? m : n;
    for(i=itmp-1; i>=0; i--){
      l = i+1;
      g = w[i];
      for(j=l; j<n; j++){ 
        a(i,j) = 0.0;
      }
      if(g!=0.0){
        g = 1.0/g;
        for(j=l; j<n; j++){
          for(s=0.0, k=l; k<m; k++){ 
            s += a(k,i)*a(k,j);
          }
          f = (s/a(i,i))*g;
          for(k=i; k<m; k++){ 
            a(k,j) += f*a(k,i);
          }
        }
        for(j=i; j<m; j++){ 
          a(j,i) *= g;
        }
      }else{
        for(j = i; j<m; j++){ 
          a(j,i) = 0.0;
        }
      }
      a(i,i) += 1.0;
    }

    for(k=n-1; k>=0; k--){
      for(its=0; its<SVD_MAX_ITER; its++){
        flag = 1;
        for(l=k; l >= 0; l--){
          nm = l-1;
          if((fabs(rv1[l])+anorm)==anorm){
            flag = 0;
            break;
          }
          if(nm>=0 && (fabs(w[nm])+anorm)==anorm) break;
        }
        if(flag){
          c = 0.0;
          s = 1.0;
          for(i=l; i<=k; i++){
            f = s*rv1[i];
            rv1[i] = c*rv1[i];
            if((fabs(f)+anorm)==anorm) break;
            g = w[i];
            h = PYTHAG_FUNC(f, g);
            w[i] = h;
            h = 1.0/h;
            c = g*h;
            s = -f*h;
            for(j=0; j<m; j++){
              y = a(j,nm);
              z = a(j,i);
              a(j,nm) = y*c+z*s;
              a(j,i) = z*c-y*s;
            }
          }
        }
        z = w[k];
        if(l==k){
          if(z<0.0){
            w[k] = -z;
            for(j=0; j<n; j++){ 
              v(j,k) = -v(j,k);
            }
          }
          break;
        }
        if(its==SVD_MAX_ITER-1){
          LOG_ERROR(ERR_LVL_ERROR, "svdcmp", \
              "svdcmp reach maxium iterations\n");
        }
        x = w[l];
        nm = k-1;
        y = w[nm];
        g = rv1[nm];
        h = rv1[k];
        f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g = PYTHAG_FUNC(f, 1.0);
        double dtmp = f>=0.0 ? g : -g;
        f = ((x-z)*(x+z)+h*(y/(f+dtmp)-h))/x;
        c = s = 1.0;
        for(j=l; j<=nm; j++){
          i = j+1;
          g = rv1[i];
          y = w[i];
          h = s*g;
          g *= c;
          z = PYTHAG_FUNC(f, h);
          rv1[j] = z;
          c = f/z;
          s = h/z;

          f = x*c+g*s;
          g = g*c-x*s;
          h = y*s;
          y *= c;
          for(jj=0;jj<n;jj++){
            x = v(jj,j);
            z = v(jj,i);
            v(jj,j) = x*c+z*s;
            v(jj,i) = z*c-x*s;
          }
          z=PYTHAG_FUNC(f,h);
          w[j]=z;
          if(z){
            z = 1.0/z;
            c = f*z;
            s = h*z;
          }

          f = c*g+s*y;
          x = c*y-s*g;

          for(jj=0; jj<m; jj++){
            y = a(jj,j);
            z = a(jj,i);
            a(jj,j) = y*c+z*s;
            a(jj,i) = z*c-y*s;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

int SVB_FUNC(STRUCT_MATRIX *am, STRUCT_MATRIX *vm, const TYPE *w, \
    const TYPE *b, TYPE *x){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Solveing A·X = B for a vector X.
      Record of revisions:
        8 Mar. 2026
      Input parameters:        
        am, a structure storing the matrix U.
        vm, a structure storing the matrix v[0...n-1][0...n-1].
        w, the diagnonal matrix.
        b, the B vector.
      Output parameters:        
        x, the solutions.
      Reference:
        Numerical Recipes in C 2nd Edition.
        Solves A·X = B for a vector X, where A is specified by the arrays 
        u[0..m-1][0..n-1], w[0..n-1], v[0..n-1][0..n-1] as returned by 
        svdcmp. m and n are the dimensions of a, and will be equal for 
        square matrices. b[1..m] is the input right-hand side. x[1..n] is 
        the output solution vector. No input quantities are destroyed, so 
        the routine may be called sequentially with different b’s.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    #define TINY 1e-50

    int m = am->ni, n = am->nj;
    TYPE *u = (TYPE *)am->data;
    TYPE *v = (TYPE *)vm->data;          
    TYPE tmp[n]; 
    TYPE s;

    #define u(i,j) u[i*n+j]
    #define v(i,j) v[i*n+j]

    for(int jj=0; jj<n; jj++){
      s = 0.0;
      if(w[jj]>TINY){
        for(int ii=0; ii<m; ii++){
          s += u(ii,jj)*b[ii]; 
        }
        s /= w[jj];
      }
      tmp[jj] = s;
    }

    for(int ii=0; ii<n; ii++){
      s = 0.0;
      for(int jj=0; jj<n; jj++){
        s += v(ii,jj)*tmp[jj]; 
      }
      x[ii] = s;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/
