
#include "RANDOM_NUMBER.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        30 Oct. 2022
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern double Ran1(long *idum){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Generate a random number between 0 and 1 (1 is not included) .
      Record of revisions:
        30 Oct. 2022 (Hao Li)
      Input parameters:
        *idum, the seed for the random number.
      Return:
        The randum number.
      Reference:
        numerical recipes in C 2ed. 
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if (*idum <= 0 || !iy) {
      if (-(*idum) < 1)
        *idum = 1;
      else
        *idum = -(*idum);
  
      for(j=NTAB+7;j>=0;j--) {
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0)
          *idum += IM;
        if (j < NTAB)
          iv[j] = *idum;
      }
      iy=iv[0];
    }
    
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    
    if (*idum < 0)
      *idum += IM;
    
    j=(int)iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    
    if ((temp=AM*iy) > RNMX)
      return RNMX;
    else
      return temp;
}

/*--------------------------------------------------------------------------------*/

extern double GASDEV(long *idum){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Generate a Gaussian distributed random number with width of 1.
      Record of revisions:
        30 Oct. 2022 (Hao Li)
      Input parameters:
        *idum, the seed for the random number.
      Return:
        The randum number.
      Reference:
        numerical recipes in C 2ed. 
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;
    
    if (*idum < 0)
      iset=0;
    
    if (iset == 0) {
      do {
        v1=2.0*Ran1(idum)-1.0;    
        v2=2.0*Ran1(idum)-1.0;
        rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);

      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
    } else {
      iset=0;
      return gset;
    }
}

/*--------------------------------------------------------------------------------*/

extern int Random_Seed(long *idum){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Generate the seed for the random number.
      Record of revisions:
        30 Oct. 2022 (Hao Li)
      Output parameters:
        *idum, the seed for the random number.
      Return:
        .
      Reference:
        numerical recipes in C 2ed. 
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    srand((unsigned int)time(NULL));

    srand(rand());    
    
    *idum = -rand();
    
    return 0;
}

/*--------------------------------------------------------------------------------*/