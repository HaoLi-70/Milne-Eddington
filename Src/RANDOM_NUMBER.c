
#include "RANDOM_NUMBER.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        06 Mar. 2026  (Hao Li)
          --- Updates:  
              Suport MPI and OpenMP. 
              Redesign the seed generator. 
              Random numbers are generated using GLS if available; 
              otherwise, the methods described in Numerical Recipes are 
              used. 

        30 Oct. 2022  (Hao Li)
          --- Initial commit 
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

void RNG_INIT(STRUCT_RNGState *State, int rank, int thread_id){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Generate seeds for the random number generator.
      Record of revisions:
        06 Mar. 2026.
      Input parameters:
        State, a structure storing the seeds.
        rank, the mpi rank
        thread_id, the openmpi id.
      Return:
        The randum number.
      Reference:
        suggestions from Chatgpt. 
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    int64_t seed = (int64_t)time(NULL);
    seed ^= ((int64_t)(rank + 1)) * 2654435761U;
    seed ^= ((int64_t)(thread_id + 1)) * 6364136223846793005ULL;

#ifdef USE_GSL
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    State->gsl = gsl_rng_alloc(T);
    gsl_rng_set(State->gsl, (unsigned long)(seed & 0xFFFFFFFFUL));
#else
    State->idum = - (int32_t)(seed & 0x7FFFFFFF);
    State->iy = 0;
    State->iset = true;
#endif

    return;
}

/*--------------------------------------------------------------------------------*/

double RNG_UNIFORM(STRUCT_RNGState *State){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Generate a uniformly distributed random number in the interval 
            [0, 1) using either GSL or the method described 
            in Numerical Recipes.
      Record of revisions:
        06 Mar. 2026.
      Input parameters:
        State, a structure storing the seeds.
      Return:
        The randum number.
      Reference:
        Numerical Recipes in C 2ed edition. 
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

#ifdef USE_GSL
    return gsl_rng_uniform(State->gsl); 
#else

    static const int IA = 16807;
    static const int IM = 2147483647;
    //static const double AM = (1.0/IM);
    static const double AM = (1.0/2147483647);
    static const int IQ = 127773;
    static const int IR = 2836;
    //static const int NDIV = (1+(IM-1)/NTAB);
    static const int NDIV = (1+(2147483647-1)/NTAB);
    static const double EPS = 1.2e-7;
    //static const double RNMX = (1.0-EPS); 
    static const double RNMX = (1.0-1.2e-7); 

    int j;
    int32_t k;
    double temp;

    if(State->idum<=0 || !State->iy){
      if (-(State->idum) < 1)
        State->idum = 1;
      else
        State->idum = -(State->idum);

      for(j=NTAB+7; j>= 0; j--){
        k = State->idum/IQ;
        State->idum = IA*(State->idum-k*IQ)-IR*k;
        if (State->idum<0)
          State->idum += IM;
        if (j<NTAB)
          State->iv[j] = State->idum;
      }
      State->iy = State->iv[0];
    }

    k = State->idum/IQ;
    State->idum = IA*(State->idum-k*IQ)-IR*k;
    if(State->idum<0)
      State->idum += IM;

    j = State->iy/NDIV;
    State->iy = State->iv[j];
    State->iv[j] = State->idum;

    temp = AM*State->iy;
    return (temp>RNMX) ? RNMX : temp;
#endif
}

/*--------------------------------------------------------------------------------*/

double RNG_GAUSS(STRUCT_RNGState *State){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Generate a Gaussian-distributed random number with unit width using 
            either GSL or the method described in Numerical Recipes.
      Record of revisions:
        06 Mar. 2026.
      Input parameters:
        State, a structure storing the seeds.
      Return:
        The randum number.
      Reference:
        Numerical Recipes in C 2nd edition / GSL. 
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

#ifdef USE_GSL
    return gsl_ran_gaussian(State->gsl, 1.0); 
#else
    double v1,v2,rsq,fac;
    if(State->iset){
      do{
        v1 = 2.0*RNG_UNIFORM(State)-1.0;
        v2 = 2.0*RNG_UNIFORM(State)-1.0;
        rsq = v1*v1+v2*v2;
      }while(rsq>= 1.0 || rsq==0.0);
      fac = sqrt(-2.0*log(rsq)/rsq);
      State->gset = v1*fac;
      State->iset = false;
      return v2*fac;
    }else{
      State->iset = true;
      return State->gset;
    }
#endif
}

/*--------------------------------------------------------------------------------*/
