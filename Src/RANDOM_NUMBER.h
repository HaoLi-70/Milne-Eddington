
#ifndef RANDOM_NUMBER_h
#define RANDOM_NUMBER_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*--------------------------------------------------------------------------------*/

extern double Ran1(long *idum);

extern double GASDEV(long *idum);

extern int Random_Seed(long *idum);

/*--------------------------------------------------------------------------------*/

#endif /* RANDOM_NUMBER_h */
