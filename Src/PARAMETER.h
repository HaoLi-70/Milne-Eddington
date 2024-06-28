
#ifndef PARAMETER_h
#define PARAMETER_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>

/*--------------------------------------------------------------------------------*/

#define M_SQR(a) ((a)*(a))
#define M_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define M_MAX(x,y) ({      \
typeof(x) _x=(x);    \
typeof(y) _y=(y);    \
(void) (&_x==&_y);   \
_x>_y?_x:_y;         \
})

#define M_MIN(x,y) ({      \
typeof(x) _x=(x);    \
typeof(y) _y=(y);    \
(void) (&_x==&_y);   \
_x>_y?_y:_x;         \
})

#define M_SWAP(x,y) ({     \
(void) (&x==&y);     \
typeof(x) _x;        \
_x=(x);              \
x=(y);               \
y=_x;                \
})

/*--------------------------------------------------------------------------------*/

// Pi Ratio of circumference to diameter          
#define Par_Pi        3.14159265358979323846264338327950288419716939937510582

// square root of Pi
#define Par_SqrtPi    1.77245385090551588191942755656782537698745727539062500
                   
// the reciprocal of Par_SqrtPi                    
#define Par_OneSqrtPi 0.56418958354775627928034964497783221304416656494140625


// Light speed[10^9 m s^-1]
//#define Par_C 2.99792458e-1
#define Par_C 299792458.0

//!> [10^5 cm^-1] to [eV]
#define Par_fktoev 1.2398419739e1

// !> Boltzmann Constant [J K^-1]
#define Par_kb 1.38064852e-23

//!> [10^5 cm^-1] to [J]
#define Par_fktoJ 1.986445824e-18



#define dopp 4.301428e-7

/*--------------------------------------------------------------------------------*/

#endif /* PARAMETER_h */
