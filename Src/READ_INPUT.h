
#ifndef READ_INPUT_h
#define READ_INPUT_h

#include <stdio.h>
#include <stdbool.h>
#include <fitsio.h>
#include "ALLOCATION.h"
#include "STR.h"
#include "ERROR.h"
#include "MPI_CTRL.h"
#include "FADDEEVA.h"

/*--------------------------------------------------------------------------------*/

enum keywordtype {KEYWORD_REQUIRED, KEYWORD_DEFAULT, KEYWORD_OPTIONAL};

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Keywords{

    // the keyword.
    char keyword[Key_Length];

    // the input parameter for the keyword.
    char line[Max_Line_Length];

    // the keyword status.
    bool Set, Required;

}STRUCT_KEYS;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_MEline{

    // center of the wavelength, effective lande factor, a precomputed 
    // coeffcient related to the Zeeman splitting.
    double Lambda0, Geff, BShift;

    // Milne-Eddington model parameters Par[1:9]
    double *Par;

}STRUCT_MELINE;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_input{

    // number of lines
    int nline;

    // structure with the line information.
    STRUCT_MELINE *Lines;

    // line configuration 
    char tmp_lines[20][Max_Line_Length];

    // verbose level, max number of the iteration, max number for 
    // redo the inversion.
    int verboselv, niter, nrun;

    // criteria for the convergence and redoing the inversion.
    double Convg_Criteria, Chisq_Criteria;

    // damping factor
    //double Damp;

    // the number of dimension, the number of pixels, x and y sizes 
    // of the data file.
    int naxis, counts, nx, ny;

    // the solution region.
    int sol_box[2][2];

    // pointer to the fits file.
    fitsfile *fptr_data, *fptr_wav, *fptr_res, *fptr_fit;

    // pointer to the cache file.
    FILE *cache;

    // the first pixel for reading the fits files. 
    long fpixel[4];

    // cache file exists or not.
    bool lcache;

    // size of the chache matrix.
    int nxcache, nycache;

    // output the fit to the profile or not.
    bool output_fit;

    // number of free parameter.
    int npar;

    // invertion flags for each parameter.
    bool inv[10], regl[10];
    double value_const[10], Regul_weight[10];

    // Path to the data, cache, verbose, result files
    char Data_Path[Max_Line_Length], Wav_Path[Max_Line_Length], \
        Cache_Path[Max_Line_Length], Verbose_Path[Max_Line_Length], \
        Result_Path[Max_Line_Length], Fit_Path[Max_Line_Length];

    // normalizing facor (hard coded)
    double step[10];   

    // matrices for the SVD decomposition
    double *W, **V;

    // thereshold for the SVD decomposition.
    double Threshold;

    // type for the input data (int, float or double), wavelength 
    // (float or double), result file ()
    enum data_type type_data, type_wav, type_res;

    int *data_int;
    float *data_flt;

    // structure with some precomputed coefficients.
    STRUCT_FADDEEVA *Fadd;

    bool Regul, Azimuth_Rotate, fastmode;

    double VCoeffi, LCoeffi, Icriteria;

    double delta_v;

    int iw;

}STRUCT_INPUT;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Stokes{

    // number of wavelength points;
    int nl;
    // nomalization factor, the weights for Stokes parameters and 
    // the squre of the weights
    double norm, Weights[4], Weights_SQR[4];
    // the wavelength, input profiles, synthesized profiles, best the fit, noise.
    double *Lambda, **prof, **syn, **fit, *noise;
    // the Jacobian used for the inversion.
    double ***Jacobian;
    // input profiles
    float ***profall;

}STRUCT_STK;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Parameter{

  // 9 model parameters, their best in the inversion, and the guess
  // Bmod, ThetaB, PhiB, V, Dopp, Damp, Eta, Src, Beta;
  double *Par, *Par_Best, *Par_tmp, Par_Guess[10];

  // limits on the parameter
  double Limits[10][2];

  // iteration number;
  int Niter;

  // chisq and best chisq
  double Chisq, Chisq_Best;

}STRUCT_PAR;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Levenberg_Marquardt{

    // the indx array for LU decomposition.
    int *indx;

    // the Hessian matrix, the b vacror and the solutions
    double **Hessian, **Hessian_new, *Jacfvec, *Sol;

    double **Regul_H, *Regul_J;

    // the damping factor, and values to update the facrot.
    double Damp, Lam_reject, Lam_accept, Penalty;

    // the limits of the damping factor.
    double Lam_Lim[2];

}STRUCT_LM;

/*--------------------------------------------------------------------------------*/

extern int RDINPUT(char Filename[], STRUCT_INPUT *Input,  \
    STRUCT_STK *Stk, STRUCT_LM *LM, STRUCT_PAR *Par, STRUCT_MPI *Mpi);

extern int Free_Ram(STRUCT_LM *LM, STRUCT_STK *Stk, STRUCT_INPUT *Input, 
    STRUCT_PAR *Par, STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/

#endif /* READ_INPUT_h */
