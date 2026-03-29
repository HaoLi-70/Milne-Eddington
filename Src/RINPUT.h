
#pragma once

#include <stdio.h>
#include <stdbool.h>
#include <fitsio.h>

#include "ALLOCATION.h"
#include "STR.h"
#include "LOG_ERROR.h"
#include "MPI_INIT.h"

/*--------------------------------------------------------------------------------*/

#define NUM_LINES 20
#define KEY_REQ(k)  {k,"",false,true}
#define KEY_DEF(k,v){k,v,false,false}

/*--------------------------------------------------------------------------------*/

typedef enum Key_Index{
    KEY_LINES,
    KEY_VERBOSE,
    KEY_NUM_ITER,
    KEY_NUM_RUN,
    KEY_CONVG_CRITERIA,
    KEY_CHISQ_CRITERIA,
    KEY_DATA_PATH,
    KEY_WAV_PATH,
    KEY_WEIGHTS,
    KEY_SOL_BOX,
    KEY_OUTPUT_FIT,
    KEY_CACHE_PATH,
    KEY_RESULT_PATH,
    KEY_FIT_PATH,
    KEY_LM_DAMP,
    KEY_LM_ACCEPT,
    KEY_LM_REJECT,
    KEY_LM_LIMITS,
    KEY_BMOD_LIMIT,
    KEY_BTHETA_LIMIT,
    KEY_BPHI_LIMIT,
    KEY_VLOS_LIMIT,
    KEY_DOPPLER_LIMIT,
    KEY_DAMPING_LIMIT,
    KEY_ETA_LIMIT,
    KEY_CONT_LIMIT,
    KEY_BETA_LIMIT,
    KEY_FIXED_BMOD,
    KEY_FIXED_BTHETA,
    KEY_FIXED_BPHI,
    KEY_FIXED_VLOS,
    KEY_FIXED_DOPPLER,
    KEY_FIXED_DAMPING,
    KEY_FIXED_ETA,
    KEY_FIXED_CONT,
    KEY_FIXED_BETA,
    KEY_REGUL_BMOD,
    KEY_REGUL_BTHETA,
    KEY_REGUL_BPHI,
    KEY_REGUL_VLOS,
    KEY_REGUL_DOPPLER,
    KEY_REGUL_DAMPING,
    KEY_REGUL_ETA,
    KEY_REGUL_CONT,
    KEY_REGUL_BETA,
    KEY_VCOEFFI,
    KEY_LCOEFFI,
    KEY_ICRITERIA,
    KEY_SVD_THRESHOLD,
    KEY_HMIREF,
    KEY_FASTMODE,
    KEY_NTHREADS,
    KEY_NFIGURES,
    KEY_STEP,
    KEY_NPROF,
    KEY_TOTAL
}KEY_INDEX;

/*--------------------------------------------------------------------------------*/

typedef struct{

    // the keyword.
    char keyword[Key_Length];
    // the input parameter for the keyword.
    char line[Max_Line_Length];

    // the keyword status.
    bool set; 
    bool required;

}STRUCT_KEYS;

/*--------------------------------------------------------------------------------*/

typedef struct __attribute__((packed)){

    // header of cache files.
    char magic[4];
    int nx, ny, ncache;

}STRUCT_CACHE;

/*--------------------------------------------------------------------------------*/

typedef struct {

    int coord[2];
    int nProf, pcounts; 

}STRUCT_SUBSET;

/*--------------------------------------------------------------------------------*/

typedef struct{

    // verbose level, max number of the iteration, max number for 
    // redo the inversion.
    int verboselv, niter, nruns;

    // number of lines
    int nline;
    // wavelength, effective G factor
    double *Lambda0, *Geff;
    // the weights for Stokes parameters
    double Weights[4];
    // accuraty for voigt function
    int nfigures;

    // thereshold for the SVD decomposition.
    double Threshold;
    // criteria for the convergence and redoing the inversion.
    double Convg_Criteria, Chisq_Criteria;
    // normalizing facor and limits for each parameter.
    double step[9];   
    double Limits[9][2];
    // the damping factor, and values to update the facrot.
    double Damp, Lam_reject, Lam_accept;
    // the limits of the damping factor.
    double Lam_Lim[2];

    // values and weights for regulations
    double value_const[9], Regul_weight[9];
    // invertion flags for each parameter.
    bool inv[9], regl[9], Regul, Subset_flg;

    // output the fit to the profile or not.
    bool output_fit;

    // structure with some precomputed coefficients.
    //STRUCT_FADDEEVA *Fadd;

    bool HMI_REF, fastmode;
    double VCoeffi, LCoeffi, Icriteria;

    // the number of dimension, the number of pixels, x and y sizes 
    // of the data file.
    int naxis, counts, nx, ny;
    // pixel number for each reading
    int nProf;
    // the solution region.
    int sol_box[2][2];
    // type for the observation data type (int, float or double)
    DATA_TYPE type_data;
    // buffer for profile reading.
    double *profbuff, *resbuff, *fitbuff;

    // cache file exists or not.
    bool lcache;
    // cache file header
    STRUCT_CACHE cache_header;
    int *cache;

    // Paths 
    char Data_Path[Max_Line_Length]; 
    char Wav_Path[Max_Line_Length];
    char Cache_Path[Max_Line_Length]; 
    char Log_Path[Max_Line_Length];
    char Result_Path[Max_Line_Length]; 
    char Fit_Path[Max_Line_Length];

    // buffer for inputs
    char tmp_lines[NUM_LINES][Max_Line_Length];

}STRUCT_INPUT;

/*--------------------------------------------------------------------------------*/

extern int RDINPUT(const char Filename[], STRUCT_INPUT *Input, STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/
